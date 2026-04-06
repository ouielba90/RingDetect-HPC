!> @brief HPC Fortran Engine for Ring Detection
!> Constructs an adjacency graph using covalent radii and executes a highly
!> parallelized Depth-First Search (DFS) to identify molecular cycles.
module ring_engine
    use iso_c_binding
    use omp_lib
    implicit none
    
    integer, parameter :: MAX_ALLOWED_RING = 100 
    integer, parameter :: BUFFER_SIZE = 65536 

contains

    subroutine find_rings(n_atoms, x, y, z, r, max_ring, sep, threads, target_rings, c_out_filename, cell, active_mask) bind(c, name="find_rings")
        integer(c_int), intent(in) :: n_atoms, max_ring, threads
        real(c_double), intent(in) :: x(n_atoms), y(n_atoms), z(n_atoms), r(n_atoms)
        character(c_char), intent(in), value :: sep
        integer(c_int), intent(in) :: target_rings(100) 
        character(kind=c_char), intent(in) :: c_out_filename(*) 
        real(c_double), intent(in) :: cell(3)
        integer(c_int), intent(in) :: active_mask(n_atoms)
        
        integer :: i, j, tid, num_t, io_status, k
        real(c_double) :: dx, dy, dz, dist_sq, rcut_sq
        integer, allocatable :: num_bonds(:), adj_list(:,:)
        integer :: total_rings(MAX_ALLOWED_RING)
        
        character(len=256) :: filename
        character(len=256) :: f_out_filename 
        character(len=BUFFER_SIZE) :: thread_buffer
        integer :: buf_ptr
        integer :: thread_path(MAX_ALLOWED_RING)
        integer :: thread_totals(MAX_ALLOWED_RING)
        integer :: safe_max_ring
        integer :: file_size 
        character(len=:), allocatable :: exact_buffer 
        
        ! Extract C-string to Fortran
        f_out_filename = " "
        k = 1
        do while (c_out_filename(k) /= c_null_char .and. k <= 256)
            f_out_filename(k:k) = c_out_filename(k)
            k = k + 1
        end do
        
        safe_max_ring = max_ring
        if (safe_max_ring > MAX_ALLOWED_RING) safe_max_ring = MAX_ALLOWED_RING
        
        if (threads > 0) call omp_set_num_threads(threads)
        num_t = omp_get_max_threads()
        
        allocate(num_bonds(n_atoms), adj_list(16, n_atoms))
        num_bonds = 0; total_rings = 0

        print *, "Fortran Engine: Building Adjacency Graph..."
        !$omp parallel do private(j, dx, dy, dz, dist_sq, rcut_sq) schedule(dynamic)
        do i = 1, n_atoms
            if (active_mask(i) == 0) cycle
            do j = 1, n_atoms
                if (i == j .or. active_mask(j) == 0) cycle
                
                dx = x(i) - x(j)
                dy = y(i) - y(j)
                dz = z(i) - z(j)
                
                ! Periodic Boundary Conditions (MIC)
                if (cell(1) > 0.0_c_double) dx = dx - cell(1) * anint(dx / cell(1))
                if (cell(2) > 0.0_c_double) dy = dy - cell(2) * anint(dy / cell(2))
                if (cell(3) > 0.0_c_double) dz = dz - cell(3) * anint(dz / cell(3))
                
                dist_sq = dx**2 + dy**2 + dz**2
                rcut_sq = (r(i) + r(j) + 0.45_c_double)**2
                
                if (dist_sq <= rcut_sq .and. num_bonds(i) < 16) then
                    num_bonds(i) = num_bonds(i) + 1
                    adj_list(num_bonds(i), i) = j
                end if
            end do
        end do
        !$omp end parallel do

        print *, "Fortran Engine: Hunting Rings & Checking Planarity..."

        !$omp parallel private(tid, filename, thread_buffer, buf_ptr, i, thread_path, thread_totals)
        tid = omp_get_thread_num()
        write(filename, '("t",I0,".tmp")') tid
        open(unit=100+tid, file=trim(filename), status="replace", access="stream")
        
        thread_buffer = ""
        buf_ptr = 1
        thread_totals = 0

        !$omp do schedule(dynamic, 4)
        do i = 1, n_atoms
            if (active_mask(i) == 0) cycle
            thread_path(1) = i
            call dfs(i, i, 1, thread_path, thread_totals, safe_max_ring, sep, 100+tid, thread_buffer, buf_ptr)
        end do
        !$omp end do
        
        if (buf_ptr > 1) write(100+tid) thread_buffer(1:buf_ptr-1)
        close(100+tid)
        
        !$omp critical
        total_rings = total_rings + thread_totals
        !$omp end critical
        !$omp end parallel

        ! Merge thread files
        open(unit=10, file=trim(f_out_filename), status="replace")
        do i = 0, num_t - 1
            write(filename, '("t",I0,".tmp")') i
            inquire(file=trim(filename), size=file_size)
            
            if (file_size > 0) then
                open(unit=30, file=trim(filename), status="old", access="stream")
                allocate(character(len=file_size) :: exact_buffer)
                read(30) exact_buffer
                write(10, '(A)', advance='no') exact_buffer
                deallocate(exact_buffer)
                close(30, status="delete")
            else
                open(unit=30, file=trim(filename), status="old", iostat=io_status)
                if (io_status == 0) close(30, status="delete")
            end if
        end do
        close(10)

        print *, "--- FINAL TOTALS ---"
        do i = 3, safe_max_ring
            if (total_rings(i) > 0) print '(I0, A, I0)', i, "-MR: ", total_rings(i)
        end do
        deallocate(num_bonds, adj_list)

    contains

        subroutine append_itoa(val, buf, ptr)
            integer, intent(in) :: val
            character(len=BUFFER_SIZE), intent(inout) :: buf
            integer, intent(inout) :: ptr
            integer :: temp, digits, m
            character(len=10) :: str
            
            if (val == 0) then
                buf(ptr:ptr) = '0'
                ptr = ptr + 1
                return
            end if
            
            temp = val
            digits = 0
            do while (temp > 0)
                digits = digits + 1
                temp = temp / 10
            end do
            
            temp = val
            do m = digits, 1, -1
                str(m:m) = char(48 + mod(temp, 10))
                temp = temp / 10
            end do
            
            buf(ptr:ptr+digits-1) = str(1:digits)
            ptr = ptr + digits
        end subroutine append_itoa

        recursive subroutine dfs(start, curr, depth, path, totals, max_r, s, u, buf, ptr)
            integer, intent(in) :: start, curr, depth, max_r, u
            integer, intent(inout) :: path(MAX_ALLOWED_RING)
            integer, intent(inout) :: totals(MAX_ALLOWED_RING), ptr
            character(c_char), intent(in) :: s
            character(len=BUFFER_SIZE), intent(inout) :: buf
            
            integer :: b, neighbor, step, m
            logical :: visited
            real(c_double) :: v1x, v1y, v1z, v2x, v2y, v2z, nx, ny, nz, n_len, px, py, pz, dist
            integer :: is_planar
            
            do b = 1, num_bonds(curr)
                neighbor = adj_list(b, curr)
                
                if (neighbor == start .and. depth >= 3) then
                    if (path(2) < curr) then
                        if (target_rings(depth) == 1) then
                            
                            v1x = x(path(2)) - x(path(1))
                            v1y = y(path(2)) - y(path(1))
                            v1z = z(path(2)) - z(path(1))
                            
                            if (depth == 3) then
                                v2x = x(curr) - x(path(1))
                                v2y = y(curr) - y(path(1))
                                v2z = z(curr) - z(path(1))
                            else
                                v2x = x(path(3)) - x(path(1))
                                v2y = y(path(3)) - y(path(1))
                                v2z = z(path(3)) - z(path(1))
                            end if
                            
                            if (cell(1) > 0.0_c_double) v1x = v1x - cell(1) * anint(v1x / cell(1))
                            if (cell(2) > 0.0_c_double) v1y = v1y - cell(2) * anint(v1y / cell(2))
                            if (cell(3) > 0.0_c_double) v1z = v1z - cell(3) * anint(v1z / cell(3))
                            if (cell(1) > 0.0_c_double) v2x = v2x - cell(1) * anint(v2x / cell(1))
                            if (cell(2) > 0.0_c_double) v2y = v2y - cell(2) * anint(v2y / cell(2))
                            if (cell(3) > 0.0_c_double) v2z = v2z - cell(3) * anint(v2z / cell(3))

                            nx = v1y*v2z - v1z*v2y
                            ny = v1z*v2x - v1x*v2z
                            nz = v1x*v2y - v1y*v2x
                            n_len = sqrt(nx**2 + ny**2 + nz**2)
                            
                            if (n_len > 1e-6) then
                                nx = nx/n_len; ny = ny/n_len; nz = nz/n_len
                            end if
                            
                            is_planar = 1
                            do m = 4, depth
                                if (m == depth) then
                                    px = x(curr); py = y(curr); pz = z(curr)
                                else
                                    px = x(path(m)); py = y(path(m)); pz = z(path(m))
                                end if
                                
                                dx = px - x(path(1))
                                dy = py - y(path(1))
                                dz = pz - z(path(1))
                                
                                if (cell(1) > 0.0_c_double) dx = dx - cell(1) * anint(dx / cell(1))
                                if (cell(2) > 0.0_c_double) dy = dy - cell(2) * anint(dy / cell(2))
                                if (cell(3) > 0.0_c_double) dz = dz - cell(3) * anint(dz / cell(3))
                                
                                dist = abs(dx*nx + dy*ny + dz*nz)
                                if (dist > 0.15_c_double) then 
                                    is_planar = 0
                                    exit
                                end if
                            end do

                            totals(depth) = totals(depth) + 1
                            call append_itoa(depth, buf, ptr)
                            
                            if (is_planar == 1) then
                                buf(ptr:ptr+13) = '-MR (PLANAR): '
                                ptr = ptr + 14
                            else
                                buf(ptr:ptr+4) = '-MR: '
                                ptr = ptr + 5
                            end if
                            
                            do m = 1, depth - 1
                                call append_itoa(path(m), buf, ptr)
                                buf(ptr:ptr) = s
                                ptr = ptr + 1
                            end do
                            call append_itoa(curr, buf, ptr)
                            
                            buf(ptr:ptr) = char(10) 
                            ptr = ptr + 1
                            
                            if (ptr > BUFFER_SIZE - 200) then
                                write(u) buf(1:ptr-1)
                                ptr = 1
                            end if
                        end if
                    end if
                else if (neighbor > start .and. depth < max_r) then
                    visited = .false.
                    do step = 1, depth
                        if (path(step) == neighbor) then
                            visited = .true.
                            exit
                        end if
                    end do
                    
                    if (.not. visited) then
                        path(depth + 1) = neighbor
                        call dfs(start, neighbor, depth + 1, path, totals, max_r, s, u, buf, ptr)
                    end if
                end if
            end do
        end subroutine dfs
    end subroutine find_rings
end module ring_engine
