/**
 * @file main.c
 * @brief HPC Ring Detector CLI & C-API
 * High-performance cycle perception engine parser. Reads common molecular
 * formats and streams multi-frame MD trajectories to Text, JSON, or CSV.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>

// Fortran Engine Signature
extern void find_rings(int *n_atoms, double *x, double *y, double *z, double *radii, 
                       int *max_ring, char sep, int *threads, int *target_rings, 
                       char *out_filename, double *cell, int *active_mask);

typedef struct {
    const char *symbol;
    double radius;
} ElementRadius;

double get_covalent_radius(const char *element) {
    // Comprehensive Periodic Table Covalent Radii (in Angstroms)
    static const ElementRadius pt[] = {
        {"H", 0.31}, {"He", 0.28}, {"Li", 1.28}, {"Be", 0.96},
        {"B", 0.84}, {"C", 0.76}, {"N", 0.71}, {"O", 0.66},
        {"F", 0.57}, {"Ne", 0.58}, {"Na", 1.66}, {"Mg", 1.41},
        {"Al", 1.21}, {"Si", 1.11}, {"P", 1.07}, {"S", 1.05},
        {"Cl", 1.02}, {"Ar", 1.06}, {"K", 2.03}, {"Ca", 1.76},
        {"Sc", 1.70}, {"Ti", 1.60}, {"V", 1.53}, {"Cr", 1.39},
        {"Mn", 1.39}, {"Fe", 1.32}, {"Co", 1.26}, {"Ni", 1.24},
        {"Cu", 1.32}, {"Zn", 1.22}, {"Ga", 1.22}, {"Ge", 1.20},
        {"As", 1.19}, {"Se", 1.20}, {"Br", 1.20}, {"Kr", 1.16},
        {"Rb", 2.20}, {"Sr", 1.95}, {"Y", 1.90}, {"Zr", 1.75},
        {"Nb", 1.64}, {"Mo", 1.54}, {"Tc", 1.47}, {"Ru", 1.46},
        {"Rh", 1.42}, {"Pd", 1.39}, {"Ag", 1.45}, {"Cd", 1.44},
        {"In", 1.42}, {"Sn", 1.39}, {"Sb", 1.39}, {"Te", 1.38},
        {"I", 1.39}, {"Xe", 1.40}, {"Cs", 2.44}, {"Ba", 2.15},
        {"La", 2.07}, {"Ce", 2.04}, {"Pr", 2.03}, {"Nd", 2.01},
        {"Pm", 1.99}, {"Sm", 1.98}, {"Eu", 1.98}, {"Gd", 1.96},
        {"Tb", 1.94}, {"Dy", 1.92}, {"Ho", 1.92}, {"Er", 1.89},
        {"Tm", 1.90}, {"Yb", 1.87}, {"Lu", 1.87}, {"Hf", 1.75},
        {"Ta", 1.70}, {"W", 1.62}, {"Re", 1.51}, {"Os", 1.44},
        {"Ir", 1.41}, {"Pt", 1.36}, {"Au", 1.36}, {"Hg", 1.32},
        {"Tl", 1.45}, {"Pb", 1.46}, {"Bi", 1.48}, {"Po", 1.40},
        {"At", 1.50}, {"Rn", 1.50}, {"Fr", 2.60}, {"Ra", 2.21},
        {"Ac", 2.15}, {"Th", 2.06}, {"Pa", 2.00}, {"U", 1.96},
        {"Np", 1.90}, {"Pu", 1.87}, {"Am", 1.80}, {"Cm", 1.69}
    };
    int num_elements = sizeof(pt) / sizeof(pt[0]);

    char clean_elem[3] = {0};
    if (element[0]) clean_elem[0] = toupper((unsigned char)element[0]);
    if (element[1]) clean_elem[1] = tolower((unsigned char)element[1]);

    for (int i = 0; i < num_elements; i++) {
        if (strcmp(clean_elem, pt[i].symbol) == 0) return pt[i].radius;
    }
    return 1.00; 
}

void print_help() {
    printf("==============================================================\n");
    printf(" HPC Ring Detector (Cycle Perception Engine)\n");
    printf("==============================================================\n");
    printf("Usage: ./ring_detector <molecule_file> [OPTIONS]\n\n");
    printf("Options:\n");
    printf("  -h, --help        Show this help message and exit.\n");
    printf("  -f <format>       Set the input file format (default: xyz).\n");
    printf("  -a <mask_list>    Active atoms mask (e.g., 1-15,30).\n");
    printf("  -c <x> <y> <z>    Set unit cell dimensions for Periodic Boundaries.\n");
    printf("  -m <number>       Set maximum ring depth to search (default: 6).\n");
    printf("  -r <list>         Search ONLY for specific ring sizes (comma-separated).\n");
    printf("  -p <threads>      Set OpenMP threads (default: max available).\n");
    printf("  -s <char>         Set output separator character (default: ' ').\n");
    printf("  -j                Output results in strict JSON format.\n");
    printf("  -v                Output results in flat CSV format (Big Data optimized).\n\n");
    printf("==============================================================\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) { print_help(); return 1; }
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(); return 0;
        }
    }

    // CLI Defaults
    char *filename = argv[1];
    char format_str[10] = "xyz"; 
    int max_ring = 6;      
    char sep = '-'; // Defaulting to '-' to ensure safe CSV parsing       
    int threads = 0;      
    int target_rings[100] = {0}; 
    int use_targets = 0;
    int json_output = 0; 
    int csv_output = 0;
    double cell[3] = {0.0, 0.0, 0.0}; 
    char active_str[512] = "";

    // Parse Arguments
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            strncpy(format_str, argv[i+1], sizeof(format_str)-1); i++;
        } else if (strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            max_ring = atoi(argv[i+1]); i++;
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            sep = argv[i+1][0]; i++;
        } else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc) {
            threads = atoi(argv[i+1]); i++;
        } else if (strcmp(argv[i], "-j") == 0) {
            json_output = 1; 
            csv_output = 0; // mutually exclusive
        } else if (strcmp(argv[i], "-v") == 0) {
            csv_output = 1;
            json_output = 0; // mutually exclusive
        } else if (strcmp(argv[i], "-c") == 0 && i + 3 < argc) {
            cell[0] = atof(argv[i+1]); cell[1] = atof(argv[i+2]); cell[2] = atof(argv[i+3]); i += 3;
        } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
            use_targets = 1;
            char rings_str[256];
            strncpy(rings_str, argv[i+1], sizeof(rings_str)-1);
            char *token = strtok(rings_str, ",");
            max_ring = 0; 
            while (token != NULL) {
                int r = atoi(token);
                if (r >= 3 && r <= 100) { target_rings[r - 1] = 1; if (r > max_ring) max_ring = r; }
                token = strtok(NULL, ",");
            }
            i++;
        } else if (strcmp(argv[i], "-a") == 0 && i + 1 < argc) {
            strncpy(active_str, argv[i+1], sizeof(active_str)-1); i++;
        }
    }

    if (max_ring > 100) max_ring = 100;
    if (!use_targets) { for (int i = 3; i <= max_ring; i++) target_rings[i - 1] = 1; }

    // Filenames
    char base_name[256], out_filename[256], temp_filename[256];
    char *slash = strrchr(filename, '/'); char *backslash = strrchr(filename, '\\');
    char *file_start = filename;
    if (slash != NULL) file_start = slash + 1;
    if (backslash != NULL && backslash > slash) file_start = backslash + 1;
    strncpy(base_name, file_start, sizeof(base_name) - 1);
    char *dot = strrchr(base_name, '.'); if (dot != NULL) *dot = '\0';
    
    if (json_output) {
        snprintf(out_filename, sizeof(out_filename), "%s.json", base_name);
        snprintf(temp_filename, sizeof(temp_filename), "%s.rings.tmp", base_name);
    } else if (csv_output) {
        snprintf(out_filename, sizeof(out_filename), "%s.csv", base_name);
        snprintf(temp_filename, sizeof(temp_filename), "%s.rings.tmp", base_name);
    } else {
        snprintf(out_filename, sizeof(out_filename), "%s.rings", base_name);
        snprintf(temp_filename, sizeof(temp_filename), "%s", out_filename); // Text writes directly
    }

    FILE *file;
    if (strcmp(filename, "-") == 0) {
        file = stdin;
    } else {
        file = fopen(filename, "r");
        if (!file) { printf("Error: Cannot open file '%s'\n", filename); return 1; }
    }
    
    FILE *out = fopen(out_filename, "w");
    if (json_output) {
        fprintf(out, "{\n  \"molecule\": \"%s\",\n  \"frames\": [\n", base_name);
    } else if (csv_output) {
        fprintf(out, "Frame,RingSize,Planar,Indices\n");
    }

    char buffer[512];
    int frame = 0;
    double total_compute_time = 0.0;

    printf("C Engine Config -> Format: %s | JSON: %s | CSV: %s | Max Depth: %d | Sep: '%c' | Threads: %d\n", 
           format_str, json_output ? "ON" : "OFF", csv_output ? "ON" : "OFF", max_ring, sep, threads);
    printf("Processing trajectory...\n");

    // --- MAIN TRAJECTORY LOOP ---
    while (1) {
        int N = 0;
        
        if (strcmp(format_str, "xyz") == 0) {
            if (!fgets(buffer, sizeof(buffer), file)) break; // EOF
            if (sscanf(buffer, "%d", &N) != 1) break; 
            fgets(buffer, sizeof(buffer), file); // skip comment line
        } else {
            if (frame > 0) break; // Multi-frame only natively handles XYZ
            if (strcmp(format_str, "mol") == 0) {
                fgets(buffer, sizeof(buffer), file); fgets(buffer, sizeof(buffer), file); fgets(buffer, sizeof(buffer), file); 
                if (!fgets(buffer, sizeof(buffer), file) || sscanf(buffer, "%d", &N) != 1) break;
                rewind(file); 
            } else if (strcmp(format_str, "pdb") == 0) {
                while (fgets(buffer, sizeof(buffer), file)) {
                    if (strncmp(buffer, "ATOM  ", 6) == 0 || strncmp(buffer, "HETATM", 6) == 0) N++;
                }
                rewind(file);
            } else {
                while (fgets(buffer, sizeof(buffer), file)) { if (strlen(buffer) > 2) N++; }
                rewind(file); 
            }
        }

        if (N == 0) break;
        frame++;

        if (frame % 100 == 0) {
            printf("\r>> Processing Frame: %d...", frame);
            fflush(stdout);
        }

        // 2. Allocate Arrays
        double *x = (double *)malloc(N * sizeof(double));
        double *y = (double *)malloc(N * sizeof(double));
        double *z = (double *)malloc(N * sizeof(double));
        double *radii = (double *)malloc(N * sizeof(double));
        int *active_mask = (int *)malloc(N * sizeof(int));
        
        int atom_idx = 0; char element[20];
        if (strcmp(format_str, "mol") == 0) { for(int k=0; k<4; k++) fgets(buffer, sizeof(buffer), file); }

        // 3. Read Coordinates
        while (atom_idx < N && fgets(buffer, sizeof(buffer), file)) {
            int len = strlen(buffer); if (len < 3) continue; 
            int success = 0;
            if (strcmp(format_str, "xyz") == 0 || strcmp(format_str, "raw") == 0) {
                success = (sscanf(buffer, "%19s %lf %lf %lf", element, &x[atom_idx], &y[atom_idx], &z[atom_idx]) == 4);
            } else if (strcmp(format_str, "csv") == 0) {
                success = (sscanf(buffer, " %19[^, \t] , %lf , %lf , %lf", element, &x[atom_idx], &y[atom_idx], &z[atom_idx]) == 4);
            } else if (strcmp(format_str, "pdb") == 0) {
                if (strncmp(buffer, "ATOM  ", 6) == 0 || strncmp(buffer, "HETATM", 6) == 0) {
                    if (len >= 54) { 
                        char x_str[9] = {0}, y_str[9] = {0}, z_str[9] = {0}, elem_str[3] = {0};
                        strncpy(x_str, buffer + 30, 8); strncpy(y_str, buffer + 38, 8); strncpy(z_str, buffer + 46, 8);
                        if (len >= 78) strncpy(elem_str, buffer + 76, 2);
                        if (elem_str[0] == ' ' || elem_str[0] == '\0') strncpy(elem_str, buffer + 12, 2); 
                        int k = 0; for (int j = 0; j < 2; j++) { if (isalpha((unsigned char)elem_str[j])) element[k++] = elem_str[j]; }
                        element[k] = '\0';
                        x[atom_idx] = atof(x_str); y[atom_idx] = atof(y_str); z[atom_idx] = atof(z_str); success = 1;
                    }
                }
            }
            if (success) { radii[atom_idx] = get_covalent_radius(element); atom_idx++; }
        }

        // 4. Parse Mask
        for (int i = 0; i < N; i++) active_mask[i] = 1; 
        if (strlen(active_str) > 0) {
            for (int i = 0; i < N; i++) active_mask[i] = 0;
            char *str_copy = strdup(active_str);
            char *token = strtok(str_copy, ",");
            while (token != NULL) {
                char *dash = strchr(token, '-');
                if (dash != NULL) {
                    *dash = '\0';
                    int start = atoi(token); int end = atoi(dash + 1);
                    for (int idx = start; idx <= end; idx++) { if (idx >= 1 && idx <= N) active_mask[idx - 1] = 1; }
                } else {
                    int idx = atoi(token);
                    if (idx >= 1 && idx <= N) active_mask[idx - 1] = 1;
                }
                token = strtok(NULL, ",");
            }
            free(str_copy);
        }

        // 5. Execute Fortran!
        double start_time = omp_get_wtime();
        
        if (json_output || csv_output) {
            find_rings(&N, x, y, z, radii, &max_ring, sep, &threads, target_rings, temp_filename, cell, active_mask);
        } else {
            fprintf(out, "=== FRAME %d ===\n", frame);
            fclose(out);
            find_rings(&N, x, y, z, radii, &max_ring, sep, &threads, target_rings, out_filename, cell, active_mask);
            out = fopen(out_filename, "a");
        }
        
        total_compute_time += (omp_get_wtime() - start_time);

        // 6. JSON / CSV Post-Processing
        if (json_output || csv_output) {
            FILE *in = fopen(temp_filename, "r");
            
            if (json_output && frame > 1) fprintf(out, ",\n");
            if (json_output) fprintf(out, "    {\n      \"frame\": %d,\n      \"total_atoms\": %d,\n      \"rings\": [\n", frame, N);
            
            char r_line[1024]; int first_item = 1;
            while (fgets(r_line, sizeof(r_line), in)) {
                int depth; char rest[512]; int is_planar = 0;
                if (sscanf(r_line, "%d-MR (PLANAR): %[^\n]", &depth, rest) == 2) is_planar = 1;
                else if (sscanf(r_line, "%d-MR: %[^\n]", &depth, rest) == 2) is_planar = 0;
                else continue; 

                if (json_output) {
                    for(int i = 0; rest[i]; i++) { if (rest[i] == sep) rest[i] = ','; }
                    if (!first_item) fprintf(out, ",\n");
                    fprintf(out, "        {\"size\": %d, \"planar\": %s, \"indices\": [%s]}", depth, is_planar ? "true" : "false", rest);
                } else if (csv_output) {
                    // Replace any accidental commas in 'rest' to avoid breaking CSV format
                    for(int i = 0; rest[i]; i++) { if (rest[i] == ',') rest[i] = '-'; }
                    fprintf(out, "%d,%d,%s,%s\n", frame, depth, is_planar ? "True" : "False", rest);
                }
                
                first_item = 0;
            }
            if (json_output) fprintf(out, "\n      ]\n    }");
            
            fclose(in);
            remove(temp_filename);
        }

        free(x); free(y); free(z); free(radii); free(active_mask);
    }

    if (json_output) fprintf(out, "\n  ]\n}\n");
    fclose(out);
    if (file != stdin) fclose(file);

    printf("\n------------------------------------------------\n");
    printf("Processed %d Frames.\n", frame);
    printf("Output generated: %s\n", out_filename);
    printf("TOTAL COMPUTE TIME: %f seconds\n", total_compute_time);
    printf("------------------------------------------------\n");

    return 0;
}
