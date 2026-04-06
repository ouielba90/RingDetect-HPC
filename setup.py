import os
import subprocess
from setuptools import setup, find_packages, Extension
from setuptools.command.build_py import build_py

class BuildWithMake(build_py):
    """
    Custom build command to compile the Fortran engine locally.
    In CI/CD environments, cibuildwheel handles compilation prior to this step.
    """
    def run(self):
        if not os.path.exists('ringdetect/libringengine.so'):
            try:
                print("Compiling Fortran Shared Library via Makefile...")
                subprocess.check_call(['make', 'lib'])
            except subprocess.CalledProcessError as e:
                print(f"Warning: Local make failed or is unavailable: {e}")
        super().run()

setup(
    name='ringdetect-hpc',
    version='1.2.0', 
    packages=find_packages(),
    # Dummy C extension to force platform-specific wheel generation (e.g., manylinux, macosx)
    ext_modules=[Extension("ringdetect._dummy", sources=["dummy.c"])],
    include_package_data=True,
    package_data={'ringdetect': ['*.so', '*.dll', '*.dylib']},
    install_requires=['numpy'],
    cmdclass={'build_py': BuildWithMake},
)
