set(CMAKE_Fortran_COMPILER "/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "14.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_COMPILER_FRONTEND_VARIANT "GNU")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/cvmfs/sft.cern.ch/lcg/releases/binutils/2.40-acaab/x86_64-el9/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/bin/gcc-ar")
set(CMAKE_RANLIB "/cvmfs/sft.cern.ch/lcg/releases/binutils/2.40-acaab/x86_64-el9/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib/gcc/x86_64-pc-linux-gnu/14.2.0/finclude;/cvmfs/sft.cern.ch/lcg/views/LCG_106a/x86_64-el9-gcc14-opt/include;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib/gcc/x86_64-pc-linux-gnu/14.2.0/include-fixed;/usr/local/include;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib/gcc/x86_64-pc-linux-gnu/14.2.0;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib/gcc;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib64;/lib64;/usr/lib64;/cvmfs/sft.cern.ch/lcg/releases/gcc/14.2.0-2f0a0/x86_64-el9/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
