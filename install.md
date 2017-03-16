## 1. Installation

The program can be downloaded from GeDI homepage, <http://github.com/BHSAI/GeDI>. There is a Linux 64 bit executable,

    gedi_v1.0/bin/gedi_linux_x86

Compilation requires a C++ compiler (e.g., g++) and GNU Scientific library [(GSL)](http://www.gnu.org/software/gsl/). 
See [here](gsl.md) for a quick recipe of GSL installation.

In most cases, compilation should work with

    $ cd gedi_v1.0/src
    $ cp makefile_srl makefile
    $ make

which will create an executable file named gedi. Please edit makefile appropriately to suite your system and, in particular, the CFLAGS and LDFLAGS option to allow the GSL library to be linked. For instance, if you installed GSL to a local directory not on the linker search path, uncomment the GSLDIR line and modify it to point to that directory, uncomment the two following lines, and comment out the next two lines. To compile a parallel (MPI) version, use the file makefile_mpi:

    $ cp makefile_mpi makefile
    $ make

Note that the mpi compiler and its options may differ depending on your system. Refer to your system parallel compilation documentation for details. The "-DMPIP" option for the CC flag turns on the parallel part of the code.


***
[Up](README.md)

[Next](usage.md)
