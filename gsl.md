### Gnu Scientific Library installation

Download the latest version from [GSL](http://www.gnu.org/software/gsl/). Untar the archive by

    $ tar xvf gsl-latest.tar
    $ cd gsl-[version#]/

Configure the installation by

    $ ./configure --prefix=/home/username/local

Replace the arguments of the "--prefix" with the target directory for the library where you have the write access. If you omit "--prefix", the installation will default to /usr/local.

Compile and install by

    $ make install

***
[Back](install.md)
