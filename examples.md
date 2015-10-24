##Examples

An example data set is provided:

     ./example/chr19sim.tfam
     ./example/chr19sim.tped

These are case-control (1000 each, 2000 total) data of 100 SNP segment around rs810188 of _C3_ gene on chr19, generated using [GWASimulator](http://biostat.mc.vanderbilt.edu/wiki/Main/GWAsimulator) and [1000 Genome](http://www.1000genomes.org) haplotypes. The SNP rs8110188 is the causal SNP. The script file

     ./example/run/gedi.sh

illustrates some of the usage based on this data set. Run this with

     $ cd ./example/run/
     
     $ chmod 775 gedi.sh
     
     $ ./gedi.sh

***
[Up](README.md)
