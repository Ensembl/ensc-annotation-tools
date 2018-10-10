# EnsC-annotation-tools
This repository contains the OBDA and Programs directory from ensc-core found at https://github.com/Ensembl/ensc-core
The binary used by the Ensembl Genome Annotation System are:
* translate
* fastasplit_random
* indicate
* RefineSolexaGenes

## Dependencies
* EnsC-core
* MySQL (recommended)
* samtools (recommended)
* libconfig (recommended)

## To compile the binaries
```
cd src
./autogen.sh
./runconfig
make
```

You can also specify parameters to runconfig like `CPPFLAGS="-I/opt/local/include" LDFLAGS="-L/opt/local/lib"`

## Previous history
All unused code have been removed, so if you want to see the full history for some of the files you need to use the code from ensc-core
