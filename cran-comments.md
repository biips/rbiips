## Test environments
* local ubuntu 16.04 64-bit, R 3.3.2
* local OS X 10.8, R 3.1.2
* local windows 7 64-bit, R 3.3.2
* win-builder (devel)

## R CMD check results
Status: 2 WARNINGs, 2 NOTEs

* checking compilation flags in Makevars ... WARNING
Non-portable flags in variable 'PKG_CXXFLAGS':
  -Wno-deprecated-declarations
  
-> To avoid warnings due to boost using deprecated 
std::binary_function in Accumulators library
cf. http://lists.boost.org/Archives/boost/2016/05/229402.php

* checking compiled code ... WARNING
File ‘rbiips/libs/rbiips.so’:
  Found ‘_ZSt4cout’, possibly from ‘std::cout’ (C++)
    Object: ‘biips/src/core/graph/Graph.o’
  Found ‘exit’, possibly from ‘exit’ (C)
    Object: ‘biips/src/compiler/compiler/scanner.o’
  Found ‘stderr’, possibly from ‘stderr’ (C)
    Objects: ‘biips/src/compiler/compiler/parser.o’,
      ‘biips/src/compiler/compiler/scanner.o’
  Found ‘stdout’, possibly from ‘stdout’ (C)
    Object: ‘biips/src/compiler/compiler/scanner.o’

Compiled code should not call entry points which might terminate R nor
write to stdout/stderr instead of to the console, nor the system RNG.

See ‘Writing portable packages’ in the ‘Writing R Extensions’ manual.

-> These calls come from the external biips library.
They should not be called.
scanner.o and parser.o are automatically generated by flex/bison.

* checking installed package size ... NOTE
  installed size is 76.2Mb
  sub-directories of 1Mb or more:
    libs  75.1Mb
    
-> The package embeds the complete sources of biips c++ library 
in src/biips (only 17Mb when compiled with -O3)

* checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements.

-> I use $(wildcard) GNU make extension because the list
of source files is too long

## Downstream dependencies
There are currently no downstream dependencies for this package.