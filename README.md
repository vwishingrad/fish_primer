To install the package dependencies, open the project (.Rproj file) in RStudio, and in the RStudio console, do:

renv::restore()

It will install everything needed, but only within the context of the project (i.e., it won't install them over whatever packages you might already have installed, it makes its own environment).
If you encounter build errors, you may need to go into renv.lock and change the version number for whichever packages are failing to their newest version (look on CRAN).
