
# First steps creating an R package ---------------------------------------

# I followed
# https://tinyheero.github.io/
# jekyll/update/2015/07/26/making-your-first-R-package.html

# Install packages
  install.packages(c("devtools", "roxygen2", "testthat", "knitr", "usethis", "git2r"))
  
# Create Framework (NAMESPACE, DESCRIPTION, R, myfirstpackage.Rproj)
  devtools::create("myfirstpackage")
  #creates automatically a new project in R

# Create a git folder
  #use_git() # Get an error somehow?

# Change DESCRIPTION text by hand
  
# Change license
  use_mit_license("Philipp Kronenberg")

# Add dependencies to package (Defaults to imports)
  usethis::use_package("MASS")
  usethis::use_package("stats")
  usethis::use_package("utils")
  usethis::use_package("ggplot2")
  usethis::use_package("openxlsx")
  usethis::use_package("reshape2")
  usethis::use_package("MCMCpack")
  usethis::use_package("Matrix")
  usethis::use_package("dplyr")
  usethis::use_package("zoo")
  usethis::use_package("tidyr")
  usethis::use_package("grDevices")

# Notes about dependencies:
  # Refer to functions with dplyr::fun()
  # Or import the package/functions for the whole file.
  # Use requireNamespace("data.table", quietly = TRUE) to test if package is
  # installed, then use data.table::fun() to refer to functions.

# Add function description
  # press: Ctrl + Shift + Alt + R for basic roxygen sceleton including all parameters.

# Update function description: Add .Rd files for all functions to man folder
  #install.packages("ps")
  devtools::document()
  # Each time you add new documentation to your R function,
  # you need to run devtools::document() again to re-generate the .Rd files.

# Create some sample data object and save it 
# in data folder and inst/extdata to provide data code, e.g.:
  # x <- c(1:10)
  # usethis::use_data(x)
  # How to call this data:
  #system.file("extdata", "x.R", package = "myfirstpackage")

# Create basic read me file for git
  use_readme_rmd()

# set up data-raw folder
  usethis::use_data_raw()

# Create tests
  use_testthat() # Creates test folder
  use_test("get_IC") # creates test file for single function
  
  
# Create vignette rmarkdown file
  # installing/loading the package:
  if(!require(installr)) { install.packages("installr"); require(installr)} 
  # Installing pandoc
  install.pandoc()
  rmarkdown::find_pandoc()
  install.packages("jsonlite")
  install.packages("rmarkdown")
  # Create new vignette (overrides old vignette file):
  usethis::use_vignette("introduction")

# Install your package
  devtools::install() # or
  devtools::load_all() # replaces installation, and attachement of package
  
  # remove.packages(vctrs)

# Install a package from github:
  # devtools::install_github("yourusername/packagename")

# To use when revising the package ---------------------------------------

# load necessary packages
  library(devtools)
  packageVersion("devtools")
  library(roxygen2)
  library(testthat)
  devtools::session_info()
  #library(tidyverse)
  library(fs)

# update documentation if changes made
  devtools::document()

# load package
  load_all()

# attach package
  library("BayesianDFM")
  
# rebuild vignette
  devtools::build_vignettes() # local I get an error because of pandoc
  
# make basic check whether the package works (including test() and run_examples())
  check() 
  # runs examples if problems occur in check()
  devtools::run_examples()
  # Check documentation
  devtools::check_man()
  # Test your functions
  test()

# Create package including Create/update Vignette
  devtools::build() # local I get an error because of pandoc

  # Render your read me file
  build_readme()


  # install the package from git
  install.packages("devtools")
  devtools::install_github("philippkronenberg/BayesianDFM")


  # create correct rda files of output data
  #use_data(dat_final,data,inventory,metadata)

  # Create package website with pkgdown():
    # Run once to configure your package to use pkgdown
    usethis::use_pkgdown()
  # Use pkgdown to update your website:

    # Run to build the website
    pkgdown::build_site()

# Open Issues -------------------------------------------------------------

# Warnings:
  #> checking data for ASCII and uncompressed saves ... OK
  #WARNING
  #‘qpdf’ is needed for checks on size reduction of PDFs

# Solution:
  # Install Rtools (is not available as a package for R 3.6, but can be
  # downloaded from website (Rtools40))

  # Questions:
  # How do I install it in the package(since it is not a regualr package)?

# Next Steps --------------------------------------------------------------

  #Done
  # Create basic structure of package
  # Create function
  # Add descriptions to functions
  # Add dependencies
  # Add sample data
  # Create read me content for git
  # Add correct contact information into description
  # Add at all functions in the functions the package name and :: or importFrom
  # nameofpackage nameoffunction
  # Creating Vignette content
  # Correct read in of data files

  #ToDo

# Add automated testing (github actions)
# Create package website with pkgdown()


