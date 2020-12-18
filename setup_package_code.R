
# First steps creating an R package ---------------------------------------

# I followed https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Install packages
  install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
  install.packages("git2r")
  install.packages("usethis")
# Create Framework (NAMESPACE, DESCRIPTION,R, myfirstpackage.Rproj)
  devtools::create("myfirstpackage")
# creates automatically a new project in R
  
  # Create a git folder
  #use_git() # Get an error somehow?  
  
  # Change DESCRIPTION text
  # Change license
  use_mit_license("Philipp Kronenberg")

# Add dependencies to package # Defaults to imports
  #install.packages(c("reshape2", "MASS", "stats", "ggplot2","openxlsx",
  #                   "MCMCpack","Matrix","dplyr","zoo","tidyr"))
  
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

  #> Adding dplyr to Imports
  #> Refer to functions with dplyr::fun()
  #usethis::use_package("data.table", "Suggests")
  #> Use requireNamespace("data.table", quietly = TRUE) to test if package is 
  #>  installed, then use data.table::fun() to refer to functions.
  
# Add function description
  # press: Ctrl + Shift + Alt + R for basic roxygen sceleton including all parameters.
  
# Update function description: Add .Rd file to man folder
  #install.packages("ps")
  devtools::document()
  # Each time you add new documentation to your R function, 
  # you need to run devtools::document() again to re-generate the .Rd files.
  
# Create some sample data object and save it in data folder and inst/extdata to provide data code
  x <- c(1:10)
  usethis::use_data(x)
  # How to call this data:
  #system.file("extdata", "x.R", package = "myfirstpackage")
  
# Create vignette rmarkdown file
  # installing/loading the package:
  if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr
  # Installing pandoc
  install.pandoc()
  rmarkdown::find_pandoc()
  
  install.packages("jsonlite")
  install.packages("rmarkdown")
  usethis::use_vignette("introduction")
  #install pandoc
  
  
# Install your package
  devtools::install()
  # or
  devtools::load_all() # replaces installation, and attachement of package
  # remove.packages(vctrs)
  library("myfirstpackage")
  # Install a package from github:
  # devtools::install_github("yourusername/myfirstpackage")
  
  # Install from github:
  # devtools::install_github("yourusername/myfirstpackage")

  library(devtools)
  packageVersion("devtools")
  library(roxygen2)
  library(testthat)
  devtools::session_info()
  library(tidyverse)
  library(fs)
  
  load_all() # or
  pkgload::load_all()
  
  # make basic check whether the package works
  check() # includes the command test() that checks also your own created tests
  
  # Test your functions
  use_testthat() # Creates test folder
  use_test("get_IC") # creates test file for single function
  test()
  
  # Create basic read me file for git
  use_readme_rmd()
  # Render your read me file
  build_readme()
  
  
# Next Steps --------------------------------------------------------------

  #Done
  # Create basic structure of package
  # Create function
  # Add descriptions to functions
  # Add dependencies
  # Add sample data
  
  #ToDo
  # Create read me content for git
# Add correct contact information into description
# Add at all functions in the functions the package name and :: or importFrom nameofpackage nameoffunction
# Creating Vignette content
# Creating more functions and descriptions
# Change metadata
# Add automated testing (github actions)
# Create a Github pages website for the package
# package down
  
  
  # Questions:
  # Do I need to install Rtools40? And how do I install it?