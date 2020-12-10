
# First steps creating an R package ---------------------------------------

# I followed https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# Install packages
  install.packages("devtools")
  install.packages("roxygen2")

# Create Framework (NAMESPACE, DESCRIPTION,R, myfirstpackage.Rproj)
  devtools::create("myfirstpackage")
# creates automatically a new project in R

# Add dependencies to package
  usethis::use_package("data.table") # Defaults to imports
  #> Adding dplyr to Imports
  #> Refer to functions with dplyr::fun()
  usethis::use_package("data.table", "Suggests")
  #> Adding dplyr to Suggests
  #> Use requireNamespace("data.table", quietly = TRUE) to test if package is 
  #>  installed, then use data.table::fun() to refer to functions.
  
# Add function description
  
# Add .Rd file to man folder
  devtools::document()
  # Each time you add new documentation to your R function, 
  # you need to run devtools::document() again to re-generate the .Rd files.
  
# Create some sample data object and save it in data folder and inst/extdata to provide data code
  x <- c(1:10)
  usethis::use_data(x)
  
# Create vignette rmarkdown file
  install.packages("rmarkdown")
  usethis::use_vignette("introduction")
  
# Install your package
  devtools::install()
  # or
  devtools::load_all()
  library("myfirstpackage")
  # Install a package from github:
  # devtools::install_github("yourusername/myfirstpackage")
  
  

# Next Steps --------------------------------------------------------------

# Creating Vignette content
# Creating more functions and descriptions
# Change metadata
# Add automated testing 
# Create a Github pages website for the package
  