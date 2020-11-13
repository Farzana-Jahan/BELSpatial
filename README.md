# Steps for working on this package

  1. Click the `BELSpatial.Rproj` file to open up a new environment.
  2. To load all of the current functions so you can test them, run `devtools::load_all()` in the console.
  3. To update the documentation after adding functions with `#'` comments above them, run `devtools::document()` or in the "Build" tab, click "More" then "Document".
  3. To check that the package can be installed, click "Check", and wait for it to check. It will fail if there are errors or warnings, so you can scroll through the output to work out the issue. Don't hesitate to send me the relevant output as well as update the git repository if you can't work out the issue.
  4. If you want to install the package you can either run `devtools::install_github("danwkenn\BELSpatial")` or `devtools::install_github("\path\to\this\directory")`.
  
There are two `.R` files in the `R/` directory containing functions. It's good practice to put similar functions in the same `.R` script file. All code needs to be put in an `.R` file in the `R/` directory. Otherwise, do what you like!

The `.R` files in there should give you an idea how to write the functions as well as how to write the documentation, using the `#'` at the start of the line.
