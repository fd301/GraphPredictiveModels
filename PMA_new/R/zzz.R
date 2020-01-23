########################################################
###
### Define .First.lib and .Last.lib functions for the package
###
########################################################

.onLoad <- function(lib,pkg){
    library.dynam("PMA.NEW",pkg,lib)
  }
  
