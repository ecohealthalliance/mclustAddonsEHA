.onAttach <- function(lib, pkg)
{
  # invisible(suppressPackageStartupMessages(
  #   requireNamespace("mclust", quietly = TRUE)
  # ))
  # startup message
  msg <- paste("Loaded package 'mclustAddons' version", packageVersion("mclustAddons"))
  packageStartupMessage(msg)      
  invisible()
}
