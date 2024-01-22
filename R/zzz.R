.onAttach <- function(lib, pkg)
{
  # invisible(suppressPackageStartupMessages(
  #   requireNamespace("mclust", quietly = TRUE)
  # ))
  # startup message
  msg <- paste("Loaded package 'mclustAddonsEHA' version", packageVersion("mclustAddonsEHA"))
  packageStartupMessage(msg)      
  invisible()
}
