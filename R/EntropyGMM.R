# theoretical Entropy for multivariate Gaussians ----
EntropyGauss <- function(sigma)
{
  sigma <- as.matrix(sigma)
  p <- ncol(sigma)
  p/2*(1 + log(2*pi)) + 0.5*log(det(sigma))
}

# generic function ----
EntropyGMM <- function(object, ...) 
{
  UseMethod("EntropyGMM")
}

# specific methods ----

EntropyGMM.densityMclust <- function(object, ...)
{
  stopifnot(inherits(object, "densityMclust"))
  -mean(log(object$density))
}

EntropyGMM.densityMclustBounded <- function(object, ...)
{
  stopifnot(inherits(object, "densityMclustBounded"))
  -mean(log(object$density))
}

EntropyGMM.Mclust <- function(object, ...)
{
  stopifnot(inherits(object, "Mclust"))
  EntropyGMM.densityMclust(as.densityMclust(object), ...)
}

EntropyGMM.data.frame <- function(object, ...)
{
  stopifnot(inherits(object, "data.frame"))
  EntropyGMM(data.matrix(object), ...)
}

EntropyGMM.matrix <- function(object, ...)
{
  data <- as.matrix(object)
  mod <- densityMclust(data, ..., plot = FALSE)
  EntropyGMM(mod)
}

# util functions ----

nats2bits <- function(x) log(exp(x), base = 2)

bits2nats <- function(x) log(2^(x), base = exp(1))

