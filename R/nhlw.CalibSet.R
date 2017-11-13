##' NHL cancer in women for districts covered by a cancer registry.
##'
##' Data on NHL cancer in women over the 2007-2011 period for the 16 districts covered by a
##' cancer registry.
##'
##' @format A data frame with  256  rows and  5  variables
##' \itemize{
##' \item \code{dist} District
##' \item \code{age} Central age of the age class
##' \item \code{C} Number of cancer incident cases
##' \item \code{H} Number of hospitalizations
##' \item \code{py} Number of person-years
##' }
##' @docType data
##' @keywords datasets
##' @name  nhlw.CalibSet
##' @details The number of cancer incident cases were provided
##' by the network of French cancer registries FRANCIM.
##' The number of cancer incident cases,
##' corresponding numbers of hospitalizations and person-years
##' are tabulated by 5 years age groups (from 0-5 to 90+,
##' variable `age` corresponds to the central age of the class) and
##' district (variable `dist`). In order to diminish the number of age
##' classes with no cancer, the number of incident cases in lower ages
##' classes were aggregated.
##' @usage data( nhlw.CalibSet )
NULL
