#' Prostate Cancer Data
#'
#' Data, which come from a study by Stamey et al. (1989), examining the
#' correlation between the level of prostate-specific antigen (PSA) and
#' various clinical measures in men who were about the receive a radical prostatectomy.
#'
#' @format A data frame with 97 rows and 9 variables:
#' \describe{
#'   \item{lcavol}{log(cancer volume (cm^3))}
#'   \item{lweight}{log(prostate weight (g))}
#'   \item{age}{age of the patient}
#'   \item{lbph}{log(amount of benign prostatic hyperplasia (cm^2))}
#'   \item{lcavol}{log(cancer volume (cm^3))}
#'   \item{lweight}{log(prostate weight (g))}
#'   \item{age}{age of the patient}
#'   \item{lbph}{log(amount of benign prostatic hyperplasia (cm^2))}
#'   \item{svi}{presence of seminal vesicle invasion (1=yes, 0=no)}
#'   \item{lcp}{log(capsular penetration (cm))}
#'   \item{gleason}{Gleason score}
#'   \item{pgg45}{percentage of Gleason scores four of five}
#'   \item{lpsa}{log(PSA (ng/mL))}
#'   }
#'
#' @source \url{https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data}
#' @references
#' Stamey, T. A., Kabalin, J. N., McNeal, J. E., Johnstone, I. M., Freiha, F., Redwine,
#' E. A., and Yang, N. (1989). Prostate specific antigen in the diagnosis and treatment of
#' adenocarcinoma of the prostate. ii. radical prostatectomy treated patients. The Journal
#' of urology, 141(5):1076-1083.
#'
"pcancer"

#' Sniffer Data
#'
#' Data examining the factors that impact the amount of hydrocarbon vapour released
#' when gasoline is pumped into a tank.
#'
#' @format A data frame with 125 rows and 5 variables:
#' \describe{
#'   \item{tanktemp}{initial tank temperature (degrees F)}
#'   \item{gastemp}{temperature of the dispensed gasoline (degrees F)}
#'   \item{tankpres}{initial vapour pressure in the tank (psi)}
#'   \item{gaspres}{vapour pressure of the dispensed gasoline (psi)}
#'   \item{y}{hydrocarbons emitted (g)}
#'   }
#'
#' @source \url{https://CRAN.R-project.org/package=alr4}
#' These data were obtained from the \code{alr4} package.
#'
#' @references
#' Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.
"sniffer"

#' Boston House Price Data
#'
#' Data, which come from a study by Harrison Jr and Rubinfeld (1978), examining
#' the association between median house prices in a particular community with
#' various community characteristics.
#'
#' @format A data frame with 506 rows and 9 variables:
#' \describe{
#'    \item{crime}{crimes committed per capita}
#'    \item{rooms}{average number of rooms per house}
#'    \item{radial}{index of accessibility to radial highways}
#'    \item{stratio}{average student-teacher ratio of schools in the community}
#'    \item{lowstat}{percentge of the population that are "lower status"}
#'    \item{lnox}{log(annual average nitrogen oxide concentration (pphm))}
#'    \item{lproptax}{log(property tax per $1000)}
#'    \item{ldist}{log(weighted distances to five employment centres in the Boston region)}
#'    \item{lprice}{log(median house price ($))}
#'    }
#'
#' @source \url{https://CRAN.R-project.org/package=wooldridge}
#' These data were obtained from the \code{wooldridge} package.
#'
#' @references
#' Harrison Jr, D. and Rubinfeld, D. L. (1978). Hedonic housing prices and the
#' demand for clean air. Journal of environmental economics and management, 5(1):81-102.
#'
#' Wooldridge, J. M. (2015). Introductory econometrics: A modern approach. Cengage learning.
#'
"bostonhouseprice"
