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
#'
#' @references
#' Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley.
"sniffer"

#' Boston House Price Data (Original)
#'
#' Original data, which come from a study by Harrison Jr and Rubinfeld (1978), examining
#' the association between median house prices in a particular community with
#' various community characteristics. See \code{\link{bostonhouseprice2}}
#' for the corrected version, with additional variables.
#'
#' @format A data frame with 506 rows and 9 variables:
#' \describe{
#'    \item{crime}{crimes committed per capita}
#'    \item{rooms}{average number of rooms per house}
#'    \item{radial}{index of accessibility to radial highways}
#'    \item{stratio}{average student-teacher ratio of schools in the community}
#'    \item{lowstat}{percentage of the population that are "lower status"}
#'    \item{lnox}{log(annual average nitrogen oxide concentration (pphm))}
#'    \item{lproptax}{log(property tax per $1000)}
#'    \item{ldist}{log(weighted distances to five employment centres in the Boston region)}
#'    \item{lprice}{log(median house price ($))}
#'    }
#'
#' @source \url{https://CRAN.R-project.org/package=wooldridge}
#'
#' @references
#' Harrison Jr, D. and Rubinfeld, D. L. (1978). Hedonic housing prices and the
#' demand for clean air. Journal of environmental economics and management, 5(1):81-102.
#'
#' Wooldridge, J. M. (2015). Introductory econometrics: A modern approach. Cengage learning.
#'
"bostonhouseprice"

#' Boston House Price Data (Corrected Version)
#'
#' Corrected data, which come from a study by Harrison Jr and Rubinfeld (1978), examining
#' the association between median house prices in a particular community with
#' various community characteristics. See \code{\link{bostonhouseprice}} for the
#' original version.
#'
#' @format A data frame with 506 rows and 13 variables:
#' \describe{
#'    \item{crim}{per capita crime rate by town}
#'    \item{zn}{proportion of residential land zoned for lots over 25,000 sq.ft}
#'    \item{indus}{proportion of non-retail business acres per town}
#'    \item{rm}{average number of rooms per dwelling}
#'    \item{age}{proportion of owner-occupied units built prior to 1940}
#'    \item{rad}{index of accessibility to radial highways}
#'    \item{ptratio}{pupil-teacher ratio by town}
#'    \item{lnox}{log(nitric oxides concentration (parts per 10 million))}
#'    \item{ldis}{log(weighted distances to five Boston employment centres)}
#'    \item{ltax}{log(full-value property-tax rate per USD 10,000)}
#'    \item{llstat}{log(percentage of lower status of the population)}
#'    \item{chast}{Charles River dummy variable (=1 if tract bounds river; 0 otherwise)}
#'    \item{lcmedv}{log(corrected median value of owner-occupied homes in USD 1000's)}
#'    }
#'
#' @source \url{https://CRAN.R-project.org/package=mlbench}
#'
#' @references
#' Harrison Jr, D. and Rubinfeld, D. L. (1978). Hedonic housing prices and the
#' demand for clean air. Journal of environmental economics and management, 5(1):81-102.
#'
#' Leisch F, Dimitriadou E (2021). mlbench: Machine Learning Benchmark Problems. R package version 2.1-3.
#'
"bostonhouseprice2"

#' Diabetes Data
#'
#' Data relating to a study of disease progression one year after baseline.
#'
#' @format A data frame with 442 rows and 11 variables:
#' \describe{
#'   \item{AGE}{age of the patient}
#'   \item{SEX}{sex of the patient}
#'   \item{BMI}{body mass index of the patient}
#'   \item{BP}{blood pressure of the patient}
#'   \item{S1}{blood serum measurement 1}
#'   \item{S2}{blood serum measurement 2}
#'   \item{S3}{blood serum measurement 3}
#'   \item{S4}{blood serum measurement 4}
#'   \item{S5}{blood serum measurement 5}
#'   \item{S6}{blood serum measurement 6}
#'   \item{Y}{quantitative measure of disease progression one year after baseline}
#'   }
#'
#' @source \url{https://CRAN.R-project.org/package=lars}
#'
#' @references
#' Efron, B., Hastie, T., Johnstone, I., Tibshirani, R., et al. (2004).
#' Least angle regression. The Annals of Statistics.
"diabetes"
