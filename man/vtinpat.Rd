\name{vtinpat}
\alias{vtinpat}
\docType{data}
\title{
Vermont Hospital Inpatient Data
}
\description{
Information on inpatients discharged from hospitals in Vermont, USA, 2012.
}
\usage{data(vtinpat)}
\format{
  A data frame with 52206 observations on the following 7 variables.
  \describe{
    \item{\code{hospital}}{a factor with levels \code{1} =  Northwestern Medical Center, \code{2} = North Country Hospital and Health Center,
      \code{3} = Northeastern Vermont Regional Hospital, \code{4} = Copley Hospital, \code{5} = Fletcher Allen Health Care,
      \code{6} = Central Vermont Hospital, \code{8} = Rutland Regional Medical Center, \code{9} = Porter Medical Center,
      \code{10} = Gifford Memorial Hospital, \code{11} = Mount Ascutney Hospital and Health Center, \code{12} = Springfield Hospital,
      \code{14} = Grace Cottage Hospital, \code{15} = Brattleboro Memorial Hospital, \code{16} = Southwestern Vermont Medical Center}
    \item{\code{admit}}{a factor with levels \code{1} = Emergency, \code{2} = Urgent, \code{3} = Elective, \code{4}, Newborn,
      \code{5} = Trauma}
    \item{\code{age.group}}{a factor with levels \code{1} = Under 1, \code{2} = 1-17, \code{3} = 18-24, \code{4} = 25-29, \code{5} = 30-34,
      \code{6} = 35-39, \code{7} = 40-44, \code{8} = 45-49, \code{9} = 50-54, \code{10} = 55-59, \code{11} = 60-64, \code{12} = 65-69,
      \code{13} = 70-74, \code{14} = 75+}
    \item{\code{sex}}{a factor with levels \code{1} = Male, \code{2} = Female}
    \item{\code{discharge}}{a factor with levels \code{1} = To another medical facility, \code{2} = home, \code{3} = against medical advice,
      \code{4} = Died, \code{5} = To court or law enforcement, \code{6} = still a patient}
    \item{\code{diagnosis}}{a factor with levels \code{1} = Brain And C.N.S., \code{2} = Eye, \code{3} = Ear, Nose & Throat,
      \code{4} = Respiratory, \code{5} = Heart & Circulatory, \code{6} = Digestive, \code{7} = Liver & Pancreas, \code{8} = Musculoskeletal,
	  \code{9} = Skin and Breast, \code{10} = Endocrine, \code{11} = Kidney & Urinary, \code{12} = Male Reproductive,
	  \code{13} = Female Reproductive, \code{14} = Pregnancy, Childbirth, \code{15} = Neonatal, \code{16} = Spleen & Blood,
	  \code{17} = Lymphatic, \code{18} = Infection, \code{19} = Mental Illness, \code{20} = Substance Abuse, \code{21} = Injury, Toxic Effects,
	  \code{22} = Burns, \code{23} = Other, \code{24} = Trauma, \code{25} = H.I.V.}
    \item{\code{los}}{a numeric vector, number of days spent in hospital}
  }
}
\details{
  This data set contains details on inpatients discharged from
  hospitals in Vermont, USA, in 2012 as part of the
  Vermont Uniform Hospital Discharge Data Set.
  The Vermont Department of Financial Regulation administers this
  program and the Vermont Department of Health manages the data set.



}

\source{
  Vermont department of Health,
  \url{https://www.healthvermont.gov/stats/systems}
  formatted into R by J. T. Gray in mid-2014.


% the URL \code{http://healthvermont.gov/research/hospital-utilization/RECENT_PU_FILES.aspx}




}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
summary(vtinpat)
}
\keyword{datasets}
% ## maybe str(vtinpat) ; plot(vtinpat) ...
