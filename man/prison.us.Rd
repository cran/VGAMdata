\name{prison.us}
\alias{prison.us}
\docType{data}
\title{
  US Prison Data
%%   ~~ data name/kind ... ~~
}
\description{
  Number of prisoners in each North American state,
  and the populations of those states from years 1977 to 2010
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(prison.us)}
\format{
  A data frame with 34 observations on the following 103 variables.
  \describe{
    \item{Year}{a numeric vector, the year}
    \item{AL.num, AL.pop}{numeric vectors}
    \item{AK.num, AK.pop, AZ.num}{numeric vectors}
    \item{AZ.pop, AR.num, AR.pop}{numeric vectors}
    \item{CA.num, CA.pop, CO.num}{numeric vectors}
    \item{CO.pop, CT.num, CT.pop}{numeric vectors}
    \item{DE.num, DE.pop, FL.num}{numeric vectors}
    \item{FL.pop, GA.num, GA.pop}{numeric vectors}
    \item{HI.num, HI.pop, ID.num}{numeric vectors}
    \item{ID.pop, IL.num, IL.pop}{numeric vectors}
    \item{IN.num, IN.pop, IA.num}{numeric vectors}
    \item{IA.pop, KS.num, KS.pop}{numeric vectors}
    \item{KY.num, KY.pop, LA.num}{numeric vectors}
    \item{LA.pop, ME.num, ME.pop}{numeric vectors}
    \item{MD.num, MD.pop, MA.num}{numeric vectors}
    \item{MA.pop, MI.num, MI.pop}{numeric vectors}
    \item{MN.num, MN.pop, MS.num}{numeric vectors}
    \item{MS.pop, MO.num, MO.pop}{numeric vectors}
    \item{MT.num, MT.pop, NE.num}{numeric vectors}
    \item{NE.pop, NV.num, NV.pop}{numeric vectors}
    \item{NH.num, NH.pop, NJ.num}{numeric vectors}
    \item{NJ.pop, NM.num, NM.pop}{numeric vectors}
    \item{NY.num, NY.pop, NC.num}{numeric vectors}
    \item{NC.pop, ND.num, ND.pop}{numeric vectors}
    \item{OH.num, OH.pop, OK.num}{numeric vectors}
    \item{OK.pop, OR.num, OR.pop}{numeric vectors}
    \item{PA.num, PA.pop, RI.num}{numeric vectors}
    \item{RI.pop, SC.num, SC.pop}{numeric vectors}
    \item{SD.num, SD.pop, TN.num}{numeric vectors}
    \item{TN.pop, TX.num, TX.pop}{numeric vectors}
    \item{UT.num, UT.pop, VT.num}{numeric vectors}
    \item{VT.pop, VA.num, VA.pop}{numeric vectors}
    \item{WA.num, WA.pop, WV.num}{numeric vectors}
    \item{WV.pop, WI.num, WI.pop}{numeric vectors}
    \item{WY.num, WY.pop        }{numeric vectors}
    \item{US.pop, US.num}{numeric vectors,
    overall counts for the whole country }
  }
}
\details{
  This is a data set of the number of prisoners in each American
  state and the populations of those states, from 1977 to 2010.
  The number of prisoners are taken from December 31st, while
  the populations are estimates taken from July 1st based on the
  previous Census, except for pop.1980, which uses exact census
  data from 1980/04/01.
%%  ~~ If necessary, more details than the __description__ above ~~


  Warning:
  a scatterplot of \code{US.pop} shows a discontinuity around 2000.


}
\source{
The prisoner data was compiled from:
Bureau of Justice Statistics,
\url{http://www.bjs.gov/index.cfm}.
Downloaded in September 2013 and formatted into R by J. T. Gray,
\code{jamsgry@gmail.com}.


The population data was compiled from:
United States Census Bureau,
\code{http://www.census.gov/popest/data}.
Downloaded in September 2013 by J. T. Gray.
This site may have become stale.
%%  ~~ reference to a publication or URL from which the data were obtained ~~



}
\examples{
summary(prison.us)
\dontrun{ # This plot shows a discontinuity around 2000.
plot(US.pop / 1e6 ~ Year, prison.us, main = "US population (millions)",
     las = 1, type = "b", col = "blue") }
}
\keyword{datasets}
