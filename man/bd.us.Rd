\name{bd.us}
\alias{bd.us}
\docType{data}
\title{
Births and Deaths of 348 Notable Americans

%%   ~~ data name/kind ... ~~
}
\description{
  A 12 x 12 table of the Births and Deaths of 348 Notable Americans.
  The rows and columns are for each month.


}
\usage{data(bd.us)}
\format{
  The format is:
 chr "bd.us"

}
\details{
  Rows denote the month of birth;
  columns for the month of death.
  These data appear as Table 1 of Phillips and Feldman (1973),
  who collected the data from Morris (1965).
  Not all of the 400 people were used because some had not died by that
  time and other individuals lacked the information.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Phillips, D. P. and Feldman, K. A. (1973)
  A dip in deaths before ceremonial occasions:
  Some new relationships between social integration and mortality,
  {American Sociological Review},
  \bold{38}, 678--696.



  Morris, R. B. (Ed.) (1965)
  \emph{Four Hundred Notable Americans}.
  Harper and Row: New York, USA.
  

}


\seealso{
  \code{\link[VGAM]{rcim}}.


}

%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
print(bd.us)
sum(bd.us)
rowSums(bd.us)
colSums(bd.us)
}
\keyword{datasets}
