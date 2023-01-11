\name{rar.df}
\alias{rar.df}
\docType{data}
\title{ Rock and Roll and Life Expectancy }

\description{
  Ages at death in rock and roll.


}
\usage{
data(rar.df)
}
\format{
  A data frame with the following variables.

  \describe{

    \item{age}{
      Numeric.
      Age when the person died, in years.
      A handful of values are approximate or missing.

      

    }
    \item{yod}{
      Numeric.
      The year that the person died.

    }

  }
}
\details{
  This type of data is very subjective because many
  definitions cannot ever be finalized.
  The Wiki website was written by an anonymous person
  and the data collected informally and there is no
  quality assurance.
  Of course, it could be debated as whether a particular person
  should be included.
  And there will be many people who should be there but
  are not.
  Nevertheless, it can be considered a fun data set 
  that should be of interest to anybody liking
  modern music or
  the history of such in musicology.



  In the source, some peoples' ages were listed as a
  decade, e.g., age was \code{40s} or \code{50--51}.
  These were either replaced by \code{NA} or the mean
  of the endpoints was taken.
  

  

}
\source{
  The lists at
\url{https://en.wikipedia.org/wiki/List_of_deaths_in_rock_and_roll}.
were downloaded and edited.



}
%\seealso{
% \code{\link[VGAMdata]{profs.nz}}.



%}
\examples{
summary(rar.df)
}
\keyword{datasets}


%
%
