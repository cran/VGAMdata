\name{profs.nz}
\alias{profs.nz}
\docType{data}
\title{
  Professors of Statistics in New Zealand


}
\description{
  This data set contains information on about 22
  past or present
  professors of statistics
  in New Zealand universities.


% who have been affiliated with New Zealand universities.


}
\usage{data(profs.nz)}
\format{
  A data frame with 22 observations on the following 7 variables.
  \describe{


  \item{\code{institution}}{character,
   with values \code{"MU"}, \code{"UA"}, 
    \code{"UC"}, \code{"UO"}, \code{"UW"}, \code{"VU"},
    the university affiliation.
    The abbreviations are for:
    Massey University,
    University of Auckland,
    University of Canterbury,
    University of Otago,
    University of Waikato and
    Victoria University Wellington.


  }

  \item{\code{initials}}{character,
   first and middle and surname initials.

  }

  \item{\code{firstyear}}{a numeric vector,
  the earliest indexed publication.


  }

  \item{\code{Surname}}{character, the surname.

  }


  \item{\code{pubtotal}}{a numeric vector,
   the total number of publications.


  }
  \item{\code{ARPtotal}}{a numeric vector, the total number of 
  author/related publications.

  }
  \item{\code{cites}}{a numeric vector, the number of citations.


  }

  \item{\code{ID}}{a numeric vector, the unique MR Author ID for each
      professor.

  }

  }
}
\details{
  This data set contains information taken from the MathSciNet database on
  professors of statistics (and some related fields)
  affiliated with New Zealand universities.


  In the future the following names may be added:
  C. F. Ansley,
  P. C. B. Phillips,
  B. S. Weir.




% The data is current as of 2014/06/04. 


}
\source{
  The data is compiled from
  \url{http://www.ams.org/mathscinet}
  by J. T. Gray in April 2014.


}

%\references{
%%  ~~ possibly secondary sources and usages ~~


%}



\examples{
profs.nz
profs.nz[order(with(profs.nz, pubtotal), decreasing = TRUE), ]
}
\keyword{datasets}



