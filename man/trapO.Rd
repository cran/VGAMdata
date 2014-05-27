\name{trapO}
\alias{trapO}
\docType{data}
\title{
  Trout Data at the Te Whaiau Trap on Lake Otamangakau

  
%%   ~~ data name/kind ... ~~
}
\description{
  Rainbow and brown trout trout trapped at the
  Te Whaiau Trap at
  Lake Otamangakau in the central North Island of New Zealand.
  The data were collected by the Department of Conservation.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(trapO)}
\format{
  A data frame with 1226 observations on the following 15 variables.
  \describe{

  \item{Date}{Date as a class \code{"Date"} variable.


  }
  \item{BFTW, BMTW, RFTW, RMTW}{
  numeric vectors,
  the number of fish trapped daily.
  \code{B}/\code{R} is for brown/rainbow trout.
  \code{F}/\code{M} is for female/male.
  \code{TW} is for the Te Whaiau trap location
  (there was another trap just off the Tongariro River).


  }
  \item{MinAT, MaxAT}{numeric vectors,
    daily minimum and maximum ambient temperatures
    in Celsius.


  }
  \item{Rain}{numeric vector,
    daily rainfall that has been scaled between 0 (none)
    and 100 (flooding situation).


  }
  \item{LevelTW}{numeric vector,
  water level of the stream  that has been scaled between 0 (none)
    and 100 (flooding situation).
  In a flooding situation it is possible that some fish going
  upstream were not caught.


  }
  \item{Year, Month, Day}{numeric vectors,
  extracted from \code{Date}.


  }
  \item{doy}{a numeric vector,
  Julian day of year. The value 1 means 1st of January,
  and so on up to 365.


  }
  \item{f.Year}{a factor vector,
  the year as a factor.

  }
  \item{fict.Year}{
  similar to \code{Date} but a fictitional year
  is used for all the data. This allows all the
  data to be plotted along one calendar year.


  }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~


These are the daily numbers of fish trapped at the Te Whaiau trap
near Lake Otamangakau, during the winter months when spawning is at
its peak.
These fish were all going upstream.
There are two species of trout, split up by males and
females, in the data set. The first is brown trout (\emph{Salmo
trutta}) and the second is rainbow trout (\emph{Oncorhynchus
mykiss}). Information on the movement patterns of brown and
rainbow trout in Lake Otamangakau and Lake Te Whaiau can be found
in Dedual et al. (2000).


Brown trout are more sedentary compared with rainbow trout,
and spawning activities of brown trout occur between May and
June whilst peak spawning of rainbow trout occurs between July
and August.  Furthermore, brown trout have been observed avoiding
water above 19 degrees Celsius and optimum temperatures for growth
are between 10--15 degrees for brown trout and  16.5--17.2 degrees
for rainbow trout.


See also \code{\link[VGAM]{lakeO}}.


  
}
\source{

  Many thanks to Dr Michel Dedual (\url{http://www.doc.govt.nz})
  for making this data available.
  Help from Simeon Pattenwise is also acknowledged.


%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{


Dedual, M. and Maxwell, I. D. and Hayes, J. W. and Strickland, R. R. (2000).
Distribution and movements of brown (\emph{Salmo trutta}) and rainbow trout
(\emph{Oncorhynchus mykiss}) in Lake Otamangakau, central North Island, New
Zealand.
\emph{New Zealand Journal of Marine and Freshwater Research},
\bold{34}: 615--627.



%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data("trapO")
summary(trapO)
}
\keyword{datasets}


