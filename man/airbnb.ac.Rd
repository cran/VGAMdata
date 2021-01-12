\name{airbnb.ac}
\alias{airbnb.ac}
\docType{data}
\title{Airbnb Accommodation in Two Sardinian Cities}
\description{
The \code{airbnb.ac} data frame has 18159 rows and 9 columns.


}
\usage{data(airbnb.ac)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{LOS}{length of stay, in days.}
    \item{discountWeek}{
Multiplicative factor for the discount for booking one week:
\eqn{1 -}
(published weekly rate) / (7 \eqn{\times}{x} published nightly rate),
e.g., 0.2 means a 20 percent savings off the regular price.

    }
    \item{NumberofReviews}{
Number of reviews on the website.

      
    }
    \item{PriceAvg}{
Average price per night per person, in Euros.

      
    }
    \item{Bedrooms}{
Number of bedrooms.

      
    }
    \item{Superhost}{
      Logical.
      Superhost?
      
    }
    \item{MinimumStay}{
      Minimum stay period, in days.

      
    }
    \item{MaxGuests}{
Maximum number of guests.
      
    }
    \item{City}{
    Character, values are \code{"Alghero"} and \code{"Cagliari"}.
      
    }
  }
}
\details{
  The data frame comprises Airbnb bookings in two
  cities located in Sardinia, Italy.
  The stays were during the whole of 2016.
  Stays of 30 days or longer
  and any rows with missing variables
  were deleted from the original source.
  Variable \code{LOS} exhibits heaping at the values 7 and 14 days.
  


}
\source{
  The data was obtained with permission
  from Luca Frigau, University of Cagliari,
  from \url{https://www.airbnb.com}.



}
\seealso{
  \code{\link[VGAM]{gaitzeta}},
  \code{\link[VGAMdata]{flamingo}}.


}

%\references{


%  Yee, T. W. and Luca, F. and Contu, G. (order to be decided on later).
%  \emph{On heaped and seeped counts,
%        with GAIT regression as a new technique
%        applied to Sardinian Airbnb data}.


%}
\examples{
\dontrun{
mytab <- with(subset(airbnb.ac, City == "Alghero"), table(LOS))
plot(prop.table(mytab), col = "blue", ylab = "Proportion")
}}
\keyword{datasets}





