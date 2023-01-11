\name{flamingo}
\alias{flamingo}
\docType{data}
\title{Flamingo Hotel in Sardinia}
\description{
The \code{flamingo} data frame has 4871 rows and 12 columns.


}
\usage{data(flamingo)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{LOS}{length of stay, in days.}
    \item{year}{year of the stay.}
    \item{arrdate}{date of arrival.}
    \item{depdate}{date of departure.}
    \item{seasonindex}{
      seasonality index of the period of the stay,
      a low value means low season,
      a high value means very high (peak) season.


%     seasonality of the period of the stay,
%     low (\code{L}),
%     medium (\code{M}),
%     high (\code{H}),
%     very high (\code{VH}).



  }
  \item{booking}{if the guests booked
      the room through the hotel's
      internet website (\code{internet}), the hotel's
      telephone (\code{direct}) or
      a tour operator (\code{agency}).
}
    \item{roomtype}{the 16 hotel's room types.}
    \item{rmtype3}{categorization in 3 groups of
      \code{roomtype}.}
    \item{guests}{number of guests.
    Does not include any childrenzz;
    see the two zz \code{kids} variables below.


  }
    \item{arrangement}{arrangement type:
      Bed and Breakfast (\code{BB}),
      Half Board (\code{HB}) and Full Board (\code{FB}).
    }
    \item{kids02}{number of kids between 0 and 2 years-old.}
    \item{kids311}{number of kids between 3 and 11 years-old.}
  }
}
\details{
  The Flamingo Hotel is located on the beach in the southern Sardinia,
  about 40 kilometers from Cagliari.
  This data concerns stays from early summer 2019 to
  late summer 2020
  with no stays during the winter period in between.
  Stays longer than 30 days were deleted from the original source.
  Variable \code{LOS} exhibits heaping at the values 7 and 14 days.
  


% https://www.hotelflamingosardinia.com/



}
\source{
  The data was obtained with permission
  from Luca Frigau, University of Cagliari.



}
\seealso{
  \code{\link[VGAM]{gaitdzeta}},
  \code{\link[VGAMdata]{airbnb.ac}}.


}

%\references{


%  Yee, T. W. and Luca, F. and Contu, G. (order to be decided on later).
%  \emph{On heaped and seeped counts,
%        with GAITD regression as a new technique
%        applied to Sardinian hotel data}.


%}
\examples{
\dontrun{ with(flamingo, spikeplot(LOS, col = "blue", ylab = "Proportion")) }}
\keyword{datasets}



% zz kids311 has values 0 1 2 3
% mytab <- with(flamingo, table(LOS))
% plot(prop.table(mytab), col = "blue", ylab = "Proportion")


