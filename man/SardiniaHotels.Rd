\name{SardiniaHotels}
\alias{SardiniaHotels}
\docType{data}
\title{
  Data from hotels in Sardinia, Italy

}

\description{
This data set contains information and satisfaction scores
appearing on the TripAdvisor website between the years 2008 and 2016
regarding hotels in Sardinia, Italy.


The satisfaction data refer to the reputation of hotel located along
Sardinian coasts, as expressed by clients with respect to different
services (e.g., breakfast, restaurant, swimming pool)
offered by the hotel.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}

\usage{data(SardiniaHotels)}

\format{
  A data frame with 518 rows and 43 columns (variables).
  Each row refers to a single hotel.


The following variables are included in the dataset:


\describe{

\item{\code{municipality}}{a factor, the municipality where
the hotel is located.

}

\item{\code{stars}}{an ordered factor with levels:
  
  \code{1OR2stars} for 1 star or 2 star hotels,

  \code{3stars} 3 star hotels, 
  
  \code{residence},
  
  \code{4stars}, 4 star hotels,
  
  \code{5starsORresort}, 5 star hotels or resorts.
    
}

    \item{\code{area}}{a factor with levels related to the
    area of the Sardinian coast where each single hotel
    is located:

\code{AlgheroSassari},
\code{CagliariVillasimius}, \code{CostaSmeralda},
\code{DorgaliOrosei}, \code{Gallura}, \code{NurraAnglona},
\code{Ogliastra}, \code{Olbia}, \code{OristanoBosa},
\code{PulaChia}, \code{Sarrabus},
\code{Sulcis}.

    }

    \item{\code{seaLocation}}{a factor with levels
    \code{yes} (if the hotel is located close to the sea)
    and \code{no} (otherwise).

    }

    \item{\code{excellent}}{a numeric vector, the number of
    people that expressed the highest level of satisfaction.

  }

  \item{\code{good}}{a numeric vector, the number
    of people that expressed a good level of satisfaction.

  }

  \item{\code{average}}{a numeric vector, the number of
    people that expressed an average level of satisfaction.

  }
  
  \item{\code{bad}}{a numeric vector, the number of
    people that expressed a bad level of satisfaction.

  }

  
    \item{\code{poor}}{a numeric vector, the number of
    people that expressed the lowest level of satisfaction.

  }
  
    \item{\code{family}}{a numeric vector, the number
    of people travelling with family.

  }
  
    \item{\code{couple}}{a numeric vector, the number
    of people travelling with
    their partner.

  }
  
\item{\code{single}}{a numeric vector, the number of
    people travelling alone.

}

\item{\code{business}}{a numeric vector, the number of
people travelling for work.

}

\item{\code{MarMay}}{a numeric vector, the number of
people travelling during the period March to May.

}

\item{\code{JunAug}}{a numeric vector, the number of
people travelling during the period
June to August.

}

\item{\code{SepNov}}{a numeric vector, the number of
people travelling during the period
September to November.

}

\item{\code{DecFeb}}{a numeric vector, the number of
people travelling during the period
December to February.

}

\item{\code{location}}{a numeric vector, the satisfaction
score expressed by tourists towards the location.

}

\item{\code{sleepQuality}}{a numeric vector, the
satisfaction score expressed by tourists towards the
sleep quality.

}

\item{\code{room}}{a numeric vector, the satisfaction score
expressed by tourists towards the comfort and quality
of the room.

}

\item{\code{services}}{a numeric vector, the satisfaction
score expressed by tourists towards the quality of
the services.

}

\item{\code{priceQualityRate}}{a numeric vector, the
satisfaction score expressed by tourists towards ratio
between price and quality.

}

\item{\code{cleaning}}{a numeric vector, the satisfaction
score expressed by tourists towards level of room and
hotel cleaning.

}

\item{\code{bt1}}{a factor with levels \code{breakfast},
\code{cleaning}, \code{location}, \code{overall},
\code{price}, \code{restaurant}, \code{room},
\code{services}, \code{staff}, \code{structure} and
\code{Wi-Fi}.

It expresses the 1st most used word in reviews for a hotel.
 
}

\item{\code{ratebt1}}{a factor with levels \code{-1} (if the
satisfaction score espressed in \code{bt1} is prevalently negative)
and \code{1} (if the satisfaction score espressed in \code{bt1}
is prevalently positive).

    
}

\item{\code{bt2}}{a factor with levels \code{breakfast},
\code{cleaning}, \code{location}, \code{overall},
\code{price}, \code{restaurant}, \code{room},
\code{services}, \code{staff}, \code{structure} and
\code{Wi-Fi}.


It expresses the 2nd most used word in reviews for a hotel.
}

\item{\code{ratebt2}}{a factor with levels \code{-1} (if the
satisfaction score espressed in \code{bt2} is prevalently negative)
and \code{1} (if the satisfaction score espressed in \code{bt2}
is prevalently positive).


}

\item{\code{bt3}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt4}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt5}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt6}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt7}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt8}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt9}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}
\item{\code{bt10}}{ similar to \code{bt1} and \code{bt2},
                    but with a corresponding different ranking.}

\item{\code{ratebt3}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt4}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt5}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt6}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt7}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt8}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt9}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}
\item{\code{ratebt10}}{similar
      to \code{ratebt1} and \code{ratebt2}, but with a
      corresponding different ranking.}

}
}

\details{


These data were manually collected during
March--June 2016 by students of the class of "Statistics for
Turism" at the University of Cagliari, Italy (Bachelor's
degree in Tourism Economics and Managment), under the
supervision of Prof. Claudio Conversano and Dr. Giulia Contu.


Many of the variables fall into several natural groups, e.g.,
[\code{municipality},
\code{stars},
\code{area},
\code{seaLocation}];
[\code{excellent},
\code{good},
\code{average},
\code{bad},
\code{poor}];
[\code{MarMay},
\code{JunAug},
\code{SepNov},
\code{DecFeb}];
[\code{family},
\code{couple},
\code{single},
\code{business}];
[\code{location},...\code{cleaning}];
[\code{bt1},...\code{bt10}];
[\code{ratebt1},...\code{ratebt10}].
}


\source{

  TripAdvisor,
  \url{https://www.tripadvisor.it/}.

}

%\references{


%}


%\seealso{
%  \code{\link[VGAM]{SardiniaHotels}}.


%}
%\section{Warning }{



%}

%\note{


%}


\examples{
data(SardiniaHotels)
summary(SardiniaHotels)
}
\keyword{datasets}


  
      
    
    
    
    
    
    
