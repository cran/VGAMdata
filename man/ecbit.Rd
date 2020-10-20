\name{ecb06.it}
\alias{ecb06.it}
\alias{ecb14.it}
\docType{data}
\title{Italian Household Data for 2006 and 2014}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Part of the data collected
at two time points (2006 and 2014)
by the Bank of Italy, as part of
the European Central Banks Eurosystem collection of statistics,
within the periodic 
sample surveys on households, businesses and selected intermediaries.


}
%\usage{
%data(ecb06.it)
%data(ecb14.it)
%}
\format{
  Data frame with the following 20 variables:
  \describe{
    \item{\code{ID}}{a numeric vector, a unique identification number of the
    household.}
    
    \item{\code{area}}{a factor with 5 levels,
    the Italian  geographic area or region
    in which the household lives:
        \code{NW} = North-West,
        \code{NE} = North-East,
        \code{C} = Center,
        \code{S} = South,
        \code{I} = Islands.
	For users wanting a North--South contrast,
	this variable might be coded as
        \code{NC} = North and Center (NW, NE and C),
        \code{SI} = South and Islands (S and I).
        }
    
	\item{\code{sex}}{a factor with 2 levels,
	  the gender of the head householder:
        \code{M} = Male,
        \code{F} = Female.
        }
    
    \item{\code{age}}{a numeric vector, age in years of the head householder.}
    
    \item{\code{marital}}{a factor with 4 levels, marital status of the head 
    householder: 
        \code{married} = Married,
        \code{single} = Single,
        \code{separated} = Separated or divorced,
        \code{widowed} = Widowed.
        }

    \item{\code{education}}{an ordered factor with 8 levels,
    the education level of the head householder:
        \code{none} = No education,
        \code{primaryschool} = Primary school,
        \code{middleschool} = Middle school,
        \code{profschool} = Professional school,
        \code{highschool} = High school,
        \code{bachelors} = Bachelors degree,
        \code{masters} = Masters degree,
        \code{higherdegree} = Higher degree.
        }
    \item{\code{job}}{a factor with 7 levels,
        the type of job done by the head householder:
        \code{worker} = Worker,
        \code{employee} = Employee,
        \code{manager} = Manager,
        \code{business} = Business person,
        \code{other} = Other kind of independent job,
        \code{retired} = Retired,
        \code{unemployed} = Unemployed.
        }
    
    \item{\code{occupants}}{a numeric vector, the number of people living in 
    the same house.}
    
    \item{\code{children}}{a numeric vector, the number of children of the head
    householder living with him/her.}
    
    \item{\code{other.children}}{a numeric vector, the number of children of the
    head householder not living with the household.}
    
    \item{\code{house.owned}}{a numeric vector, the ownership of the house in 
    which the householder lives;
        \code{0} = The house is not owned,
        \code{1} = The house is owned.
        }
    
    \item{\code{houses}}{a numeric vector, the number of houses owned 
    by the family, including the one in which the family lives.}
    
    \item{\code{earners}}{a numeric vector, the number of people in the house 
    who have salary or some kind of earnings.}

    \item{\code{accounts}}{a numeric vector, the number of bank accounts 
     collectively owned by the household.}
    
    \item{\code{ccards}}{a numeric vector, the number of credit cards 
    collectively owned by the household.}
    
    \item{\code{tot.income},
          \code{dep.income},
          \code{pens.income},
          \code{self.income},
          \code{cap.income}}{numeric vectors, the amount of income 
    (in Euros) collectively earned by the household through different 
    activities. The variables can be negative if the household has installments. 
    In order, they are the total amount of income,
    the amount of income earned through dependent working,
    the amount of income earned by the household through pensions,
    the amount of income earned by the household through self-dependent working,
    the amount of income earned by the household through capital investments.
    }
  }

}
\details{
  %%  ~~ If necessary, more details than the __description__ above ~~
  The European Central Banks (ECB) Eurosystem requests and
  helps organize
  each country within the European Union to routinely collect statistical
  information via their central banks.
  These data frames are a subset from data collected by the
  Bank of Italy.
  Each year can be considered a cross-sectional study, although
  there are some people common to each year. Hence the data
  collectively can be considered partly a longitudinal study too.



}
\source{
Data was downloaded at \code{https://www.bancaditalia.it} in May 2016 by Lucia 
Pilleri.


%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
\emph{Supplements to the Statistical Bulletin,
Sample Surveys, Household Income and Wealth in 2006},
New series, Volume XVIII, Number 7--28, January 2008.
Banca D'Italia, Centro Stampa, Roma, Pubbl. Mensile,
\code{https://www.bancaditalia.it}.



}
\examples{
data(ecb06.it); data(ecb14.it)
summary(ecb06.it)
summary(ecb14.it)
\dontrun{
with(ecb14.it, table(house.owned))
with(ecb14.it, barplot(table(education), col = "lightblue"))
}
}
\keyword{datasets}
