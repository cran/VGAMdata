\name{rugby}
\alias{rugby}
\alias{rugby.ties}
\docType{data}
\title{Wins, Losses and Draws Between the Top 10 Rugby Teams
%%   ~~ data name/kind ... ~~
}
\description{
  The number of wins, losses and draws for each of the top 10 
  rugby teams agains each other
  
  
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
data(rugby)
data(rugby.ties)
}

\format{
  The format is as two matrices.


}
\details{
  The first matrix is of the number of games won by each team 
  against each of the other teams.
  The other matrix is the number of draws (ties) between each team.
  This data is current as of 2013-10-07.



%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
The match statistics are compiled from
\url{http://www.rugbydata.com/}
on 2013-10-07 by J. T. Gray,
\code{jamsgry@gmail.com}.


The top ten teams are determined by the International Rugby Board
world rankings, \url{http://www.irb.com}.

%%  ~~ reference to a publication or URL from which the data were obtained ~~

}
\examples{
data(rugby); data(rugby.ties)
rugby
rugby.ties
}
\keyword{datasets}



