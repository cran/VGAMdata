\name{pirates2}
\alias{pirates2}
\docType{data}
\title{
  Personal data of the crew of
  the ship Ranger, associated with the pirate
  Edward Low

  %%   ~~ data name/kind ... ~~
}
\description{
  A data frame containing the age, name, birthplace and
  verdict of 35 members of a pirate ship associated
  with Edward Low, who were taken to trial on 10th to 12th July, 1723.


  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(pirates2)}
\format{
  The variables \code{age} and \code{name} are
  analagous to \code{\link{pirates1}}.
  The variable \code{guilty} is binary and 1 means yes,
  0 means not guilty. Guilty crew members were executed
  except for two: John Brown and Patrick Cunningham;
  they were respited for one year and recommended to the
  King's favour.


}
\details{
  Starting on the 10th July, 1723,
  the crew of the \emph{Ranger} were judged.
  The captain of the ship was Charles Harris,
  and this ship was one of two pirate ships under Captain Edward Low. 
  Their personal data (name, age, place of birth) and verdicts
  are recorded in Johnson (1955).
  This data was was constructed from pp.295--296 of that book
  and includes those who were not found guilty (and therefore
  were not executed).
  The execution of the 25 men were performed on 19 July
  near Newport, Rhode Island, USA.
  The notorious pirate Edward Low himself was brought to
  trial in 1724 under different
  circumstances and was hanged in Martinique.
  
  
  %%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  Same as \code{\link{pirates1}}.
  This data was entered into R by Lucia Pilleri.



}
\seealso{
  \code{\link{pirates1}}.
  
}

\examples{
summary(pirates2)
}
\keyword{datasets}



