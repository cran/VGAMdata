\name{belcap}
\alias{belcap}
\docType{data}
\title{
  BELCAP Dental Data


%%   ~~ data name/kind ... ~~
}
\description{
  A prospective data set containing the DMFT index of children in
  Belo Horizonte at the beginning and end of the BELCAP study.


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(belcap)}
\format{
  A data frame with 797 observations on the following 5 variables.
  \describe{
    \item{\code{dmftb}}{a numeric vector, DMFT-Index at the beginning of the study.}
    \item{\code{dmfte}}{a numeric vector, DMFT-Index at the end of the study.}
    \item{\code{gender}}{a factor with levels \code{0} = female, \code{1} = male.}
    \item{\code{ethnic}}{a factor with levels \code{1} = dark, \code{2} = white, 
  \code{3} = black.}
    \item{\code{school}}{the kind of prevention measure.
  A factor with levels \code{1} = 
  oral health education, \code{2} = all four methods, \code{3} = control group, 
  \code{4} = enrichment of the school diet with ricebran, \code{5} = mouthrinse 
  with 0.2\% NaF-solution, \code{6} = oral hygiene.}
  }

}
\details{
  This data set is from the Belo Horizonte Caries Prevention (BELCAP) study. The
  data is collected from children in Belo Horizonte (Brazil) aged seven years at 
  the start of the study. In in order to determine which method(s) were best for 
  preventing tooth decay, six treatments were randomized to six separate schools.
  The measure used is the decayed, missing and filled teeth (DMFT) index - a well
  known and important measure of a persons dental health. Only the eight deciduous
  molars are considered, so the lowest value is 0, and the highest is 8. 
 
 
 
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{

\code{https://onlinelibrary.wiley.com/}
contains the data file
(a supplement of the JRSSA article).
Downloaded in January 2014 and formatted into R by J. T. Gray,
\code{jamsgry@gmail.com}.



% url{http://onlinelibrary.wiley.com} contains the data file



%%  ~~ reference to a publication or URL from which the data were obtained ~~



}
\references{


Bohning, D., D. Ekkehart, P. Schlattmann, L. Mendonca, and U. Kircher (1999). 
The Zero-Inflated Poisson Model and the Decayed, Missing and Filled Teeth Index
in Dental Epidemiology,
\emph{Journal of the Royal Statistical Society, A} \bold{162}(2), 195--209.



%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(belcap)
## maybe str(belcap) ; plot(belcap) ...
}
\keyword{datasets}



% \url{http://onlinelibrary.wiley.com/journal/10.1111/%28ISSN%291467-985X/homepage/162_2.htm}
