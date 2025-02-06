\name{birds.df}
\alias{birds.df}
\docType{data}
\title{
A Bird Data Set
}
\description{
  An undocumented data set concerning birds.


  
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("birds.df")}
\format{
  A data frame with 67 observations on the following 15 variables.
  \describe{
    \item{\code{Species}}{a character vector}
    \item{\code{Status}}{a numeric vector}
    \item{\code{Length}}{a numeric vector}
    \item{\code{Mass}}{a numeric vector}
    \item{\code{Range}}{a numeric vector}
    \item{\code{Migr}}{a numeric vector}
    \item{\code{Insect}}{a numeric vector}
    \item{\code{Diet}}{a numeric vector}
    \item{\code{Clutch}}{a numeric vector}
    \item{\code{Broods}}{a numeric vector}
    \item{\code{Wood}}{a numeric vector}
    \item{\code{Upland}}{a numeric vector}
    \item{\code{Water}}{a numeric vector}
    \item{\code{Release}}{a numeric vector}
    \item{\code{Indiv}}{a numeric vector}
  }
}
\details{
  This is a data set used in a book (or in a
  draft version of a book).
  Unfortunately very few details seem
  available concerning the background and
  the variables.
  The purpose of its inclusion in
  \pkg{VGAMdata} is to show that
  a logistic regression applied to these data,
  mimicking the book,
  has estimates suffering from the
  Hauck--Donner effect (HDE;
  \code{\link[VGAM]{hdeff}}).
  It is called \code{model.final} in the
  book.



}
\source{
  \code{http://rcompanion.org/documents/RCompanionBioStatistics.pdf},
  around pages 249--.



% \url{http://rcompanion.org/documents/RCompanionBioStatistics.pdf},



  
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(birds.df)
str(birds.df)

model.final.copy <-   # Try mimic the book example
  vglm(Status ~ Upland + Migr + Mass + Indiv + Insect + Wood,
       binomialff, data = birds.df, trace = TRUE)
summary(model.final.copy, hde.NA = FALSE)  # Wald table
hdeff(model.final.copy)
hdeff(model.final.copy, deriv = 2)
}
\keyword{datasets}


