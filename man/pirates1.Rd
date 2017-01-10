\name{pirates1}
\alias{pirates1}
\docType{data}
\title{
  Personal data of the executed pirates
  associated with Bartholomew Roberts
  
  %%   ~~ data name/kind ... ~~
}
\description{
  The age, names and habitation of 52
  pirates who were found guilty of piracy and executed, after
  the ships associated with Bartholomew Roberts were captured.


% executed, according to their sentence, without the
% Gates of Cape Corso Castle, within the flood-marks, viz.


  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(pirates1)}
\format{
  A data frame with the following 3 variables.


  \describe{
    \item{\code{age}}{a numeric vector,
      their age in years at the time of trial.
      Bartholomew Roberts himself was 39 years old at his death.


    }
    \item{\code{name}}{character.

    }
    \item{\code{habitation}}{character.

    }
  }
}
\details{
  According to Wiki,
  in February 1722 Captain Ogle was sent by
the British Government to find and capture
the notorious pirate
Bartholomew Roberts
(real name: John Roberts, but also known
later as Black Bart).
When his warship
caught up with the \emph{Royal Fortune} he attacked and
Bartholomew Roberts was the first to fall, followed by
2 others.
The remaining pirates surrendered soon afterwards.
A total of 272 men were captured, and
of these, 65 were
black, and they were sold into slavery. The remainder were taken to Cape
Coast Castle, apart from those who died on the voyage back.
The trial was held in April, 1722, and
54 were condemned to death, of whom 52 were hanged and two were reprieved.
Of those executed, their personal data (name, age, habitation)
were recorded.


  
  %%  ~~ If necessary, more details than the __description__ above ~~
}
\source{ Pages 248--249 of
  Johnson, Captain Charles, (1955)
  (Editor: Arthur L. Hayward).
  \emph{A General History of the Robberies and Murders of the Most Notorious
    Pirates},
  London: Routledge and Kegan Paul Ltd.
  This edition was first published in 1926.
  The earliest manuscript of the book dates back to 1724.



  This data was entered into R by Lucia Pilleri.


% \url{https://books.google.it/books?id=z81WnQEACAAJ}




  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%  \emph{}.
  
  
  %%  ~~ possibly secondary sources and usages ~~
%}
\seealso{
  \code{\link{pirates2}}.
  
}

\examples{
summary(pirates1)
}
\keyword{datasets}



