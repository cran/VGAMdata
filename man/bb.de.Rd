\name{bb.de}
\alias{bb.de}
\docType{data}
\title{
Battle of Britain Data (a Luftwaffe subset)
  %%   ~~ data name/kind ... ~~
}
\description{
  Luftwaffe losses during a subset of the Battle of Britain.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(bb.de)}
\format{
  The format is a \eqn{3}-dimensional array.
  The first dimension is the event
  (in order:
   shot down or failed to return,
   written off,
   seriously damaged),
  the second dimension is the day,
  the third is the aircraft type.
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  This is a 
  Battle of Britain data set of Luftwaffe losses during operations
  26 August--31 August 1940 continued on to 1--7 September 1940.
  The aircraft types are prefixed
  \code{Bf} for Messerschmitt (Bayerische Flugzeugwerke),
  \code{Do} for Dornier,
  \code{He} for Heinkel,
  \code{Ju} for Junkers.


  Note that p.151 and p.165 of Bowyer (1990)
  contain tables (during the first week of September) and
  almost the same data;
  however, the former is labelled "shot down" whereas
  the latter is "shot down or failed to return".
  The latter is used here.
  Also, there are some other small discrepancies.

 
  Yet to do: add the data available at other dates,
  and include the RAF data.



}
\source{

Bowyer, M. J. F. (1990)
\emph{The Battle of Britain: 50 years On}.
Patrick Stephens Limited, Northamptonshire, U.K.


%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(bb.de)
bb.de[,, "Bf109"]

\dontrun{
plot(bb.de["sdown",, "Bf109"] ~ as.Date(dimnames(bb.de)[[2]]),
     type = "h", col = "blue", las = 1, lwd = 3,
     ylab = "Frequency", xlab = "1940",
     main = "Numbers shot down (Bf 109)")
abline(h = c(5, 10, 15, 20), lty = "dashed", col = "grey")
points(bb.de["sdown",,"Bf109"] ~ as.Date(dimnames(bb.de)[[2]]), col = "blue")
}
}
\keyword{datasets}


