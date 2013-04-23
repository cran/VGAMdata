\name{VGAMdata-package}
\alias{VGAMdata-package}
\alias{VGAMdata}
\docType{package}
\title{
Data for the VGAM package
}
\description{
  \pkg{VGAMdata} is an assortment of larger data sets which
  are a useful resource for the \pkg{VGAM} package.


}
\details{

This package contains some larger data sets originally
distributed with the \pkg{VGAM} package.
Ideally both packages can be installed and loaded to be
fully functional.
The main intent was to limit the size of \pkg{VGAM} to
a bare essential.
Many data sets in my monograph will refer to data sets
in either package.



}
\author{
Thomas W. Yee, \email{t.yee@auckland.ac.nz}.


Maintainer: Thomas Yee \email{t.yee@auckland.ac.nz}.

}

%\section{Warning}{
%  This package is undergoing continual development and improvement.
%  Until my monograph comes out and this package is released as version 1.0-0
%
%
%}


\references{


Yee, T. W.
Vector Generalized Linear and Additive Models.
\emph{Monograph in preparation}.


}


\keyword{ package }
\keyword{models}
\keyword{regression}
\seealso{
  \code{\link[VGAM]{VGAM-package}}.



%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# Example 1; xs.nz
head(xs.nz)
summary(xs.nz)

# Example 2; ugss
head(ugss)
summary(ugss)
}


