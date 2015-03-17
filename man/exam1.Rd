\name{exam1}
\alias{exam1}
\docType{data}
\title{
  Examination data
%%   ~~ data name/kind ... ~~
}
\description{
  Exam results of 35 students on 18 questions.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~


}
\usage{data(exam1)}
\format{
  A data frame with 35 observations on the following 18 variables.
  \describe{
%   \item{name}{Name of the student }
    \item{q01, q02, q03, q04, q05, q06}{binary response}
    \item{q07, q08, q09, q10, q11, q12}{binary response}
    \item{q13, q14, q15, q16, q17, q18}{binary response}
  }
}
\details{
  For each question, a 1 means correct, a 0 means incorrect.
  A simple Rasch model may be fitted to this dataframe using
  \code{\link[VGAM]{rcim}} and \code{\link[VGAM]{binomialff}}.


}
\source{

Taken from William Revelle's \emph{Short Guide to R},
\url{http://www.unt.edu/rss/rasch_models.htm},
\url{http://www.personality-project.org/r/}.
Downloaded in October 2013.


}
\examples{
summary(exam1)  # The names of the students are the row names

# Fit a simple Rasch model.
# First, remove all questions and people who were totally correct or wrong
exam1.1 <- exam1  [, colMeans(exam1  ) > 0]
exam1.1 <- exam1.1[, colMeans(exam1.1) < 1]
exam1.1 <- exam1.1[rowMeans(exam1.1) > 0, ]
exam1.1 <- exam1.1[rowMeans(exam1.1) < 1, ]
Y.matrix <- rdata <- exam1.1

\dontrun{ # The following needs: library(VGAM)
rfit <- rcim(Y.matrix, family = binomialff(multiple.responses = TRUE),
             trace = TRUE)

coef(rfit)  # Row and column effects
constraints(rfit, matrix = TRUE)  # Constraint matrices side-by-side
dim(model.matrix(rfit, type = "vlm"))  # 'Big' VLM matrix
}

\dontrun{ # This plot shows the (main) row and column effects
par(mfrow = c(1, 2), las = 1, mar = c(4.5, 4.4, 2, 0.9) + 0.1)
saved <- plot(rfit, rcol = "blue", ccol = "orange",
              cylab = "Item effects", rylab = "Person effects",
              rxlab = "", cxlab = "")

names(saved@post)  # Some useful output put here
cbind(saved@post$row.effects)
cbind(saved@post$raw.row.effects)
round(cbind(-saved@post$col.effects), dig = 3)
round(cbind(-saved@post$raw.col.effects), dig = 3)
round(matrix(-saved@post$raw.col.effects, ncol = 1,  # Rename for humans
             dimnames = list(colnames(Y.matrix), NULL)), dig = 3)
}
}
\keyword{datasets}

% coef(rfit, matrix = TRUE)


