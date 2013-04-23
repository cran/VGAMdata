\name{students.tw}
\alias{students.tw}
\docType{data}
\title{
  Taiwanese students answer a multiple response question

}
\description{
  This data is a subset from a survey of 49609 first-year college
  students in Taiwan collected in the year 2003 about their
  preferences for college study.


}
\usage{data(students.tw)}
\format{
  A data frame with 49609 observations on the following 12 response
  variables.
  For binary variables, a \code{"1"} means \code{yes},
  and \code{"0"} means \code{no}.
  See below for exact wording (translated from the original Chinese).

  
  \describe{
    \item{\code{ID}}{a numeric vector, a unique identification number
      for each student in the survey.

    }
    \item{\code{read}}{Read Chinese and foreign classics.    

    }
    \item{\code{t.travel}}{Travel around Taiwan.      

    }
    \item{\code{conference}}{Present academic papers in conferences.

    }
    \item{\code{act.leader}}{Lead large-scale activities.

    }
    \item{\code{team}}{Be on a school team.

    }
    \item{\code{stu.leader}}{Be a student association leader.

    }
    \item{\code{intern}}{Participate internship programs.

    }
    \item{\code{love}}{Fall in love.

    }
    \item{\code{sex}}{Have sexual experience.

    }
    \item{\code{o.travel}}{Travel around the world.

    }
    \item{\code{friends}}{Make many friends.

    }
    \item{\code{other}}{Other experience which is not included in the survey.

    }
}
}
    
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

This data frame is a subset of a larger data set
where any student with any missing value was deleted.
The remaining data set contains of 32792 students. 
Unfortunately, other variables such as age and sex were not
made available.



Each student was asked the following
multiple response question.


Question : What kind of experience do you expect to receive during 
the period of college study? (Select at least one response)

       1. Read Chinese and foreign classics

       2. Travel around Taiwan

       3. Present academic papers in conferences

       4. Lead large-scale activities

       5. Be on a school team

       6. Be a student association leader

       7. Participate internship programs

       8. Fall in love

       9. Have sexual experience

      10. Travel around the world

      11. Make many friends

      12. Other


}



\source{

  Originally, the data set for was downloaded from a survey center
  of Academia Sinica \url{https://srda.sinica.edu.tw/news}.
  It now seems unavailable.


% the data set for Bayesian ranking responses was downloaded from a survey


}
\references{

  Wang, H. and Huang, W. H. (2013)
  Bayesian Ranking Responses in Multiple Response Questions.
  \emph{Journal of the Royal Statistical Society, Series A},
  \bold{(to appear)}.


%  \emph{wang@stat.nctu.edu.tw}.



  Help from Viet Hoang Quoc is gratefully acknowledged.


}


\examples{
data(students.tw)
summary(students.tw)

with(students.tw, table(love, sex))
\dontrun{
plot(jitter(sex) ~ jitter(love), data = students.tw, col = "blue",
     main = "Taiwanese students")
}
}
\keyword{datasets}



