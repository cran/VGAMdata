\name{crim.nz}
\alias{crim.nz}
\docType{data}
\title{New Zealand Conviction and
  Sentencing Data
  Subset 2001--2022


% for Fairness Analysis

%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~

These data were collected from the
New Zealand Ministry of Justice and
comprises selected years from
2001 to 2022.



%This data frame came
%from the



}
\usage{data("crim.nz")}
\format{
  A data frame with 257600 observations
  on the following 7 variables.
  \describe{
    \item{\code{prison}}{a logical vector;
  \code{FALSE} if and only if the sentence was either
  \code{"Imprisonment"},
  \code{"Imprisonment sentences"} or
  \code{"Life Imprisonment"}.
  So the other sentences were not imprisonment.


    }
    \item{\code{agegp}}{an ordered factor with
      levels \code{17-19} < \code{20-24} <
      \code{25-29} < \code{30-39} <
      \code{40+}.
   The age of the offender.
      
    }
    \item{\code{offence}}{a factor with
      levels
      \code{abduction},
      \code{injury},
      \code{endanger},
      \code{fraud},
      \code{homicide},
      \code{drugs},
      \code{miscoff},
      \code{antigovt},
      \code{weapons},
      \code{property},
      \code{order},
      \code{robbery},
      \code{sexoff},
      \code{theft},
      \code{burglary}.
      The \emph{main} offence.
      In more detail,
      the levels correspond to
\code{"Abduction, harassment and other offences against the person"},
\code{"Acts intended to cause injury"},
\code{"Dangerous or negligent acts endangering persons"},
\code{"Fraud, deception and related offences"},
\code{"Homicide and related offences"},
\code{"Illicit drug offences"},
\code{"Miscellaneous offences"},
\code{"Offences against justice procedures, gov. security and gov. operations"},
\code{"Prohibited and regulated weapons and explosives offences"},
\code{"Property damage and environmental pollution"},
\code{"Public order offences"},
\code{"Robbery, extortion and related offences"},
\code{"Sexual assault and related offences"},
\code{"Theft and related offences"},
\code{"Unlawful entry with intent/burglary, break and enter"},             
respectively.



  }
  \item{\code{gender}}{a factor with levels 
    \code{F} and \code{M}.


  }
  \item{\code{ethnicity}}{a factor with levels
    \code{Asian},
    \code{European},
    \code{Maori},
    \code{Other},
    \code{Polynesian}.


%   \code{Polynesian} (Pacifica).



  }
  \item{\code{year}}{a numeric vector.
    The calendar year when the crime occurred.


  }
  \item{\code{sentence}}{a factor with levels
    \code{ComDet},
    \code{ComSent},
    \code{ComWork},
    \code{HomeDet},
    \code{Prison},
    \code{PrisonSent},
    \code{IntSupv},
    \code{Life},
    \code{Money},
    \code{Other},
    \code{PrevDet},
    \code{Supv}.
    In more detail,
    the levels correspond to
    \code{"Community Detention"},
    \code{"Community sentences"},
    \code{"Community Work"},
    \code{"Home Detention"},
    \code{"Imprisonment"},
    \code{"Imprisonment sentences"},
    \code{"Intensive Supervision"},
    \code{"Life Imprisonment"},
    \code{"Monetary "},
    \code{"Other"},
    \code{"Preventive Detention"},
    \code{"Supervision"} respectively.


 
  }
  }
}
\details{

The data were collected in late 2023 and is
described in detail in Garcia et al. (2023).
Almost all the information here comes from
that document.
The original data comprised each year from
2001 to 2022 inclusive, however only the years
2001, 2010, 2019 and 2022 are included here
for brevity.

  




  Variables with values
  \code{"Unknown/Organisation"}
  were treated as missing and deleted.
  The \code{\link[base]{rownames}} were
  stripped off to make the data frame
  much smaller.

  

A matching
New Zealand Police data set on proceedings
against offenders was also collected but
not included in \pkg{VGAMdata} because of
its size.



Offences are categorised using the Australian
and New Zealand Standard Offence Classification
(ANZSOC).  There are sixteen top level
`divisions' which are further subdivided into
`subdivisions' and `groups'.  One of the sixteen
main divisions (Traffic Offences) are excluded.
Only the most serious offence by the offender
is recorded here.


Each row contains the count of sentences
received by offenders with particular
demographic information in the given period,
along with the ANZSOC offence code. However,
since ANZSOC is a hierarchical classification, a
unique offender is represented in three separate
rows: one for the ANZSOC Group of the offence,
one for the Subdivision of the offence, and
one for the Division of the offence.

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~


\url{https://datainfoplus.stats.govt.nz/}




}
\references{
%%  ~~ possibly secondary sources and usages ~~


Fairness analysis and machine learning model
implementation using New Zealand prosecution
and conviction data.
Garcia, F. and Omidyar, P. and Rodrigues, L.
and Xie, J. (2023).
Postgraduate Report,
Computer Science Department,
University of Auckland.



Only ANZSOC Divisions 01--16 were collected.
More details can be found at
\url{https://www.police.govt.nz/about-us/publications-statistics/data-and-statistics/policedatanz/proceedings-offender-demographics}
and
\url{https://www.abs.gov.au/statistics/classifications/australian-and-new-zealand-standard-offence-classification-anzsoc/latest-release}




%\url{https://datainfoplus.stats.govt.nz/item/nz.govt.stats/a10413bf-f78a-4f17-a9c1-55e7717ab91d}.






  


}
\examples{
data(crim.nz)
summary(crim.nz)
}
\keyword{datasets}
