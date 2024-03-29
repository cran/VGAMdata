\name{xs.nz}
\alias{xs.nz}
\docType{data}
\title{
  Cross-sectional Data from the New Zealand Population

}
\description{
  A cross-sectional data set of a workforce company, plus
  another health survey, in New Zealand during the 1990s,


%%  ~~ A concise (1-5 lines) description of the dataset. ~~
%%  zz on the following 64zz  variables.... adjust according to
%   the true number of columns.



}
\usage{data(xs.nz)}
\format{
  A data frame with 10529 observations on the
  following 64 variables.
  For binary variables, a \code{"1"} or \code{TRUE}
  means \code{yes},
  and \code{"0"} or \code{FALSE} means \code{no}.
  Also, \code{"D"} means don't know,
  and \code{"-"} means not applicable.
  The pregnancy questions were administered to women only.


  \describe{
    \item{\code{regnum}}{a numeric vector,
      a unique registration number.
      This differs from their original registration number,
      and the rows are sorted by their new registration number.

    }
    \item{\code{study1}}{a logical vector,
      Study 1 (workforce) or Study 2?

    }
    \item{\code{age}}{a numeric vector, age in years.

    }
    \item{\code{sex}}{a factor with levels \code{F} and \code{M}.

    }
    \item{\code{pulse}}{a numeric vector, beats per minute.

    }
    \item{\code{sbp}}{a numeric vector,
      systolic blood pressure (mm Hg).

    }
    \item{\code{dbp}}{a numeric vector,
      diastolic blood pressure
      (mm Hg).

    }
    \item{\code{cholest}}{a numeric vector, cholesterol (mmol/L).

    }
    \item{\code{height}}{a numeric vector, in m.

    }
    \item{\code{weight}}{a numeric vector, in kg.

    }
    \item{\code{fh.heartdisease}}{a factor
      with levels \code{0}, \code{1},
      \code{D}.
      Has a family history of heart disease
      (heart attack, angina, or
      had a heart bypass operation) within the immediate
      family (brother, sister, father or mother,
      blood relatives only)?
      Note that \code{D} means: do not know.


    }
    \item{\code{fh.age}}{a factor, following
      from \code{fh.heartdisease},
      if yes, how old was the family member when
      it happened (if
      more than one family member, give the age of the
      youngest person)?


    }
    \item{\code{fh.cancer}}{a factor with levels \code{0}, \code{1},
      \code{D}.
      Has a family history of cancer within the immediate
      family (blood relatives only)?
      Note that \code{D} means: do not know.


    }
    \item{\code{heartattack}}{a numeric vector,
      have you ever been told by
      a doctor that you have had a heart attack ("coronary")?


    }
    \item{\code{stroke}}{a numeric vector, have you ever been told by
      a doctor that you have had a stroke?

    }
    \item{\code{diabetes}}{a numeric vector,
      have you ever been told by
      a doctor that you have had diabetes?


    }
    \item{\code{hypertension}}{a numeric vector,
      have you ever been told by
      a doctor that you have had high blood pressure
      (hypertension)?

    }
    \item{\code{highchol}}{a numeric vector, have you ever been
      told by a doctor that you have had high cholesterol?

% 20140709


    }
    \item{\code{asthma}}{a numeric vector,
      have you ever been told by
      a doctor that you have had asthma?


    }
    \item{\code{cancer}}{a numeric vector,
      have you ever been told by
      a doctor that you have had cancer?

    }
    \item{\code{acne}}{a numeric vector, have you ever
      received treatment from a doctor for acne (pimples)?


    }
    \item{\code{sunburn}}{a numeric vector, have you ever
      received treatment from a doctor for sunburn?

    }
    \item{\code{smokepassive}}{a numeric vector, on average,
      how many hours each week (at work and at home) would you
      spend near someone who is smoking?
      (put \code{"0"} if none)


% 20220202


    }
    \item{\code{smokeever}}{a numeric vector, have you ever
      smoked tailor-made or roll-you-own cigarettes once a
      week or more?
      A 1 means yes and 0 means no.


    }
    \item{\code{smokenow}}{a numeric vector,
      do you smoke tailor-made or roll-you-own cigarettes now?
      A 1 means yes and 0 means no.


    }
    \item{\code{smokeagequit}}{a factor,
      if no to \code{smokenow}, how old were you when
      you stopped smoking?
      Using \code{as.numeric(as.character(smokeagequit))}
      will work for those values which are not
      \code{as.character(smokeagequit) == "-"}.


    }
    \item{\code{smokeyears}}{a numeric vector,
      if yes to \code{smokeever}, for how many years altogether
      have you smoked tailor-made or roll-you-own cigarettes?


    }
    \item{\code{smoketailormade}}{a numeric vector,
      how many tailor-made cigarettes do you smoke each day?



% 20211218



    }
    \item{\code{smokeweekpack}}{a numeric vector, how many
      packets of \emph{roll-your-own} tobacco do you use
      each week?
      (put \code{"0"} if none)

% 20220202

    }
    \item{\code{smokepacketsize}}{a numeric vector,
      what size packets of \emph{roll-your-own} tobacco do you
      usually buy?
      (\code{"0"} means don't smoke \emph{roll-your-owns},
      else 25g or 30g or 35g or 50g)

% 20220202

    }
    \item{\code{drinkmonth}}{a numeric vector,
      do you drink alcohol once a month or more?


% 20140709 was alcmonth


    }
    \item{\code{drinkfreqweek}}{a numeric vector,
      if yes to \code{drinkmonth}, about how often do you
      drink alcohol (days per week)?
      Note: 0.25 is once a month,
      0.5 is once every two weeks,
      1 is once a week,
      2.5 is 2-3 days a week,
      4.5 is 4-5 days a week,
      6.5 is 6-7 days a week.


      Further note:
      1 can, small bottle or
      handle of beer or home brew = 1 drink,
      1 quart bottle of beer = 2 drinks,
      1 jug of beer = 3 drinks,
      1 flagon/peter of beer = 6 drinks,
      1 glass of wine, sherry = 1 drink,
      1 bottle of wine = 6 drinks,
      1 double nip of spirits = 1 drink.


    }
    \item{\code{drinkweek}}{a numeric vector,
      how many drinks per week, on average.
      This is the average daily amount of drinks multiplied
      by the frequency of drinking per week.
      See \code{drinkfreqweek} on what constitutes a 'drink'.



    }
    \item{\code{drinkmaxday}}{a numeric vector,
      in the last three months, what is the largest number of
      drinks that you had on any one day?
      Warning: some values are considered
      unrealistically excessive.


    }
    \item{\code{eggs}}{a numeric vector,
      how many eggs do you eat a week
      (raw, boiled, scrambled, poached, or in quiche)?


    }
    \item{\code{chocbiscuits}}{a numeric vector,
      how many chocolate biscuits
      do you usually eat in a week?


    }
    \item{\code{pregnant}}{a factor,
      have you ever been pregnant for more than 5 months?


    }
    \item{\code{pregfirst}}{a factor, if
      yes to \code{pregnant}, how old were you when your first
      baby was born (or you had a miscarriage after 5 months)?

    }
    \item{\code{preglast}}{a factor, how old
      were you when your last
      baby was born (or you had a miscarriage after 5 months)?


    }
    \item{\code{babies}}{numeric,
      how many babies have you given birth to?


    }
    \item{\code{moody}}{a numeric vector,
      does your mood often go up or down?


    }
    \item{\code{miserable}}{a numeric vector,
      do you ever feel 'just miserable' for no reason?


    }
    \item{\code{hurt}}{a numeric vector,
      are your feelings easily hurt?

    }
    \item{\code{fedup}}{a numeric vector,
      do you often feel 'fed up'?

    }
    \item{\code{nervous}}{a numeric vector,
      would you call yourself a nervous person?

    }
    \item{\code{worrier}}{a numeric vector,
      are you a worrier?

    }
    \item{\code{worry}}{a numeric vector,
      do you worry about awful things that might happen?

    }
    \item{\code{tense}}{a numeric vector,
      would you call yourself tense or 'highly strung'?

    }
    \item{\code{embarrassed}}{a numeric vector,
      do you worry too long after an embarrassing
      experience?

    }
    \item{\code{nerves}}{a numeric vector,
      do you suffer from 'nerves'?

    }
    \item{\code{nofriend}}{a numeric vector,
      do you have a friend or family member that you
      can talk to about problems or worries that you may have?
      The value 1 effectively means \code{"no"},
      i.e., s/he has no friend or friends.


% 20140709



    }
    \item{\code{depressed}}{a numeric vector,
      in your lifetime, have you ever had two weeks or more
      when nearly every day you felt sad or depressed?

    }
    \item{\code{exervig}}{a numeric vector,
      how many hours per week would you do any vigorous
      activity or exercise either at work or away from
      work that makes you breathe hard and sweat?
      Values here ought be be less than 168.


    }
    \item{\code{exermod}}{a numeric vector,
      how many hours per week would you do any moderate
      activity or exercise such as brisk walking, cycling or
      mowing the lawn?
      Values here ought be be less than 168.


    }
    \item{\code{feethour}}{a numeric vector,
      on an average work day, how long would you spend on your
      feet, either standing or moving about?


% 20140709


    }
    \item{\code{ethnicity}}{a factor with 4 levels,
      what ethnic group do you belong to?
      \code{European} = European (NZ European or
      British or other European),
      \code{Maori} = Maori,
      \code{Polynesian} = Pacific Island Polynesian,
      \code{Other} = Other (Chinese, Indian, Other).


% 20140709
%     0 = European (NZ European or British or other European),
%     1 = Maori,
%     2 = Pacific Island Polynesian,
%     3 = Other (Chinese, Indian, Other).


    }
    \item{\code{sleep}}{a numeric vector,
      how many hours do you usually sleep each night?


    }
    \item{\code{snore}}{a factor with levels \code{0}, \code{1},
      \code{D}.
      Do you usually snore?
      Note that \code{D} means: do not know.


    }
    \item{\code{cat}}{a numeric vector,
      do you have a household pet cat?


    }
    \item{\code{dog}}{a numeric vector,
      do you have a household pet dog?


    }
    \item{\code{hand}}{a factor with levels
      \code{right} = right,
      \code{left} = left,
      \code{both} = either.
      Are you right-handed,
      left-handed,
      or no preference for left or right?


%     \code{0} = right,
%     \code{1} = left,
%     \code{2} = either.



    }
    \item{\code{numhouse}}{an ordered factor with 4 levels:
      \code{1} = 1,
      \code{2} = 2,
      \code{3} = 3,
      \code{4+} = four or more;
      how many people (including yourself)
      usually live in your house?


    }
    \item{\code{marital}}{a factor with 4 levels:
         \code{single} = single,
         \code{married} = married or living with a partner,
         \code{separated} = separated or divorced,
         \code{widowed} = widowed.



%        \code{1} = single,
%        \code{2} = married or living with a partner,
%        \code{3} = separated or divorced,
%        \code{4} = widowed.


  }
  \item{\code{educ}}{an ordered factor with 4 levels:
          \code{primary} = Primary school,
          \code{secondary} = High school/secondary school,
          \code{polytechnic} = Polytechnic or similar,
          \code{university} = University.
    What was the highest level of education you received?


%         Primary school = \code{1},
%         High school/secondary school = \code{2},
%         Polytechnic or similar = \code{3},
%         University =  \code{4}.


  }
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~

The data frame is a subset of the entire data set which was
collected from a confidential self-administered questionnaire
administered in a large New Zealand workforce observational
study conducted during 1992--3. The data were augmented
by a second study consisting of retirees. The data can be
considered a reasonable representation of the white male New
Zealand population in the early 1990s.  There were physical,
lifestyle and psychological variables that were measured.
The psychological variables were headed
"Questions about your feelings".



Although some data cleaning was performed and logic checks
conducted, anomalies remain. Some variables, of course,
are subject to a lot of measurement error and bias.  It is
conceivable that some participants had poor reading skills!
In particular, the smoking variables contain a small
percentage of conflicting values, and when \code{NA}s are taken
into account then there would be several different ways
the data might be cleaned.
If \code{smokeever == 0} then strictly speaking, only
\code{smokepassive} is the other variable---the other
smoking variables should either be \code{NA} or \code{0}.




}



\source{

  Originally,
  Clinical Trials Research Unit,
  University of Auckland, New Zealand,
  \code{http://www.ctru.auckland.ac.nz}.
  Originally much of the error checking and formatting was
  performed by Stephen Vander Hoorn.
  Lately (2014), more changes and error checks were made to the
  data by James T. Gray.


  %%  ~~ reference to a publication or URL from
  %% which the data were obtained ~~

%  MacMahon, S., Norton, R., Jackson, R., Mackie, M. J.,
%  Cheng, A., Vander Hoorn, S., Milne, A., McCulloch, A. (1995)
%  Fletcher Challenge-University of Auckland Heart &
%  Health Study: design and baseline findings.
%  \emph{New Zealand Medical Journal},
%  \bold{108}, 499--502.



}
\references{

  MacMahon, S., Norton, R., Jackson, R., Mackie, M. J.,
  Cheng, A., Vander Hoorn, S., Milne, A., McCulloch, A. (1995).
  Fletcher Challenge-University of Auckland Heart &
  Health Study: design and baseline findings.
  \emph{New Zealand Medical Journal},
  \bold{108}, 499--502.

}

\seealso{
  \code{\link[VGAM]{chest.nz}}.



}
\section{Warning }{
  More variables may be added in the future and these may
  be placed in any column position. Therefore references
  such as \code{xs.nz[, 12]} are dangerous.
  Also, variable names may change in the future as well as
  their format or internal structure,
  e.g., \code{factor} versus \code{numeric}.


}


\note{
  More error checking are needed for the
  pregnancy
  and
  smoking variables.



}


\examples{
data(xs.nz)
summary(xs.nz)
}
\keyword{datasets}


% Warning: some values are considered unrealistically excessive.


%# Handling of factors requires care
%is.factor(xs.nz$babies)  # TRUE
%summary(xs.nz$babies)  # Note the "-"s
%charbabies <- as.character(xs.nz$babies)
%summary(as.numeric(charbabies))  # "-"s converted to NAs+warning
%table(as.numeric(charbabies))  # Ditto



