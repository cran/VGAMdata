\name{smqP}
\alias{smqP}
\docType{data}
\title{ Smoking Data in the NHANES 2017--2020 Cycle }

\description{
  Selected variables mainly focussed on the smoking
  questionnaire of the 
  National Health and Nutrition Examination Survey
  collected during 1.5 cycles just prior to Covid-19.

}
\usage{
data(smqP)
}
\format{
  A data frame with the following variables.
  \describe{

    \item{SEQN}{
      Identifier for individuals that can be used
      to merge in other data sets from the same cycle.


    }
    \item{TC30}{    % cigarettes
      Numeric.
      Tobacco consumption:
      average number of cigarettes/day during past 30 days.
      Aka \code{SMD650} with some preprocessing.


    }
%   \item{year_quit_smoking}{   % year_quit_smoking
%     Aka \code{SMQ050Q} with some preprocessing.

%   }
    \item{SCA}{
      Smoking cessation age (SCA).
      Is missing for those who have not quit.
      Computing this involved variables
      \code{SMQ050Q}, \code{SMQ050U} and \code{age}.
      The variable \code{SMQ050Q} is for
      \code{How long since quit smoking cigarettes?}
      and
      \code{SMQ050U} are the units
      (e.g., years, months, days).
      The variable \code{SMQ050Q} is right-censored
      at 50 years (66666 means 50 or more years)
      and for such people SCA was set to \code{NA}.


    }
    \item{TCquit}{
%   cigarettes_smoked_per_day_when_quit
Tobacco consumption:
      number of cigarettes smoked per day when quit.
      Aka \code{SMD057} with some preprocessing.

    }
    \item{SIA}{   % age_started_smoking
      Smoking initiation age (SIA):
      age when individuals started smoking cigarettes regularly.
      Aka \code{SMD030} with some preprocessing.

    }
    \item{age, gender}{
      Age (\code{RIDAGEYR}) when surveyed;
      the value 80 is right-censored.
      Gender is \code{RIAGENDR} .

      
    }
    \item{race, race2}{
      Race (\code{RIDRETH1}) and a
      binary simplification of race
      (\code{"Non-Hispanic White"} versus
       \code{"Others"}).

    }
    \item{educ, marital}{
      Education (\code{DMDEDUC2}) and
      marital status (\code{DMDMARTZ}).
      Both variables were not collected for those
      aged 12--19, hence those are \code{NA}s.



    }
    \item{poverty_ratio}{
      Ratio of family income to poverty (\code{INDFMPIR}).

    }
    \item{diet, meals_fastfood_pizza, readytoeat_foods}{
      For example, the number of meals from a
      fast food or pizza place
      (\code{DBD900}).
      May be described later in more detail.



    }
    \item{frozen_meals}{
      To be described later.


    }
    \item{people_fam_smoking, people_home_smoking}{
      Numeric.
For \code{people_fam_smoking}, this is
    Number of people who live here smoke tobacco?
    (\code{SMD460}). For this,
0 is: No one in household is a smoker;
1 is: 1 household member is a smoker;
2 is: 2 or more household members are smokers.



For \code{people_home_smoking}, this is
the number of people who smoke inside this home?
(\code{SMD470}). For this,
0 is: No one smokes inside the house;
1 is: 1 household member smokes inside the house;
2 is: 2 or more household members smoke inside the house.



    }
    \item{use_cigarettes, use_pipes, use_cigars}{
      Use during the last 5 days of smoking variants.
The codes are \code{SMQ690A}, \code{SMQ690B}, \code{SMQ690C}.

    }
    \item{use_hookah, use_ecigarettes}{
The codes are \code{SMQ690G}, \code{SMQ690H}.



    }
    \item{use_snuff, use_otherSmokeless}{
The codes are \code{SMQ690E}, \code{SMQ690K}.



    }
    \item{passiveSmoke_job, passiveSmoke_rest}{
      Binary 0 (no) or 1 (yes)
      measuring exposure to passive smoke at
      certain places in the past 7 days,
e.g.,
While you were working at a job or business outside
of the home,
did someone else smoke cigarettes or other tobacco
products indoors?
The codes are \code{SMQ858}, \code{SMQ862}.


% Variable Name: SMQ858 


    }
    \item{passiveSmoke_bar, passiveSmoke_car}{
      See above.
The codes are \code{SMQ868}, \code{SMQ872}.


% Variable Name: SMQ868


    }
    \item{passiveSmoke_home, passiveSmoke_other}{
      See above.
The codes are \code{SMQ876}, \code{SMQ880}.


    }
    \item{passiveSmoke_ecigarettes}{
      See above.
The code is \code{SMQ940}.


    }

  }
}
\details{
  The National Health and Nutrition Examination Survey
  (NHANES)
  is a well-known longitudinal study located in USA
  (a country just north of Mexico).
  This data frame shares a selection of variables mainly
  to do with the smoking questionnaire
  (codeword: \code{SMQ});
  some demographic and anthropometric variables
  might also be included and/or added later.
  The significance of \code{P} is that the
  2019--2020 cycle was not completed due to Covid-19,
  hence this data concerns 2017--2020 and is a merging of
  the 2017--2018 cycle (data codenamed \code{J}) with
  further data collected just prior to the pandemic.



  The original data has been preprocessed and/or simplified.
  For example, \code{"Don't know"} and \code{"Refused"}
  usually have been converted to a \code{NA}.



  This data frame is subject to change,
  especially with the addition of new variables.



}
\source{
  The URL
\url{https://wwwn.cdc.gov/Nchs/Nhanes/}
provides an entry point from which many
data sets may be downloaded.
Comprehensive documentation is available there too.


%\url{https://wwwn.cdc.gov/Nchs/Nhanes}





  The data was downloaded in late 2022
  by Luca Frigau, University of Cagliari,
  and some subsequent edits were made by Thomas Yee.


}
%\seealso{
% \code{\link[VGAMdata]{profs.nz}}.



%}
\examples{
summary(smqP)
}
\keyword{datasets}


%
%
