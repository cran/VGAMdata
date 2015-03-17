\name{tube10}
\alias{tube10}
\docType{data}
\title{
  London Underground (Tube) Passenger Counts for November 2010
}
\description{
  The data set contains counts of the number of passengers entering
  London Underground stations
  (also known as \emph{the Tube})
  on a typical day of November 2010 in quarter-hour time blocks. 
}
\usage{data("tube10")}
\format{
  A data frame with 100 observations on the following 270 variables.
  \describe{
    \item{\code{ActonTown}}{a numeric vector}
    \item{\code{Aldgate}}{a numeric vector}
    \item{\code{AldgateEast}}{a numeric vector}
    \item{\code{Alperton}}{a numeric vector}
    \item{\code{Amersham}}{a numeric vector}
    \item{\code{Angel}}{a numeric vector}
    \item{\code{Archway}}{a numeric vector}
    \item{\code{ArnosGrove}}{a numeric vector}
    \item{\code{Arsenal}}{a numeric vector}
    \item{\code{BakerStreet}}{a numeric vector}
    \item{\code{Balham}}{a numeric vector}
    \item{\code{BankAndMonument}}{a numeric vector}
    \item{\code{Barbican}}{a numeric vector}
    \item{\code{Barking}}{a numeric vector}
    \item{\code{Barkingside}}{a numeric vector}
    \item{\code{BaronsCourt}}{a numeric vector}
    \item{\code{Bayswater}}{a numeric vector}
    \item{\code{Becontree}}{a numeric vector}
    \item{\code{BelsizePark}}{a numeric vector}
    \item{\code{Bermondsey}}{a numeric vector}
    \item{\code{BethnalGreen}}{a numeric vector}
    \item{\code{Blackfriars}}{a numeric vector}
    \item{\code{BlackhorseRoad}}{a numeric vector}
    \item{\code{BondStreet}}{a numeric vector}
    \item{\code{Borough}}{a numeric vector}
    \item{\code{BostonManor}}{a numeric vector}
    \item{\code{BoundsGreen}}{a numeric vector}
    \item{\code{BowRoad}}{a numeric vector}
    \item{\code{BrentCross}}{a numeric vector}
    \item{\code{Brixton}}{a numeric vector}
    \item{\code{BromleyByBow}}{a numeric vector}
    \item{\code{BuckhurstHill}}{a numeric vector}
    \item{\code{BurntOak}}{a numeric vector}
    \item{\code{CaledonianRoad}}{a numeric vector}
    \item{\code{CamdenTown}}{a numeric vector}
    \item{\code{CanadaWater}}{a numeric vector}
    \item{\code{CanaryWharf}}{a numeric vector}
    \item{\code{CanningTown}}{a numeric vector}
    \item{\code{CannonStreet}}{a numeric vector}
    \item{\code{CanonsPark}}{a numeric vector}
    \item{\code{ChalfontAndLatimer}}{a numeric vector}
    \item{\code{ChalkFarm}}{a numeric vector}
    \item{\code{ChanceryLane}}{a numeric vector}
    \item{\code{CharingCross}}{a numeric vector}
    \item{\code{Chesham}}{a numeric vector}
    \item{\code{Chigwell}}{a numeric vector}
    \item{\code{ChiswickPark}}{a numeric vector}
    \item{\code{Chorleywood}}{a numeric vector}
    \item{\code{ClaphamCommon}}{a numeric vector}
    \item{\code{ClaphamNorth}}{a numeric vector}
    \item{\code{ClaphamSouth}}{a numeric vector}
    \item{\code{Cockfosters}}{a numeric vector}
    \item{\code{Colindale}}{a numeric vector}
    \item{\code{ColliersWood}}{a numeric vector}
    \item{\code{CoventGarden}}{a numeric vector}
    \item{\code{Croxley}}{a numeric vector}
    \item{\code{DagenhamEast}}{a numeric vector}
    \item{\code{DagenhamHeathway}}{a numeric vector}
    \item{\code{Debden}}{a numeric vector}
    \item{\code{DollisHill}}{a numeric vector}
    \item{\code{EalingBroadway}}{a numeric vector}
    \item{\code{EalingCommon}}{a numeric vector}
    \item{\code{EarlsCourt}}{a numeric vector}
    \item{\code{EastActon}}{a numeric vector}
    \item{\code{EastFinchley}}{a numeric vector}
    \item{\code{EastHam}}{a numeric vector}
    \item{\code{EastPutney}}{a numeric vector}
    \item{\code{Eastcote}}{a numeric vector}
    \item{\code{Edgware}}{a numeric vector}
    \item{\code{EdgwareRoadBak}}{a numeric vector}
    \item{\code{EdgwareRoadCir}}{a numeric vector}
    \item{\code{ElephantAndCastle}}{a numeric vector}
    \item{\code{ElmPark}}{a numeric vector}
    \item{\code{Embankment}}{a numeric vector}
    \item{\code{Epping}}{a numeric vector}
    \item{\code{Euston}}{a numeric vector}
    \item{\code{EustonSquare}}{a numeric vector}
    \item{\code{Fairlop}}{a numeric vector}
    \item{\code{Farringdon}}{a numeric vector}
    \item{\code{FinchleyCentral}}{a numeric vector}
    \item{\code{FinchleyRoad}}{a numeric vector}
    \item{\code{FinsburyPark}}{a numeric vector}
    \item{\code{FulhamBroadway}}{a numeric vector}
    \item{\code{GantsHill}}{a numeric vector}
    \item{\code{GloucesterRoad}}{a numeric vector}
    \item{\code{GoldersGreen}}{a numeric vector}
    \item{\code{GoldhawkRoad}}{a numeric vector}
    \item{\code{GoodgeStreet}}{a numeric vector}
    \item{\code{GrangeHill}}{a numeric vector}
    \item{\code{GreatPortlandStreet}}{a numeric vector}
    \item{\code{GreenPark}}{a numeric vector}
    \item{\code{Greenford}}{a numeric vector}
    \item{\code{Gunnersbury}}{a numeric vector}
    \item{\code{Hainault}}{a numeric vector}
    \item{\code{HammersmithDis}}{a numeric vector}
    \item{\code{HammersmithHC}}{a numeric vector}
    \item{\code{Hampstead}}{a numeric vector}
    \item{\code{HangerLane}}{a numeric vector}
    \item{\code{Harlesden}}{a numeric vector}
    \item{\code{HarrowAndWealdstone}}{a numeric vector}
    \item{\code{HarrowOnTheHill}}{a numeric vector}
    \item{\code{HattonCross}}{a numeric vector}
    \item{\code{HeathrowTerminals123}}{a numeric vector}
    \item{\code{HeathrowTerminal4}}{a numeric vector}
    \item{\code{HeathrowTerminal5}}{a numeric vector}
    \item{\code{HendonCentral}}{a numeric vector}
    \item{\code{HighBarnet}}{a numeric vector}
    \item{\code{HighStreetKensington}}{a numeric vector}
    \item{\code{HighburyAndIslington}}{a numeric vector}
    \item{\code{Highgate}}{a numeric vector}
    \item{\code{Hillingdon}}{a numeric vector}
    \item{\code{Holborn}}{a numeric vector}
    \item{\code{HollandPark}}{a numeric vector}
    \item{\code{HollowayRoad}}{a numeric vector}
    \item{\code{Hornchurch}}{a numeric vector}
    \item{\code{HounslowCentral}}{a numeric vector}
    \item{\code{HounslowEast}}{a numeric vector}
    \item{\code{HounslowWest}}{a numeric vector}
    \item{\code{HydeParkCorner}}{a numeric vector}
    \item{\code{Ickenham}}{a numeric vector}
    \item{\code{Kennington}}{a numeric vector}
    \item{\code{KensalGreen}}{a numeric vector}
    \item{\code{KensingtonOlympia}}{a numeric vector}
    \item{\code{KentishTown}}{a numeric vector}
    \item{\code{Kenton}}{a numeric vector}
    \item{\code{KewGardens}}{a numeric vector}
    \item{\code{Kilburn}}{a numeric vector}
    \item{\code{KilburnPark}}{a numeric vector}
    \item{\code{KingsCrossStPancras}}{a numeric vector}
    \item{\code{Kingsbury}}{a numeric vector}
    \item{\code{Knightsbridge}}{a numeric vector}
    \item{\code{LadbrokeGrove}}{a numeric vector}
    \item{\code{LambethNorth}}{a numeric vector}
    \item{\code{LancasterGate}}{a numeric vector}
    \item{\code{LatimerRoad}}{a numeric vector}
    \item{\code{LeicesterSquare}}{a numeric vector}
    \item{\code{Leyton}}{a numeric vector}
    \item{\code{Leytonstone}}{a numeric vector}
    \item{\code{LiverpoolStreet}}{a numeric vector}
    \item{\code{LondonBridge}}{a numeric vector}
    \item{\code{Loughton}}{a numeric vector}
    \item{\code{MaidaVale}}{a numeric vector}
    \item{\code{ManorHouse}}{a numeric vector}
    \item{\code{MansionHouse}}{a numeric vector}
    \item{\code{MarbleArch}}{a numeric vector}
    \item{\code{Marylebone}}{a numeric vector}
    \item{\code{MileEnd}}{a numeric vector}
    \item{\code{MillHillEast}}{a numeric vector}
    \item{\code{MoorPark}}{a numeric vector}
    \item{\code{Moorgate}}{a numeric vector}
    \item{\code{Morden}}{a numeric vector}
    \item{\code{MorningtonCrescent}}{a numeric vector}
    \item{\code{Neasden}}{a numeric vector}
    \item{\code{NewburyPark}}{a numeric vector}
    \item{\code{NorthActon}}{a numeric vector}
    \item{\code{NorthEaling}}{a numeric vector}
    \item{\code{NorthGreenwich}}{a numeric vector}
    \item{\code{NorthHarrow}}{a numeric vector}
    \item{\code{NorthWembley}}{a numeric vector}
    \item{\code{Northfields}}{a numeric vector}
    \item{\code{Northolt}}{a numeric vector}
    \item{\code{NorthwickPark}}{a numeric vector}
    \item{\code{Northwood}}{a numeric vector}
    \item{\code{NorthwoodHills}}{a numeric vector}
    \item{\code{NottingHillGate}}{a numeric vector}
    \item{\code{Oakwood}}{a numeric vector}
    \item{\code{OldStreet}}{a numeric vector}
    \item{\code{Osterley}}{a numeric vector}
    \item{\code{Oval}}{a numeric vector}
    \item{\code{OxfordCircus}}{a numeric vector}
    \item{\code{Paddington}}{a numeric vector}
    \item{\code{ParkRoyal}}{a numeric vector}
    \item{\code{ParsonsGreen}}{a numeric vector}
    \item{\code{Perivale}}{a numeric vector}
    \item{\code{PiccadillyCircus}}{a numeric vector}
    \item{\code{Pimlico}}{a numeric vector}
    \item{\code{Pinner}}{a numeric vector}
    \item{\code{Plaistow}}{a numeric vector}
    \item{\code{PrestonRoad}}{a numeric vector}
    \item{\code{PutneyBridge}}{a numeric vector}
    \item{\code{QueensPark}}{a numeric vector}
    \item{\code{Queensbury}}{a numeric vector}
    \item{\code{Queensway}}{a numeric vector}
    \item{\code{RavenscourtPark}}{a numeric vector}
    \item{\code{RaynersLane}}{a numeric vector}
    \item{\code{Redbridge}}{a numeric vector}
    \item{\code{RegentsPark}}{a numeric vector}
    \item{\code{Richmond}}{a numeric vector}
    \item{\code{Rickmansworth}}{a numeric vector}
    \item{\code{RodingValley}}{a numeric vector}
    \item{\code{RoyalOak}}{a numeric vector}
    \item{\code{Ruislip}}{a numeric vector}
    \item{\code{RuislipGardens}}{a numeric vector}
    \item{\code{RuislipManor}}{a numeric vector}
    \item{\code{RussellSquare}}{a numeric vector}
    \item{\code{SevenSisters}}{a numeric vector}
    \item{\code{ShepherdsBushCen}}{a numeric vector}
    \item{\code{ShepherdsBushHC}}{a numeric vector}
    \item{\code{SloaneSquare}}{a numeric vector}
    \item{\code{Snaresbrook}}{a numeric vector}
    \item{\code{SouthEaling}}{a numeric vector}
    \item{\code{SouthHarrow}}{a numeric vector}
    \item{\code{SouthKensington}}{a numeric vector}
    \item{\code{SouthKenton}}{a numeric vector}
    \item{\code{SouthRuislip}}{a numeric vector}
    \item{\code{SouthWimbledon}}{a numeric vector}
    \item{\code{SouthWoodford}}{a numeric vector}
    \item{\code{Southfields}}{a numeric vector}
    \item{\code{Southgate}}{a numeric vector}
    \item{\code{Southwark}}{a numeric vector}
    \item{\code{StJamessPark}}{a numeric vector}
    \item{\code{StJohnsWood}}{a numeric vector}
    \item{\code{StPauls}}{a numeric vector}
    \item{\code{StamfordBrook}}{a numeric vector}
    \item{\code{Stanmore}}{a numeric vector}
    \item{\code{StepneyGreen}}{a numeric vector}
    \item{\code{Stockwell}}{a numeric vector}
    \item{\code{StonebridgePark}}{a numeric vector}
    \item{\code{Stratford}}{a numeric vector}
    \item{\code{SudburyHill}}{a numeric vector}
    \item{\code{SudburyTown}}{a numeric vector}
    \item{\code{SwissCottage}}{a numeric vector}
    \item{\code{Temple}}{a numeric vector}
    \item{\code{TheydonBois}}{a numeric vector}
    \item{\code{TootingBec}}{a numeric vector}
    \item{\code{TootingBroadway}}{a numeric vector}
    \item{\code{TottenhamCourtRoad}}{a numeric vector}
    \item{\code{TottenhamHale}}{a numeric vector}
    \item{\code{TotteridgeAndWhetstone}}{a numeric vector}
    \item{\code{TowerHill}}{a numeric vector}
    \item{\code{TufnellPark}}{a numeric vector}
    \item{\code{TurnhamGreen}}{a numeric vector}
    \item{\code{TurnpikeLane}}{a numeric vector}
    \item{\code{Upminster}}{a numeric vector}
    \item{\code{UpminsterBridge}}{a numeric vector}
    \item{\code{Upney}}{a numeric vector}
    \item{\code{UptonPark}}{a numeric vector}
    \item{\code{Uxbridge}}{a numeric vector}
    \item{\code{Vauxhall}}{a numeric vector}
    \item{\code{Victoria}}{a numeric vector}
    \item{\code{WalthamstowCentral}}{a numeric vector}
    \item{\code{Wanstead}}{a numeric vector}
    \item{\code{WarrenStreet}}{a numeric vector}
    \item{\code{WarwickAvenue}}{a numeric vector}
    \item{\code{Waterloo}}{a numeric vector}
    \item{\code{Watford}}{a numeric vector}
    \item{\code{WembleyCentral}}{a numeric vector}
    \item{\code{WembleyPark}}{a numeric vector}
    \item{\code{WestActon}}{a numeric vector}
    \item{\code{WestBrompton}}{a numeric vector}
    \item{\code{WestFinchley}}{a numeric vector}
    \item{\code{WestHam}}{a numeric vector}
    \item{\code{WestHampstead}}{a numeric vector}
    \item{\code{WestHarrow}}{a numeric vector}
    \item{\code{WestKensington}}{a numeric vector}
    \item{\code{WestRuislip}}{a numeric vector}
    \item{\code{WestbournePark}}{a numeric vector}
    \item{\code{Westminster}}{a numeric vector}
    \item{\code{WhiteCity}}{a numeric vector}
    \item{\code{Whitechapel}}{a numeric vector}
    \item{\code{WillesdenGreen}}{a numeric vector}
    \item{\code{WillesdenJunction}}{a numeric vector}
    \item{\code{Wimbledon}}{a numeric vector}
    \item{\code{WimbledonPark}}{a numeric vector}
    \item{\code{WoodGreen}}{a numeric vector}
    \item{\code{WoodLane}}{a numeric vector}
    \item{\code{Woodford}}{a numeric vector}
    \item{\code{WoodsidePark}}{a numeric vector}
    \item{\code{Total}}{a numeric vector; the total over all stations.}
    \item{\code{mins24}}{a numeric vector; minutes on a 24 hour clock,
       e.g., 0 is midnight, 120 is 2am. }
  }
}
\details{
  Each cell contains the number of passengers entering a station
  during a quarter-hour period of time on a 
  typical day during November 2010.
  The column names of the data frame are the station names and
  the most of the rows are 
  the start time of each quarter-hour time block given in 24 hour time,
  e.g., 2215 means 10:15pm to 10:29pm. 
  The last four rows are
  \code{"Total"},
  \code{"A.M.Peak"},
  \code{"Interpeak"},
  \code{"P.M.Peak"}.
%  which are self-explanatory.
 


  The data is adjusted to remove the effect of abnormal
  circumstances that many affect passeger numbers such 
  as industrial action.



}
\source{
  The data comes from the UK Government Transport for London website
  \url{http://www.tfl.gov.uk/info-for/open-data-users/our-feeds?intcmp=3671}.
  Downloaded in December 2013 and formatted into R by J. T. Gray
  (and slightly edited by T. W. Yee).


}
%%\references{
%%}
\examples{
\dontrun{
data(tube10)
waterloo <- tube10[1:(4*24), "Waterloo"]  # Omit the totals and the peaks
barplot(unlist(waterloo))
barplot(log10(1 + unlist(waterloo)), col = "lightblue",
        ylab = "log10(1+.)", las = 1)
}
}
\keyword{datasets}



% ## maybe str(tube10) ; plot(tube10) ...
% summary(tube10)
% N <- nrow(tube10)
% waterloo <-tube10[-((N-3):N), "Waterloo"]

