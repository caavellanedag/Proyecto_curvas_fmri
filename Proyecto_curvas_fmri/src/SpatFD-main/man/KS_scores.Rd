\name{KS_scores}
\alias{KS_scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fuctional kriging using simple kriging
     %%  ~~function to do ... ~~
}
\description{
        Fuctional kriging using simple kriging
     %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
     KS_scores(SFD, newcoords,model,vari=NULL,fill.all=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
     \item{SFD}{Is an object of class SpatFD.}
     \item{newcoords}{A data fame with two columns indicating X and Y coordinates.}
     \item{model}{An object of class variogramModel or a list of objects of class variogramModel(see vgm from gstat package).}
     \item{vari}{A numeric value to specifying the position of the variable or a character value specifying the name of the variable. NULL by default.}
     \item{fill.all}{A boolean value to fit the same model to all scores. TRUE by default. If TRUE, model must be class variogramModel. If FALSE, model must be a list of objects of class variogramModel.}
}

\value{
     An object of class scores_pred, which inherits from list. Its components are:
          \item{X}{the object of class scores_pred where are specified the functional predictions.}

     %%  ~Describe the value returned
     %%  If it is a LIST, use
     %%  \item{comp1 }{Description of 'comp1'}
     %%  \item{comp2 }{Description of 'comp2'}
     %% ...
}

\author{ Angie Villamil, Martha Bohorquez
     %%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

     \seealso{ \code{\link{SpatFD}}
     %%  \code{\link{FKCK}}, \code{\link{FCOK}},
          %% ~~objects to See Also as \code{\link{help}}, ~~~
     }
\examples{
     ##---- Should be DIRECTLY executable !! ----
     ##-- ==>  Define data, use random,
     ##--	or do  help(data=index)  for the standard data sets.

     ## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Functional Kriging}% use one of  RShowDoc("KEYWORDS")
