\name{contributeData}
\alias{contributeData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
contributeData
}
\description{
}
Main function used to facilitate data contributions to the Synapse Commons Repository.  More information on
this project  is available here:

https://sagebionetworks.jira.com/wiki/display/SCR/Home

\usage{
contributeData(contribution, project = NULL, keepLocal = FALSE, logFile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contribution}{
Either a list containing properties or annotations for a single study / dataset, a dataframe with each row consisting of a single data entity to 
contribute, or a character string corresponding to a file name.  Within the file should be a matrix with a header row, with each subsequent row consisting of
a study / data entity to contribute.
}
  \item{project}{
The synapse ID corresponding to the project the user wishes to contribute data to.
}
  \item{keepLocal}{
A logical flag that determines if the data file will be deleted from the local file system once the entity has been created in Synapse.
}
  \item{logFile}{
A character string corresponding to a filename or set equal to TRUE in which case the log file will be set to a random file name.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Either an entity or list of entities corresponding to the contributed data in synapse.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brig Mecham <brig.mecham@sagebase.org>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (contribution, project = NULL, keepLocal = FALSE, logFile = TRUE) 
{
    if (logFile == "TRUE") {
        logFile <- .setUpLogFile(logFile)
        cat("log file for this workflow: ", logFile, "\n")
    }
    project <- handleProject(project)
    translatedContribution <- handleContribution(project, contribution, 
        logFile)
    if (class(translatedContribution[[1]]) == "character") {
        layer <- contribute(translatedContribution, keepLocal, 
            logFile)
    }
    else {
        layer <- sapply(translatedContribution, function(x) {
            contribute(x, keepLocal, logFile)
        })
    }
    return(layer)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
