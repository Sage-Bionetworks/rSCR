\name{updateScrWithTcgaCrawlerOutput}
\alias{updateScrWithTcgaCrawlerOutput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Accepts output from TCGA crawler and adds any new data entities to SCR public or private. 
}
\description{
This function builds new TCGA data entities into the SCR.  
}
\usage{
updateScrWithTcgaCrawlerOutput(tcgaCrawlerOutput, private = FALSE, update = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tcgaCrawlerOutput}{
Data frame containing output from crawlTcga() function.
}
  \item{private}{
Logical value.  If TRUE then new data entities will be added to the private project, otherwise they will be added to the public project. 
}
  \item{update}{
Logical value.  If TRUE then the data entities are contributed using the contribute data function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The newly created entities. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#tcgaCrawlerOutput <- crawlTcga()
#newEntities <- contributeData(tcgaCrawlerOutput)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
