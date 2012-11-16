\name{crawlNcbiGeo}
\alias{crawlNcbiGeo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crawls the NCBI GEO repository.
}
\description{
Crawls the NCBI GEO repository to identify available data.  Stores
the output as an entity in Synapse.
}
\usage{
crawlNcbiGeo()
}
%- maybe also 'usage' for other objects documented here.
\value{
Returns a synapse entity containing the output of the crawler.

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
function () 
{
    ncbi.input.file <- paste(file.path(.path.package("rSCR"), 
        "extdata"), "/ncbiGPLIDs.txt", sep = "")
    ncbi.perl.crawler <- paste(file.path(.path.package("rSCR"), 
        "Perl"), "/getGSEs.pl", sep = "")
    system(paste("perl", ncbi.perl.crawler, ncbi.input.file))
    geo <- read.table("all.GSEs.txt", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE, row.names = 1)
    geo <- geo[which(geo$layer.url != "FALSE"), ]
    geo <- .storeGeoInSynapse(geo)
    return(geo)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
