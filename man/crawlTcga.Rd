\name{crawlTcga}
\alias{crawlTcga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function calls the public or private TCGA http data portals.
}
\description{
This function crawls the TCGA public and private data portals.  It returns a data frame in which
each row is a unique data entity that should be made available in either the public or private
SCR.  
}
\usage{
crawlTcga(user = "anonymous", pwd = "anonymous", private = FALSE, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{user}{
Username required to crawl the private portal.
}
  \item{pwd}{
Password required to crawl the private portal.
}
  \item{private}{
Logical value.  If TRUE then it will crawl the private portal.
}
  \item{writeOut}{
Logical value. If true the data frame is written out to a file in the current working directory named
tcga.txt.
}
  \item{debug}{
Logical value.  If TRUE then it will only crawl the GBM data.  Used only for testing.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (debug = FALSE, user = "anonymous", pwd = "anonymous", 
    private = FALSE) 
{
    header <- c("study.name\tacronym\ttissue\trepository\tdisease\tspecies\tstudy.numPatients\tdata.url\tdata.name\tdata.md5\tdata.type\tdata.platformName\tdata.status\tdata.add\tdata.lastUpdate\tdata.numSamples\tdata.tcgaLevel")
    if (isTRUE(private)) {
        url <- "https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/"
        write.table(header, file = "tcgaPrivate.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    else {
        url <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
        write.table(header, file = "tcga.txt", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, append = FALSE)
    }
    links <- .getLinks(url, debug = debug, user = user, pwd = pwd)
    tcga <- read.table("tcga.txt", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    return(tcga)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
