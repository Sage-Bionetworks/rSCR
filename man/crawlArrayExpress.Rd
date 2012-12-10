\name{crawlArrayExpress}
\alias{crawlArrayExpress}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Crawler for Array Express Repository.
}
\description{
This function crawls the array express repository via their API.
It accesses their JSON web services to find studies that were directly
submitted to Array Express and use either Affymetrix or Agilent platforms.
It returns a list containing information about available studies.
}
\usage{
crawlArrayExpress()
}
%- maybe also 'usage' for other objects documented here.
\value{
A synapse entity containing a list of available studies. 
The list has the following 6 items for each study.

\item{name} The study name.

\item{species} The species explored by the study.

\item{numSamples} The number of samples

\item{description} A description for the study.

\item{platform} The platform(s) used by the experiment.  Mapped to the SCR platform names.

\item{urls} List of URLs for the raw data.
}
\author{
Brig Mecham <brig.mecham@sagebase.org>
}

\examples{
% arrayExpress <- crawlArrayExpress()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
