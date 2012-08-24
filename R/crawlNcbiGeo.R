
crawlNcbiGeo <- function() {
	# Gather information on available data from the repository
	ncbi.input.file <- paste(file.path(.path.package("rSCR"),"extdata"),"/ncbiGPLIDs.txt",sep="")
	ncbi.perl.crawler  <- paste(file.path(.path.package("rSCR"),"Perl"),"/getGSEs.pl",sep="")
	system(paste("perl",ncbi.perl.crawler,ncbi.input.file))
	# Translate the information so we can contribute it to Synapse
	all.gses <- read.table("all.GSEs.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)
	names(all.gses) <- c('data.lastUpdate', 'species', 'description', 'data.url', 'numSamples', 'platform')
	all.gses[is.na(all.gses)] <- ""
	all.gses[,"study.name"] <- rownames(all.gses)
	all.gses[,"data.type"] <- "E"
	all.gses[,"data.status"] <- "raw"
	all.gses[,"repository"] <- "NCBI GEO"
  all.gses[,"data.name"] <- paste(all.gses[,"study.name"], "Raw Data Layer from NCBI GEO")
	all.gses[,"numSamples"] <- as.numeric(all.gses[,"numSamples"])
#	all.gses <- all.gses[which(all.gses$rawDataAvailable == "TRUE"),]
#	all.gses <- all.gses[,which(colnames(all.gses) != "rawDataAvailable")]
	
	return(all.gses)
}

