
crawlNcbiGeo <- function() {
	# Gather information on available data from the repository
	ncbi.input.file <- paste(file.path(.path.package("rSCR"),"extdata"),"/ncbiGPLIDs.txt",sep="")
	ncbi.perl.crawler  <- paste(file.path(.path.package("rSCR"),"Perl"),"/getGSEs.pl",sep="")
	system(paste("perl",ncbi.perl.crawler,ncbi.input.file))
	# Translate the information so we can contribute it to Synapse
	geo <- read.table("geo.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE,row.names=1)
	geo <- storeGeoInSynapse(geo)
	return(geo)
}

.storeGeoInSynapse <- function(geo, private=FALSE, fname='geo.txt'){ 
	# Adds the output of the GEO crawler to Synapse
	entName <- 'geoPublicCrawlerOutput'
	if(private=="TRUE"){
		entName <- 'geoPrivateCrawlerOutput'
	}
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	if(!is.null(qry)){
		geoCrawlerEntity <- getEntity(qry$entity.id)
	}else{
		geoCrawlerEntity <- Data(list(name = entName, parentId = 'syn1452692'))
	}
	geoCrawlerEntity <- addFile(geoCrawlerEntity,fname)
	geoCrawlerEntity <- addObject(geoCrawlerEntity,geo)
	if(is.null(qry)){
		geoCrawlerEntity <- createEntity(geoCrawlerEntity)		
	}else{
		geoCrawlerEntity <- updateEntity(geoCrawlerEntity)
	}
	geoCrawlerEntity <- storeEntity(geoCrawlerEntity)
	return(geoCrawlerEntity)
}
