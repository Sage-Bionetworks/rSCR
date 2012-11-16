
crawlArrayExpress <- function() {
	# Function crawls the array express repository.
	# Returns a list, each element is a valid study.
	
	# Load valid array express platforms.
	aeInputFile <- paste(file.path(.path.package("rSCR"),"extdata"),"/arrayExpressPlatforms.txt",sep="")
	aePlatformMap <- read.table(aeInputFile,stringsAsFactors=FALSE, row.names=1, sep="\t", header=TRUE)
	
	# Access the array express database through the JSON API, then parse the object.
	arrayExpressJSON <- .getArrayExpressJSON()
	parsedArrayExpressJSON <- .parseAEJsonObject(arrayExpressJSON,aePlatformMap)
	
	# Now we add the URLs for the raw data layers of each study
	arrayExpress <- .getArrayExpressRawDataURLs(parsedArrayExpressJSON)
	arrayExpress <- .storeArrayExpressInSynapse(arrayExpress)
	# Return the list. 
	return(arrayExpress)
}

.storeArrayExpressInSynapse <- function(arrayExpress){ 
	# Adds the output of the GEO crawler to Synapse
	entName <- 'arrayExpressPublicCrawlerOutput'
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	if(!is.null(qry)){
		arrayExpressCrawlerEntity <- loadEntity(qry$entity.id)
	}else{
		arrayExpressCrawlerEntity <- Data(list(name = entName, parentId = 'syn1452692'))
	}
	arrayExpressCrawlerEntity <- addObject(arrayExpressCrawlerEntity,arrayExpress)
	if(is.null(qry)){
		arrayExpressCrawlerEntity <- createEntity(arrayExpressCrawlerEntity)		
	}else{
		arrayExpressCrawlerEntity <- updateEntity(arrayExpressCrawlerEntity)
	}
	arrayExpressCrawlerEntity <- storeEntity(arrayExpressCrawlerEntity)
	return(arrayExpressCrawlerEntity)
}

.getArrayExpressRawDataURLs <- function(output) {
	# Function gets the URLs for the available raw data for each study. 
	sapply(1:nrow(output), function(i){
				cat("\r",i, "\t", output$accession[i])
				if(output$rawDataAvailable[i]  == FALSE){
					urls[[i]] <- "Not Available"
					return(NA)
				}else{
					url <- try(getURL(paste("http://www.ebi.ac.uk/arrayexpress/experiments/", output$accession[i], sep = "")))
					rawDataLinks <- regmatches(url,
							gregexpr('<a href=\\\"[^"]+raw\\.\\d+\\.zip\\\">[^<]+<\\/a>',url))[[1]]
					baseUrl <- 'http://www.ebi.ac.uk/arrayexpress/files/'
					
					links <- unlist(regmatches(rawDataLinks,
									gregexpr('<a href=\\\"[^"]+\\\">[^<]+<',rawDataLinks)))
					sapply(strsplit(rawDataLinks, '>'), function(x){ 
								thisUrl <- paste(baseUrl, output$accession[i], '/', gsub("<\\/a", "", x[[2]], perl=TRUE), sep="")
							}) -> allUrls
					return(allUrls)
				}
			}) -> urls
	res <- list()
	for(i in 1:length(urls)){
		a <- list(name=output$accession[i],
				species=output$species[i],
				numSamples=output$samples[i],
				description=output$description[i],
				platform=output$platform[i],
				urls = urls[[i]])
		res[[i]] <- a
	}
	names(res) <- output$accession
	return(res)
}



.translatePlatform <- function(map, platform){
	# Simply maps the array express internal id to our common platform name
	if(platform %in% names(map)){
		return(as.character(unlist(map[platform])))
	}else{
		return(platform)
	}
}

.parseAEJsonObject <- function(arrayExpressJSON,aePlatformMap){
	# Build map from array express platform identifiers to SCR platform names
	aeIdentifierToString <- split(aePlatformMap$String, aePlatformMap$GPL_ID)
	
	# Get the study name, url, description, 
	species <- sapply(arrayExpressJSON, function(x){ paste(unlist(x$species), collapse=",")})
	accession <- sapply(arrayExpressJSON, function(x){ x$accession })
	samples <- sapply(arrayExpressJSON, function(x){ x$samples})
	arraydesign <- sapply(arrayExpressJSON, function(x){ x$arraydesign })
	description <- sapply(arrayExpressJSON, function(x){ paste(unlist(x$description), collapse=" ")}) 
	rawdatafiles <- sapply(arrayExpressJSON, function(x){ x$rawdatafiles })
	# Parse out the multiple platforms 
	sapply(1:length(arraydesign), function(i){
				x <- arraydesign[[i]]
				lngth <- length(unlist(x[[1]]))
				if(lngth == 1){
					x <- unlist(x)
					nms <- names(x)
					id <- which(nms=="accession")
					return(.translatePlatform(aeIdentifierToString, x[id]))
				}else{
					# Multiple levels. Return each one
					sapply(x, function(y){
								y <- unlist(y)
								nms <- names(y)
								id <- which(nms=="accession")
								.translatePlatform(aeIdentifierToString, y[id])
							}) -> z
					paste(unlist(z), collapse=",")
				}
			}) -> arrayAccessions
	
	# Join into a single data frame
	obj <- data.frame(accession = accession,
			species=species, 
			samples=samples,
			description=description,
			platform=arrayAccessions, 
			rawDataAvailable = rawdatafiles, stringsAsFactors=FALSE)
	obj <- obj[which(obj$rawDataAvailable == "TRUE"),]
	return(obj)
}

.getArrayExpressJSON <- function(){
	# The following downloads the JSON object
	ae.json <- getURL("http://www.ebi.ac.uk/arrayexpress/json/v2/experiments?directsub=true&array=A-AFFY*")
	obj.affy <- fromJSON(ae.json)
	ae.json <- getURL("http://www.ebi.ac.uk/arrayexpress/json/v2/experiments?directsub=true&array=A-GEOD*")
	obj.geod <- fromJSON(ae.json)
	ae.json <- getURL("http://www.ebi.ac.uk/arrayexpress/json/v2/experiments?directsub=true&array=A-AGIL*")
	obj.agil<- fromJSON(ae.json)
	return(c(obj.affy[[1]][[6]], obj.geod[[1]][[6]], obj.agil[[1]][[6]]))
}
