
crawlArrayExpress <- function() {
# Load valid array express platforms.
	aeInputFile <- paste(file.path(.path.package("rSCR"),"extdata"),"/arrayExpressPlatforms.txt",sep="")
	aePlatformMap <- read.table(aeInputFile,stringsAsFactors=FALSE, row.names=1, sep="\t", header=TRUE)
	
	# Access the array express database through the JSON API, then parse the object.
	arrayExpressJSON <- .getArrayExpressJSON()
	parsedArrayExpressJSON <- .parseAEJsonObject(arrayExpressJSON,aePlatformMap)
	
	# Now we get all the URLs for the raw data layers of each study
	# Unfortunately an indeterminate # of raw data layers exist for each dataset.
	# So, we are forced to figure out by brute force.
	arrayExpressURLs <- .getArrayExpressRawDataURLs(parsedArrayExpressJSON)
	
	# Now create the arrayExpress matrix that we want to write to a file.
	arrayExpress <- list()
	total.urls <- sum(sapply(arrayExpressURLs,length))
	arrayExpress <- vector("list", total.urls)
	counter <- 1	
	for(i in 1:length(arrayExpressURLs)){
		x <- arrayExpressURLs[[i]]
		dataset <- names(arrayExpressURLs[i])
		for(j in 1:length(x)){
			arrayExpress[counter] <- parsedArrayExpressJSON[dataset]
			arrayExpress[[counter]]$data.url <- x[j]
			arrayExpress[[counter]]$data.name <- paste(parsedArrayExpressJSON[[dataset]]$dataset.name,"raw data layer",j)
			counter <- counter+1
		}
	}
	arrayExpress <- t(sapply(arrayExpress, function(x){ unlist(x)}))
	colnames(arrayExpress) <- gsub("layer","data",colnames(arrayExpress))
	colnames(arrayExpress) <- gsub("dataset","study",colnames(arrayExpress))
	arrayExpress <- arrayExpress[,-2]
	arrayExpress[,"description"] <- gsub("\n","",arrayExpress[,"description"])
	arrayExpress <- cbind(arrayExpress, rep("Array Express",nrow(arrayExpress)))
	colnames(arrayExpress)[ncol(arrayExpress)] <- "repository"
	return(arrayExpress)
}

.getArrayExpressRawDataURLs <- function(output) {
	urls <- vector("list", length(output))
	names(urls) <- names(output)
	for(i in 1:length(output)){
		cat("\r",i, "\t", names(output)[i])
		obj <- output[[i]]
		if(obj$rawDataAvailable  == FALSE){
			urls[[i]] <- "Not Available"
		}else{
			url <- try(getURL(paste("http://www.ebi.ac.uk/arrayexpress/json/v2/files/", obj$dataset.name, sep = "")))
			if(class(url)=="try-error"){
				urls[[i]] <- 'Not Available'
				next;
			}
			files <- fromJSON(url)
			# If multiple file exist, find the slot that corresponds to just these samples.
			experiment <- NA
			if(files$files$"total-experiments" > 1){ 
				experiment <- files$files$experiment[[which(sapply(files$files$experiment, function(x){ x$accession}) == obj$dataset.name)]]
			}else{
				experiment <- files$files$experiment
			}
			urls[[i]] <- sapply(experiment$file, function(x) {
						if (!is.null(x$kind)) {
							if (x$kind == "raw") {
								x$url
							}
							else {
								NA
							}
						}
						else {
							NA
						}
					})
			urls[[i]] <- urls[[i]][which(!is.na(urls[[i]]))]
		}
	}
	
return(urls)
}

.getProvider <- function(obj){
	accession <- sapply(obj, function(x){ x$accession})
	sapply(obj, function(x){
				contact=''
				if(length(x$provider) == 3){ 
					contact <- paste(x$provider,collapse=",")
				}else{
					prov <- grep("submitter", x$provider)
					inve <- grep("submitter", x$investigator)
					ret <- NA
					if(length(prov) > 0 & length(inve) > 0) {
						contact <- paste(paste(x$provider[prov],collapse=","),
								paste(x$provider[inve],collapse=","))
					}else if(length(prov) > 0){
						contact <- paste(x$provider[prov],collapse=",")
					}else if(length(inve) > 0){
						contact <- paste(x$provider[inve],collapse=",")
					}
				}
				contact
			}) -> res
	names(res) <- accession	
}

.parseAEJsonObject <- function(arrayExpressJSON,aePlatformMap){
	
	sort(unique(unlist(sapply(arrayExpressJSON, function(x){ names(x)}))))
	[1] "accession"          "arraydesign"        "assays"            
	[4] "bibliography"       "bioassaydatagroup"  "description"       
	[7] "experimentalfactor" "experimentdesign"   "experimenttype"    
	[10] "fgemdatafiles"      "id"                 "lastupdatedate"    
	[13] "miamescores"        "minseqescores"      "name"              
	[16] "protocol"           "provider"           "rawdatafiles"      
	[19] "releasedate"        "sampleattribute"    "samples"           
	[22] "secondaryaccession" "species"            "submissiondate"    
	
	
	# Get the study name, url, description, 
	species <- sapply(arrayExpressJSON, function(x){ paste(unlist(x$species), collapse=",")})
	accession <- sapply(arrayExpressJSON, function(x){ x$accession })
	samples <- sapply(arrayExpressJSON, function(x){ x$samples})
	arraydesign <- sapply(arrayExpressJSON, function(x){ x$arraydesign })
	description <- sapply(arrayExpressJSON, function(x){ paste(unlist(x$description), collapse=" ")}) 
	sapply(list(species, accession, samples, arraydesign, description), length)
	
	rawDataAvailable <- unlist(sapply(obj, 
					function(x){ x$rawdatafiles}))
	not.available <- which(! rawDataAvailable )
	available <- which( rawDataAvailable )
	accession <- sapply(obj, function(x){ x$accession})
	
	sapply(obj, function(x){ 
				if(length(x)>1 & is.null(names(x$arraydesign))){ 
					paste(sapply(x$arraydesign, function(y){ y$accession}),collapse=",")
				}else{
					x$arraydesign$accession
				}
			}) -> arrayDesignAccession
	
	sapply(obj, function(x){ 
				if(length(x$arraydesign)>1 & is.null(names(x$arraydesign))){ 
					vec <- sapply(x$arraydesign, function(y){ y$accession})
				}else{
					vec <- x$arraydesign$accession
				}
				ids <- which(vec %in% aE.platformMap[,4])
				plats <- vec
				plats[ids] <- aE.platformMap[match(vec[ids],aE.platformMap[,4]),2]
				if(length(plats) > 1){ 
					paste(plats,collapse=", ")
				}else{
					plats
				}
			}) -> platforms
	
	lastupdate <- sapply(obj, function(x){ x$lastupdatedate})
	
	sapply(obj, function(x){ 
				if(length(x$species)>1){ 
					paste(sapply(x$species, function(y){ y}),collapse=",")
				}else{
					x$species
				}
			}) -> species
	
	study <- sapply(obj, function(x){ x$name})
	
	rawdata <- unlist(sapply(obj, function(x){ x$rawdatafiles}))
	
	sapply(obj, function(x){ 
				if(is.list(x$samples)){ 
					sapply(x$samples, function(y) { 
								y$samples
							})
				}else{
					x$samples
				}
			}) -> samples
	
	emails <- sapply(obj, function(x){ 
				if(is.list(x)){
					x$provider[[1]]["email"]
				}else{
					x$provider["email"]
				}})
	
	sapply(obj, function(x){ 
				if(length(x$arraydesign)==2){ 
					paste(sapply(x$arraydesign, function(y){ gsub("\\]","",gsub("[^\\[]+\\[","",y$name,perl=TRUE),perl=TRUE)}),collapse=",")
				}else{
					gsub("\\]","",gsub("[^\\[]+\\[","",x$arraydesign$name,perl=TRUE),perl=TRUE)
					
				}
			}) -> arrayDesignName
	
#	arrayExpress <- matrix("NA",nr=length(obj), nc=8)
#	for( i in 1:length(obj)) {
#		arrayExpress[i,1] <- accession[i]  
#		arrayExpress[i,2] <- ifelse(is.null(arrayDesignAccession[[i]]),"NA", arrayDesignAccession[[i]])
#		arrayExpress[i,3] <- ifelse(is.null(lastupdate[[i]]),"NA", lastupdate[[i]])
#		arrayExpress[i,4] <- species[[i]]
#		arrayExpress[i,5] <- study[[i]][1]
#		arrayExpress[i,6] <- rawdata[[i]]
#		arrayExpress[i,7] <- samples[[i]]
#		arrayExpress[i,8] <- ifelse(is.null(emails[[i]][[1]]),"NA", emails[[i]])
#		arrayExpress[i,8] <- ifelse(is.null(arrayDesignName[[i]]),"NA", arrayDesignName[[i]])
#	}
	
	output <- list()
	lapply(1:length(obj), function(i) {  
				b <- list(
						dataset.name = accession[i],
						layer.url= paste("http://www.ebi.ac.uk/arrayexpress/experiments/",accession[i],sep=""),
						lastUpdate = ifelse(is.null(lastupdate[[i]]),"NA", lastupdate[[i]]),
						species = species[[i]],
						description = study[[i]],
						rawDataAvailable = rawdata[[i]],
						dataset.numSamples = samples[[i]], 
						contact = ifelse(is.null(emails[[i]]),"Not Available", emails[[i]]),
						platform = platforms[i],
						layer.status="raw",
						layer.type="E")
				b
			}) -> output
	names(output) <- accession
	lapply(output, function(x){ 
				x[which(is.na(x))] <- "Not Available"
				x
			}) -> output
	output
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
