
updateArrayExpress <- function(id){ 
# run and return output from the GEO crawler. 
	aeCrawlerOutput <- loadEntity('syn1488304');
	ae <- aeCrawlerOutput$objects$arrayExpress
	
# Filter out the ones we've already created.
	alreadyCreatedFolders <- synapseQuery('select id, name from entity where entity.parentId=="syn1488308"')
	if(!is.null(alreadyCreatedFolders )){
		ae <- ae[setdiff(names(ae), alreadyCreatedFolders$entity.name)]
	}
	
# Anything that's new should be added
	for(i in 1:length(ae)){
		cat("\n\n", i)
		res <- try(.contributeArrayExpressStudy(ae, i), silent=TRUE)
		numRetries <- 1;
		while(class(res) == "try-error") {
			res <- try(.contributeArrayExpressStudy(ae, i), silent=TRUE)
			numRetries <- numRetries + 1
			if(numRetries > 4){ 
				#stop("Tried 10 times to build the entity with no success!")
				break;
			}
		}
	}
}

.contributeArrayExpressStudy <- function(ae, i, maxFileSize=(1 * 1073741824)){
	# For each ae study
	# Create folder with content from crawler
	folder <- Folder(list(parentId='syn1488308',
					name=names(ae)[i],
					description = ae[[i]]$description))
	annotValue(folder, 'numSamples') = as.numeric(ae[[i]]$numSamples)
	annotValue(folder, 'platform') = strsplit(ae[[i]]$platform, ',')[[1]]
	annotValue(folder, 'species') = strsplit(ae[[i]]$species, ',')[[1]]
	annotValue(folder, 'repository') <- 'Array Express'
	folder <- createEntity(folder)
	
	# Create the raw folder
	rawFolder <- Folder(list(parentId=propertyValue(folder, 'id'),
					name='raw data',
					description = ae[[i]]$description))
	annotValue(rawFolder, 'numSamples') = as.numeric(ae[[i]]$numSamples)
	annotValue(rawFolder, 'platform') = strsplit(ae[[i]]$platform, ',')[[1]]
	annotValue(rawFolder, 'species') = strsplit(ae[[i]]$species, ',')[[1]]
	annotValue(rawFolder, 'repository') <- 'Array Express'
	rawFolder <- createEntity(rawFolder)
	
	urls <- unlist(strsplit(ae[[i]]$urls, ','))
	for(j in 1:length(urls)){
		# Create a raw data entity for each of these external links
		rawDataEntity <- Data(list(name=basename(urls[j]), 
						parentId=propertyValue(rawFolder, 'id'),
						platform = strsplit(ae[[i]]$platform, ',')[[1]],
						species = strsplit(ae[[i]]$species, ',')[[1]]))
		annotValue(rawDataEntity, 'repository') <- 'NCBI GEO'
		fileSizeInBytes <- .getFileSize(urls[j])
		annotValue(rawDataEntity, 'fileSize') <- .prettifyFileSize(fileSizeInBytes)
		rawDataEntity <- .addExternalLocationToAeRawDataEntity(rawDataEntity, urls[j])
		rawDataEntity <- createEntity(rawDataEntity)
		cat(propertyValue(rawDataEntity,'id'), "\n")
	}
	return(folder)
}

.addExternalLocationToAeRawDataEntity <- function(entity, url, maxFileSize=(5 * 1073741824), numRetries=20){
	fileSize <- .getFileSize(url)
	if(fileSize > maxFileSize){
		# File is more than 20GB large.
		var <- paste("File too large.  Will not download files larger than ", .prettifyFileSize(maxFileSize),".\n");
		class(var) <- "try-error"
		return(var)
	}
	destfile =  tempfile()
	for(i in 1:numRetries) { # Try numRetries times to download.    
		download <- try(synapseClient:::.curlWriterDownload(url,destfile))
		if(class(download) != "try-error"){
			break; # download successful.  
		}else{
			file.remove(destfile) # If it fails, remove whatever file was there and try again
		}
	}
	if(i == numRetries & class(download) != "try-error"){ # Throw error if the download breaks.
		var <- "Could not download data\n"
		class(var) <- "try-error"
		return(var)
	}
	# Download successful.  Add it and return the entity	
	propertyValue(entity, 'md5') <- as.character(tools::md5sum(destfile))
	propertyValue(entity, 'locations') <- list(list(path=url, type="external"))
	file.remove(destfile)
	entity
}

.prettifyFileSize <- function(fileSize){
	if(fileSize > 1073741824){
		# Its bigger than a gigabyte, so translate to something GB
		return(paste(round(fileSize / 1073741824,1), 'GB'))
	}else{
		size <- round(fileSize / 1048576,1)
		if(size < 0.1){ 
			size <- 0.1
		}
		size <- paste(size, 'MB',sep=" ")
		return(size)
	}
}

.getFileSize <- function(url){ 
	curlHandle = getCurlHandle() 
	h = basicTextGatherer()
	response <- curlPerform(URL=url, 
			.opts=list(header=TRUE, nobody=TRUE),  
			writefunction=h$update, 
			curl=curlHandle) 
	curlInfo <- getCurlInfo(curlHandle)
	return(curlInfo$content.length.download)
}

