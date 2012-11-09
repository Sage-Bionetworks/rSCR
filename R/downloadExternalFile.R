downloadExternalFile <- function(entity, url, md5, maxFileSize=(20 * 1073741824), numRetries=20, store=FALSE, private=FALSE, user='anonymous', pwd='anonymous'){
	# Check if the request if on a private url.  If so, add user, pwd to url
	if(private=="TRUE"){
		url <- gsub("https://", paste("https://", user,':', pwd, '@', sep=""), url)
	}
	# This function adds the file at the specified url.
	parsedURL <- .parseURL(url)
	# Does file already exist in local cache?
	if(file.exists(parsedURL)[1]){
		# File exists, so we add file to the entity and return it
		cat("File exists, adding file\n")
		fileSize <- file.info(parsedURL[1])$size
		annotValue(entity,'fileSize') <- .prettifyFileSize(fileSize)
		if(private != "TRUE"){
			entity <- addFile(entity, parsedURL[1])
			entity <- storeEntity(entity)
		}else{
			entity <- createEntity(entity)
			propertyValue(entity, 'md5') <- as.character(tools::md5sum(parsedURL[1]))
			propertyValue(entity, 'locations') <- list(list(path=tcga$data.url[i], type="external"))
			entity <- updateEntity(entity)
		}
		return(entity)
	}
	cat("File does not exist.  Downloading...\n")
	#If not, then try to download.
	fileSize <- .getFileSize(url)
	if(fileSize > maxFileSize){
		# File is more than 20GB large.
		var <- paste("File too large.  Will not download files larger than ", .prettifyFileSize(maxFileSize),".\n");
		class(var) <- "try-error"
		return(var)
	}
	if(!file.exists(parsedURL[2])){
		dir.create(parsedURL[2], recursive = TRUE) # Create local directory if it doesn't already exist
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
	file.link(destfile, parsedURL[1])	
	annotValue(entity,'fileSize') <- .prettifyFileSize(fileSize)
	if(private != "TRUE"){
		entity <- addFile(entity, parsedURL[1])
		entity <- storeEntity(entity)
	}else{
		entity <- createEntity(entity)
		propertyValue(entity, 'md5') <- as.character(tools::md5sum(parsedURL[1]))
		propertyValue(entity, 'locations') <- list(list(path=tcga$data.url[i], type="external"))
		entity <- updateEntity(entity)
	}
	return(entity)
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

.parseURL <- function(url){ 
	# Returns vector. Element 1 is the path to file in the synapse cache. 
	# Element 2 is its directory.
	# Useful to determine if a file already has been downloaded
	parsedUrl <- synapseClient:::.ParsedUrl(url)
	destfile <- file.path(synapseCacheDir(), gsub("^/", "", parsedUrl@path))
	destfile <- path.expand(destfile)
	destdir <- gsub(parsedUrl@file, "", destfile, fixed = TRUE)
	destdir <- gsub("[\\/]+$", "", destdir)			
	c(destfile, destdir)
}