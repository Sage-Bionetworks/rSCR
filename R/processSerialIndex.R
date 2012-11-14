processSerialIndex <- function(id, user, pwd){
	# Download the synapse id
	.downloadEntity(id, user, pwd)
	
	
}

.downloadEntity <- function(id, user, pwd, numRetries=10) {
	entity <- getEntity(id)
	private <- annotValue(entity, 'private')
	url <- propertyValue(entity, 'locations')[[1]][1]
	if(private=="TRUE"){
		url <- gsub("https://", paste("https://", user,':', pwd, '@', sep=""), url)
	}
	parsedURL <- .parseURL(url)
	if(file.exists(parsedURL[1])){
		# File exists locally
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
	
}