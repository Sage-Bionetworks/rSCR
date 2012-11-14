downloadPrivateTcgaFiles <- function(entity, user=user,pwd=pwd){
	fname2 <- propertyValue(entity, 'name')
	cat(fname,"\n")
	url <- propertyValue(entity, 'locations')[[1]][1]
	fullPath <- paste(synapseCacheDir(), '/', paste(strsplit(url,'/')[[1]][c(-1,-2,-3)],collapse="/"),sep="")
	if(isTRUE(file.exists(fullPath))){
		# File already exists in synapse cache
	}else{
		# File does not exist yet.  Download it and move it to the appropriate place
		cat("Downloading ", fname)
		string <- gsub("https://", paste("https://", user,':', pwd, '@', sep=""), url)
		download <- try(synapseClient:::.curlWriterDownload(string,fullPath))
	}
	entity <- loadEntity(propertyValue(entity,'id'))
}
