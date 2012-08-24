downloadPrivateTcgaFiles <- function(user=user,pwd=pwd,entitiesToDownload){
	lEnt <- entitiesToDownload[grep("SNP",entitiesToDownload$entity.name),]
	urls <- c(0)
	for(i in 1:nrow(lEnt)) {	
		cat(i,"\n")
		entity <- getEntity(lEnt[i,"entity.id"])
		fname <- paste(annotValue(entity,"acronym"),'/',lEnt[i,'entity.name'],sep="")
		fname2 <- as.character(unlist(lEnt$entity.name[i]))
		cat(fname,"\n")
		url <- propertyValue(entity, 'locations')[[1]]$path
		fullPath <- paste(synapseCacheDir(), '/', paste(strsplit(url,'/')[[1]][c(-1,-2,-3)],collapse="/"),sep="")
		if(isTRUE(file.exists(fullPath))){
			# File already exists in synapse cache
			urls[i] <- NA
		}else if(isTRUE(file.exists(fname)) ){
			# Move the file from local cache to synapse cache
			cat(fname," exists... moving...\n")
			dir.create(dirname(fullPath),recursive=TRUE)
			file.copy(fname, dirname(fullPath))
			urls[i] <- NA
		}else if(isTRUE(file.exists(fname2))){
			cat(fname," exists... moving...\n")
			dir.create(dirname(fullPath),recursive=TRUE)
			file.copy(fname2, dirname(fullPath))
			urls[i] <- NA
		}else{
			# File does not exist yet.  Download it and move it to the appropriate place
			cat("Downloading ", fname)
			string <- gsub("https://", paste("https://", user,':', pwd, '@', sep=""), url)
			download <- try(synapseClient:::.curlWriterDownload(string,fullPath))
			urls[i] <- string
		}
	}
#	urlsToDownload <- urls[which(!is.na(urls))]
	# Writes out groups of curl calls to files.  This is how I'm multithreading... which can be improved upon.
#	splits <- rep(1:8,times=7)[1:length(urlsToDownload)]
#	sapply(1:max(splits), function(i){
#				lUrls <- urlsToDownload[which(splits==i)]
#				write.table(lUrls, quote=FALSE, col.names=FALSE, row.names=FALSE, file=paste("cmdsNext",i,".sh",sep=""))
#				lUrls
#			}) -> eme
	return(NA)
}
