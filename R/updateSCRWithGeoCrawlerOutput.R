
updateGEO <- function(){ 
# run and return output from the GEO crawler. 
	geoCrawlerOutput <- loadEntity('syn1468623');# crawlNcbiGeo()
	geo <- geoCrawlerOutput$objects$geo	
# Run a query to get the list of studies already created
	alreadyCreatedStudies <- synapseQuery('select id, name from study where study.parentId=="syn1124722"')
# Filter out the ones we've already created.
	alreadyCreatedFolders <- synapseQuery('select id, name from entity where entity.parentId=="syn1491485"')
	geo <- geo[setdiff(rownames(geo), alreadyCreatedFolders$entity.name),]
	
# Anything that's new should be added
	for(i in 1:nrow(geo)){
		cat("\n\n", i, "\t", rownames(geo)[i])
		res <- try(.contributeGeoStudy(geo, i, maxFileSize=(50 * 1073741824), alreadyCreatedStudies=alreadyCreatedStudies), silent=TRUE)
		numRetries <- 1;
		while(class(res) == "try-error") {
			res <- try(.contributeGeoStudy(geo, i, maxFileSize=(50 * 1073741824), alreadyCreatedStudies=alreadyCreatedStudies), silent=TRUE)
			numRetries <- numRetries + 1
			if(numRetries > 10){ 
				#stop("Tried 10 times to build the entity with no success!")
				break;
			}
		}
	}
}

.createFolderForStudy <- function(geo, i){
	# If folder exists, return the id
	qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="syn1491485" and folder.name=="',rownames(geo)[i] ,'"',sep=""))
	if(!is.null(qry)){
		return(getEntity(qry$folder.id))
	}
	# Create folder with content from crawler
	folder <- Folder(list(parentId='syn1491485',
					name=rownames(geo)[i],
					description = geo$description[i]))
	annotValue(folder, 'numSamples') = as.numeric(geo$numSamples[i])
	annotValue(folder, 'platform') = strsplit(geo$platform[i], ';')[[1]]
	annotValue(folder, 'species') = strsplit(geo$species[i], ';')[[1]]
	annotValue(folder, 'repository') <- 'NCBI GEO'
	
# Inherit annotations / properties from study
# If folder is a study in SCR
	id <- match(propertyValue(folder, 'name'), alreadyCreatedStudies$study.name)
	if(!is.na(id)){
		studyEntity <- getEntity(alreadyCreatedStudies$study.id[id])
		# Move metadata when available
	}
# Create the folder
	folder <- createEntity(folder)
	return(folder)	
}

.createRawFolder <- function(geo, i, folder){
	qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="',propertyValue(folder, 'id'),'" and folder.name=="raw data"',sep=""))
	if(!is.null(qry)){
		return(getEntity(qry$folder.id))
	}
	
	rawFolder <- Folder(list(parentId=propertyValue(folder, 'id'),
					name='raw data',
					description = geo$description[i]))
	annotValue(rawFolder, 'numSamples') = as.numeric(geo$numSamples[i])
	annotValue(rawFolder, 'platform') = strsplit(geo$platform[i], ';')[[1]]
	annotValue(rawFolder, 'species') = strsplit(geo$species[i], ';')[[1]]
	annotValue(rawFolder, 'repository') <- 'NCBI GEO'
	rawFolder <- createEntity(rawFolder)
	return(rawFolder)
}

.createRawDataEntity <- function(geo,i,folder,rawFolder){
	rawDataEntity <- Data(list(name= paste(rownames(geo)[i],'_RAW.tar',sep=""),
					parentId=propertyValue(rawFolder, 'id'),
					'numSamples' = as.numeric(geo$numSamples[i]),
					'platform' = strsplit(geo$platform[i], ';')[[1]],
					'species' = strsplit(geo$species[i], ';')[[1]]))
	annotValue(rawDataEntity, 'repository') <- 'NCBI GEO'
	annotValue(rawDataEntity, 'status') <- 'raw'
	annotValue(rawDataEntity, 'study') <- rownames(geo)[i]
	fileSizeInBytes <- .getFileSize(geo$layer.url[i])
	annotValue(rawDataEntity, 'fileSize') <- .prettifyFileSize(fileSizeInBytes)
	return(rawDataEntity)
}

.addExternalLocation <- function(geo,i, rawDataEntity, maxFileSize=(1 * 1073741824), alreadyCreatedStudies){
	id <- match(rownames(geo)[i], alreadyCreatedStudies$study.name)
	fileSizeInBytes <- .getFileSize(geo$layer.url[i])
	if(!is.na(id)){
		# Study exists, so raw data might be there already.
		qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="', alreadyCreatedStudies$study.id[id], '"',sep=""))
		# Create raw data link, check to see if it already exists, if so, simply move it 
		if(any(grepl("raw", tolower(qry$entity.name)))){
			id <- which(grepl("raw", tolower(qry$entity.name)))
			existingRawDataEntity <- getEntity(qry$entity.id[id])
			propertyValue(rawDataEntity, 'md5') <- propertyValue(existingRawDataEntity, 'md5')
			propertyValue(rawDataEntity, 'locations') <- propertyValue(existingRawDataEntity, 'locations')
			return(rawDataEntity)
		}else{
			# New data, so download to calculate md5.
			if(fileSizeInBytes < maxFileSize){
				rawDataEntity <- .addExternalLocationToGeoRawDataEntity(rawDataEntity, geo$layer.url[i], maxFileSize)
			}else{
				return(NA)
			}
		}
	}else{
		# New data, so download to calculate md5.
		if(fileSizeInBytes < maxFileSize){
			rawDataEntity <- .addExternalLocationToGeoRawDataEntity(rawDataEntity, geo$layer.url[i], maxFileSize)
		}else{
			return(NA)
		}
	}
}

.contributeGeoStudy <- function(geo, i, maxFileSize=(1 * 1073741824), alreadyCreatedStudies=NULL){
	# For each geo study
	folder <- .createFolderForStudy(geo, i)
	
	# Create the raw folder
	rawFolder <- .createRawFolder(geo, i, folder)
	
	# Create the raw data entity, but don't store in Synpase just yet.
	rawDataEntity <- .createRawDataEntity(geo,i,folder,rawFolder)
	
	# Add the external locations, checking first to see if we've previously stored it in the old SCR.
	tmp <- .addExternalLocation(geo,i,rawDataEntity, maxFileSize=maxFileSize, alreadyCreatedStudies)
	rawDataEntity <- tmp$entity
	destfile <- tmp$file
	if(class(rawDataEntity) != "Data"){
		return(NA)
	}
	# Create it in synapse
	rawDataEntity <- createEntity(rawDataEntity)
	# Move the downloaded entity into the cache so we don't have to download it again.
	path <- paste(gsub("/archive.zip_unpacked", "", rawDataEntity$cacheDir), "/", basename(geo$layer.url[i]), sep="")
	file.copy(destfile, path, overwrite=TRUE)
	cat(propertyValue(rawDataEntity,'id'), "\n")
	return(rawDataEntity)
}

.addExternalLocationToGeoRawDataEntity <- function(entity, url, maxFileSize=(1 * 1073741824), numRetries=20, keep=TRUE){
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
#	if(isTRUE(keep)){
#		file.remove(destfile)
#	}
	list(entity=entity, file=destfile)
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



.inherit <- function(study, folder) {
	if('disease' %in% names(properties(study))){
#		annotValue(folder, tolower(propertyValue(study, 'disease'))) <- TRUE
		annotValue(folder, 'disease') <- tolower(propertyValue(study, 'disease'))
	}
	
	if('tissueType' %in% names(properties(study))){
#		annotValue(folder, tolower(propertyValue(study, 'tissueType'))) <- TRUE
		annotValue(folder, 'tissueType') <- tolower(propertyValue(study, 'tissueType'))	
	}
	return(folder)
}


.inherit2 <- function(data1, data2) {
	if('disease' %in% names(properties(data1))){
#		annotValue(folder, tolower(propertyValue(study, 'disease'))) <- TRUE
		propertyValue(data2, 'disease') <- tolower(propertyValue(data1, 'disease'))
	}
	
	if('tissueType' %in% names(properties(data1))){
#		annotValue(folder, tolower(propertyValue(study, 'tissueType'))) <- TRUE
		propertyValue(data2, 'tissueType') <- tolower(propertyValue(data1, 'tissueType'))	
	}
	return(data2)
}

