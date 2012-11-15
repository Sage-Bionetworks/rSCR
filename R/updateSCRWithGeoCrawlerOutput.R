source("~/login.R"); 
source("~/dataActLogin.R")

# run and return output from the GEO crawler. 
geoCrawlerOutput <- loadEntity('syn1468623');# crawlNcbiGeo()
geo <- geoCrawlerOutput$objects$geo

# Run a query to get the list of studies already created
alreadyCreatedStudies <- synapseQuery('select id, name from study where study.parentId=="syn1124722"')

# Filter out the ones we've already created.
alreadyCreatedFolders <- synapseQuery('select id, name from entity where entity.parentId=="syn1453669"')
geo <- geo[setdiff(rownames(geo), alreadyCreatedFolders$entity.name),]

# Anything that's new should be added
for(i in 2:nrow(geo)){
	cat("\n\n", i)
	res <- try(.contributeGeoStudy(geo, i), silent=TRUE)
	numRetries <- 1;
	while(class(res) == "try-error") {
		res <- try(.contributeGeoStudy(geo, i), silent=TRUE)
		numRetries <- numRetries + 1
		if(numRetries > 4){ 
			#stop("Tried 10 times to build the entity with no success!")
			break;
		}
	}
}

.contributeGeoStudy <- function(geo, i, maxFileSize=(1 * 1073741824)){
	# For each geo study
	# Create folder with content from crawler
	folder <- Folder(list(parentId='syn1453669',
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
		folder <- .inherit(studyEntity, folder)
		# Move metadata when available
	}
	# Create the folder
	folder <- createEntity(folder)
	
	# Set the fileSize annotation
	rawDataEntity <- Data(list(name='rawData', 
					parentId=propertyValue(folder, 'id'),
					'numSamples' = as.numeric(geo$numSamples[i]),
					'platform' = strsplit(geo$platform[i], ';')[[1]],
					'species' = strsplit(geo$species[i], ';')[[1]]))
	annotValue(rawDataEntity, 'repository') <- 'NCBI GEO'
	fileSizeInBytes <- .getFileSize(geo$layer.url[i])
	annotValue(rawDataEntity, 'fileSize') <- .prettifyFileSize(fileSizeInBytes)
	if(!is.na(id)){
		# Study exists, so raw data might be there already.
		qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="', propertyValue(studyEntity, 'id'),'"',sep=""))
		# Create raw data link, check to see if it already exists, if so, simply move it 
		if(any(grepl("raw", tolower(qry$entity.name)))){
			id <- which(grepl("raw", tolower(qry$entity.name)))
			existingRawDataEntity <- getEntity(qry$entity.id[id])
			propertyValue(rawDataEntity, 'md5') <- propertyValue(existingRawDataEntity, 'md5')
			propertyValue(rawDataEntity, 'locations') <- propertyValue(existingRawDataEntity, 'locations')
			rawDataEntity <- .inherit2(existingRawDataEntity, rawDataEntity)
		}else{
			# New data, so download to calculate md5.
			if(fileSizeInBytes < maxFileSize){
				rawDataEntity <- .addExternalLocationToGeoRawDataEntity(rawDataEntity, geo$layer.url[i])
			}else{
				return(folder)
			}
		}
	}else{
		# New data, so download to calculate md5.
		if(fileSizeInBytes < maxFileSize){
			rawDataEntity <- .addExternalLocationToGeoRawDataEntity(rawDataEntity, geo$layer.url[i])
		}else{
			return(folder)
		}
	}
	rawDataEntity <- createEntity(rawDataEntity)
	cat(propertyValue(rawDataEntity,'id'), "\n")
	return(rawDataEntity)
}

.addExternalLocationToGeoRawDataEntity <- function(entity, url, maxFileSize=(1 * 1073741824), numRetries=20){
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

