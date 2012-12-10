contributeTcga <- function(user='anonymous', pwd='anonymous') {
	parentId2='syn1450029'
	entName <- 'tcgaPublicCrawlerOutput'
	qry1 <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	entName <- 'tcgaPrivateCrawlerOutput'
	qry2 <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	tcgaPublic <- loadEntity(qry1$entity.id)$objects$tcga
	tcgaPrivate <- loadEntity(qry2$entity.id)$objects$tcga
	
	# Remove any entity from private if its also in public
	ids <- match(tcgaPublic$data.name, tcgaPrivate$data.name)
	tcgaPrivate <- tcgaPrivate[-ids[!is.na(ids)],]
	tcgaPrivate$private <- TRUE
	tcgaPublic$private <- FALSE

	# Merge into one data frame
	tcga <- merge(tcgaPrivate, tcgaPublic, all=TRUE)	

	# Remove data which has been superseeded
	ids <- which(tcga$data.platform == "bcgsc.ca")
	if(any(!is.na(ids))){
		tcga <- tcga[-ids[!is.na(ids)],]
	}
	
	# Clean up platform and study names
	tcga$data.platform[grep("agilentg450", tcga$data.platform)] <- 'agilentg4502a_07'
	tcga$data.platform[grep("illuminadnamethylation", tcga$data.platform)] <- 'illuminadnamethylation'
	tcga$study.name[grep("TCGA_SARC", tcga$study.name)] <- 'TCGA_Sarcoma'
	tcga$tissueType[grep("TCGA_Sarcoma", tcga$study.name)] <- 'Unknown'
	tcga$disease[grep("KICH", tcga$acronym)] <- 'Cancer'
	tcga$tissueType[grep("TCGA_Controls", tcga$study.name)] <- 'Unknown'
	
	# Make studies
	tmp <- .makeStudies(parentId2, tcga)
	
	# Remove any entities that have already been created	
	createdByPrincipalId <- 273975
	qry <- .bigQuery('select id, name from entity where entity.benefactorId == "syn1450028" and entity.repository == "TCGA"', limit=1000)
	ids <- match(qry$entity.name, tcga$data.name)
	if(any(!is.na(ids))){
		tcga <- tcga[-ids[!is.na(ids)],]
	}
	if(nrow(tcga) == 0){
		return("Everything up to date!")
	}
	
	# Build the entities
	newEntities <- vector(mode = "character", length = nrow(tcga))
	for(i in 1:nrow(tcga)){
		cat("\r", i)
		res <- try(.contributeData(tcga, i, user, pwd, newEntities), silent=TRUE)
		numRetries <- 1;
		while(class(res) == "try-error") {
			res <- try(.contributeData(tcga, i, user, pwd, newEntities), silent=TRUE)
			numRetries <- numRetries + 1
			if(numRetries > 10){ 
				stop("Tried 10 times to build the entity with no success!")
			}
		}
	}
	return(ent)
}


.bigQuery <- function(qry, limit=100) 
{
	allQryRes <- matrix()
	offset <- 1
	for(i in 1:100000){
		cat("\rQuery",i)
		qryRes <- synapseQuery(paste(qry, 'limit', limit, 'offset', offset))
		offset <- offset + limit 
		if(i == 1){
			allQryRes <- qryRes
		}else{
			allQryRes <- merge(allQryRes, qryRes, all=TRUE)
		}
		if(nrow(qryRes) != limit){
			break;
		}
	}
	allQryRes
}


.contributeData <- function(tcga,i, user, pwd, newEntities){
	# Get folder id and entity
	qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="', 
					parentId2, '" and folder.name=="', tcga$study.name[i], '"', sep=""))
	folderId <- qry$folder.id
	folderEnt <- getEntity(folderId)
	# Get the subfolder entity
	subFolderEnt <- .makeSubFolders(tcga, folderId, folderEnt,i)
	subFolderId <- propertyValue(subFolderEnt, 'id')
	# Add the entity if its new, otherwise do nothing
	qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="', 
					subFolderId, '" and entity.name=="',
					tcga$data.name[i], '"', sep=""))
	if(is.null(qry)){
		# Create the data entities
		dataEntity <- Data(list(name=tcga$data.name[i], 
						parentId=subFolderId, 
						platform=tcga$data.platform[i],
						disease="cancer",
						tissueType=tcga$tissueType[i],
						species="Homo sapiens"))
		if(tcga$data.type[i] != "C"){
			propertyValue(dataEntity, 'numSamples') <- tcga$data.numSamples[i]
		}
		annotValue(dataEntity, 'repository') <- 'TCGA'
		annotValue(dataEntity, 'tcgaLevel') <- tcga$data.tcgaLevel[i]
		annotValue(dataEntity, 'lastUpdate') <- tcga$data.lastUpdate[i]
		annotValue(dataEntity, 'private') <- tcga$private[i]
		dataEntity <- .addBatchInfo(dataEntity)
		if(!is.na(tcga$data.md5[i])){
			# Add url as an external location or pointer to the files on the tcga ftp site
			dataEntity <- .addExternalLocation(dataEntity, tcga$data.url[i], 
					tcga$data.md5[i])
		}else if(tcga$data.add[i] == "TRUE"){
			# Here you have to download the file, then add the file to the entity.
			dataEntity <- downloadExternalFile(dataEntity, tcga$data.url[i], 
					tcga$data.md5[i], maxFileSize=(20 * 1073741824), 
					numRetries=20, store=FALSE, private=tcga$private[i], user=user, pwd=pwd)
		}else{
			newEntities[i] <- "ERROR"
		}
		newEntities[i] <- propertyValue(dataEntity,'id')
	}else{
		# It already exists then do nothing
		newEntities[i] <- NA
	}
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

.addExternalLocation <- function(entity, url, md5){
	fileSize <- .getFileSize(url)
	annotValue(entity,'fileSize') <- .prettifyFileSize(fileSize)
	entity <- createEntity(entity)
	propertyValue(entity, 'md5') <- md5
	propertyValue(entity, 'locations') <- list(list(path=url, type="external"))
	entity <- updateEntity(entity)
	return(entity)	
}

.makeSubFolders <- function(tcga, folderId, folderEnt, i){
	# Get the subfolders, build if necessary
	# If its a clinical entity, then build clinicalFiles/, return clinical folder id.
	level <- tcga$data.tcgaLevel[i]
	if((level=="clinical" | is.na(tcga$data.platform[i])) & !grepl("mage-tab", tcga$data.name[i]) )  {
		qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="', folderId, '" and folder.name=="clinicalFiles"', sep=""))
		if(is.null(qry)){
			# Build the clinical folder
			clinFolder <- Folder(list(name="clinicalFiles", 
							parentId=folderId))
			annotValue(clinFolder, 'acronym') <- annotValue(folderEnt, 'acronym')
			annotValue(clinFolder, 'tissueType') <- annotValue(folderEnt, 'tissueType')
			annotValue(clinFolder, 'disease') <- annotValue(folderEnt, 'disease')
			annotValue(clinFolder, 'repository') <- annotValue(folderEnt, 'repository')
			clinFolder <- createEntity(clinFolder)
		}else{
			clinFolder <- getEntity(qry$folder.id)
		}
		return(clinFolder)
	}else if(grepl("mage-tab", tcga$data.name[i])){
		# Its a mage-tab file
		qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="', folderId, '" and folder.name=="mageTabFiles"', sep=""))
		if(is.null(qry)){
			# Build the mageTabFiles folder
			mageTabFolder <- Folder(list(name="mageTabFiles", 
							parentId=folderId))
			annotValue(mageTabFolder, 'acronym') <- annotValue(folderEnt, 'acronym')
			annotValue(mageTabFolder, 'tissueType') <- annotValue(folderEnt, 'tissueType')
			annotValue(mageTabFolder, 'disease') <- annotValue(folderEnt, 'disease')
			annotValue(mageTabFolder, 'repository') <- annotValue(folderEnt, 'repository')
			mageTabFolder <- createEntity(mageTabFolder)
		}else{
			mageTabFolder <- getEntity(qry$folder.id)
		}
		return(mageTabFolder)
	}else{
		qry <- synapseQuery(paste('select id, name, platform from folder where folder.parentId=="', folderId, '" and folder.name=="',tcga$data.platform[i],'"', sep=""))
		if(is.null(qry)){
			# The platform folder does not exist.  So create it.
			platformFolder <- Folder(list(name=tcga$data.platform[i], 
							parentId=folderId))
			annotValue(platformFolder, 'acronym') <- annotValue(folderEnt, 'acronym')
			annotValue(platformFolder, 'tissueType') <- annotValue(folderEnt, 'tissueType')
			annotValue(platformFolder, 'disease') <- annotValue(folderEnt, 'disease')
			annotValue(platformFolder, 'platform') <- tcga$data.platform[i]
			annotValue(platformFolder, 'repository') <- annotValue(folderEnt, 'repository')
			platformFolder <- createEntity(platformFolder)
		}else{
			platformFolder <- getEntity(qry$folder.id)
		}
		qry2 <- synapseQuery(paste('select id, name, platform from folder where folder.parentId=="', propertyValue(platformFolder,'id'), '" and folder.name=="',level,'"', sep=""))
		if(is.null(qry2)){
			# level folder does not exist.  So create it.
			levelFolder <- Folder(list(name=level, 
							parentId=propertyValue(platformFolder, 'id')))
			annotValue(levelFolder, 'acronym') <- annotValue(folderEnt, 'acronym')
			annotValue(levelFolder, 'tissueType') <- annotValue(folderEnt, 'tissueType')
			annotValue(levelFolder, 'disease') <- annotValue(folderEnt, 'disease')
			annotValue(levelFolder, 'repository') <- annotValue(folderEnt, 'repository')
			annotValue(levelFolder, 'platform') <- tcga$data.platform[i]
			levelFolder <- createEntity(levelFolder)				
		}else{
			levelFolder <- getEntity(qry2$folder.id)
		}
		return(levelFolder)
	}
}

.makeStudies <- function(parentId2, tcga){
	uniqueStudies <- unique(tcga$study.name)
	qry <- synapseQuery(paste('select id, name from folder where folder.parentId=="', parentId2,'"',sep=""))
	uniqueStudies <- setdiff(uniqueStudies, qry$folder.name)
	tmp <- tcga[match(uniqueStudies, tcga$study.name),]
	if(nrow(tmp)==0){
		return(NULL)
	}
	for(i in 1:nrow(tmp)){
		cat("\r", i)
		folderEnt <- Folder(list(name=tmp$study.name[i],
						parentId=parentId2))
		annotValue(folderEnt, 'disease') <- tmp$disease[i] 
		annotValue(folderEnt, 'tissueType') <- tmp$tissueType[i] 
		annotValue(folderEnt, 'repository') <- tmp$repository[i]
		annotValue(folderEnt, 'acronym') <- tmp$acronym[i]
		annotValue(folderEnt, 'lastUpdate') <-  tmp$data.lastUpdate[i]
		annotValue(folderEnt, 'numPatients') <-  tmp$study.numPatients[i]
		folderEnt <- createEntity(folderEnt)		
	}
}


.getBatchInfo <- function(name) {
	if(!grepl('Level',name) & !grepl("mage-tab", name)){
		return(name)
	}
	domain <- strsplit(name,"_")[[1]][1]
	if(grepl('Level',name)){
		m <- regexpr("Level_\\d\\.\\d+\\.\\d+\\.\\d+", name,perl=TRUE)
	}else{
		m <- regexpr("mage-tab\\.\\d+\\.\\d+\\.\\d+", name,perl=TRUE)
	}
	if(m[1] == -1){
		return(c(NA,NA))
	}
	mtch <- regmatches(name, m)
	mtch <-strsplit(mtch, '\\.')[[1]]
	serialIndex <- mtch[2]
	revisionNumber <- mtch[3]
	versionNumber <- mtch[4]
	c(serialIndex, revisionNumber, domain, versionNumber)
}

.addBatchInfo <- function(dataEntity){
	info <- .getBatchInfo(propertyValue(dataEntity,'name'))
	if(length(info) == 1){
		return(dataEntity)
	}
	annotValue(dataEntity, 'serialIndex') <- info[1]
	annotValue(dataEntity, 'revisionNumber') <- info[2]
	annotValue(dataEntity, 'institution') <- info[3]
	annotValue(dataEntity, 'tcgaVersion') <- info[4]
	return(dataEntity)
}


