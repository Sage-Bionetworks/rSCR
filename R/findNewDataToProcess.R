

results <- findNewDataToProcess()

findNewDataToProcess <- function(...){
	qry <- synapseQuery('select id, name from folder where folder.parentId=="syn1450029"')
	qry <- qry[-1,]
	write.table(paste("entity.platform", "entity.id", "entity.name", "entity.numSamples",sep="\t"), file="results.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
	for(i in 1:nrow(qry)){
		cat(i,"\t", qry$folder.name[i],"\n")
		qry2 <- synapseQuery(paste('select id, name from folder where folder.parentId=="', qry$folder.id[i], '"', sep=""))
		for(j in 1:nrow(qry2)){
			qry3 <- synapseQuery(paste('select id, name from folder where folder.name=="batchProcessed" and folder.parentId=="', qry2$folder.id[j],'"', sep=""))
			if(is.null(qry3)){
				next;
			}
			newData <- .findNewDataToProcessWithinBatchProcessedFolder(qry3$folder.id)
			if(!is.na(newData[1])){
				sapply(newData, function(x){
							ent <- getEntity(x)
							res <- paste(qry2$folder.name[j], propertyValue(ent, 'id'), 
									propertyValue(ent, 'name'), propertyValue(ent, 'numSamples'),sep="\t")
							write.table(res, file="results.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
							cat(res,"\n")
						})
			}
		}
	}
	results <- read.table("results.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
	
	# Adds new TCGA layers to the SCR
	entName <- 'newDataToProcess'
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	if(!is.null(qry)){
		newTcgaDataEntity<- getEntity(qry$entity.id)
	}else{
		newTcgaDataEntity<- Data(list(name = entName, parentId = 'syn1452692'))
	}
	newTcgaDataEntity<- addFile(newTcgaDataEntity,'results.txt')
	newTcgaDataEntity<- addObject(newTcgaDataEntity,results)
	if(is.null(qry)){
		newTcgaDataEntity<- createEntity(newTcgaDataEntity)		
	}else{
		newTcgaDataEntity<- updateEntity(newTcgaDataEntity)
	}
	newTcgaDataEntity<- storeEntity(newTcgaDataEntity)
	
	
}

.findNewDataToProcessWithinBatchProcessedFolder <- function(inputFolderId){
	# Accepts as input a batchProcessed folder. This folder needs to have the fromLevel annotation set properly.	
	inputFolder <- getEntity(inputFolderId)
	alreadyProcessedBatches <- synapseQuery(paste('select * from entity where entity.parentId=="', inputFolderId, '"',sep=""))
	fromLevel <- annotValue(inputFolder, 'fromLevel')
	fromFolderId <- synapseQuery(paste('select id, name, platform from entity where entity.parentId=="', propertyValue(inputFolder, 'parentId'),
					'" and entity.name=="', fromLevel,'"',sep=""))
	
	# Get the most recent batches.
	mostRecentBatches <- .getMostRecentBatches(fromFolderId$entity.id)
	if(is.null(alreadyProcessedBatches)){
		return(mostRecentBatches$entity.id)
	}
	# Handle the agilentg450_2 case
	if(alreadyProcessedBatches$entity.platform[1] == "agilentg4502a_07"){
		platformVersion <- rep(1, nrow(alreadyProcessedBatches))
		platformVersion[grep('agilentg4502a_07_2', tolower(alreadyProcessedBatches$entity.name))] <- 2
		platformVersion[grep('agilentg4502a_07_3', tolower(alreadyProcessedBatches$entity.name))] <- 3
		alreadyProcessedBatches$entity.serialIndex <- as.character(platformVersion + as.numeric(alreadyProcessedBatches$entity.serialIndex)/10)
	}
	
	# Are there new batches to process
	newSerialIndices <- setdiff(mostRecentBatches$entity.serialIndex, 
			alreadyProcessedBatches$entity.serialIndex)
	newBatches <- NA
	
	if(length(newSerialIndices) > 0){
		# New batches to process
		ids <- match(newSerialIndices, mostRecentBatches$entity.serialIndex)
		newBatches <- mostRecentBatches$entity.id[ids]
	}
	
# Align the two queries so data from the same serial index is in the same row
	inCommon <- intersect(mostRecentBatches$entity.serialIndex, alreadyProcessedBatches$entity.serialIndex)
	alreadyProcessedBatches <- alreadyProcessedBatches[match(inCommon, alreadyProcessedBatches$entity.serialIndex),]
	mostRecentBatches <- mostRecentBatches[match(inCommon, mostRecentBatches$entity.serialIndex),]
	
# Are there any updated serial indices?
	updatedBatches <- NA
	if(any(as.numeric(mostRecentBatches$entity.revisionNumber) > as.numeric(alreadyProcessedBatches$entity.revisionNumber))){
		# Some batch has been updated
		ids <- which(as.numeric(mostRecentBatches$entity.revisionNumber) > as.numeric(alreadyProcessedBatches$entity.revisionNumber))
		updatedBatches <- mostRecentBatches$entity.id[ids]
	}
	
	if(is.na(newBatches)[1] & is.na(updatedBatches)){
		return(NA)
	}else{
		batchesToProcess <- c(newBatches, updatedBatches)
		batchesToProcess <- batchesToProcess[!is.na(batchesToProcess )]
		return(	batchesToProcess)
	}
}



.getMostRecentBatches <- function(folderId){
	qry <- synapseQuery(paste('select id,name, revisionNumber, serialIndex, tcgaVersion, numSamples, platform from entity where entity.parentId=="', folderId,'"',sep=""))
	if(qry$entity.platform == 'agilentg4502a_07'){
		# Handle situation where agilentg4502a platform has three versions.  
		platformVersion <- rep(1, nrow(qry))
		platformVersion[grep('agilentg4502a_07_2', tolower(qry$entity.name))] <- 2
		platformVersion[grep('agilentg4502a_07_3', tolower(qry$entity.name))] <- 3
		oldSerialIndices <- qry$entity.serialIndex
		qry$entity.serialIndex <- as.character(platformVersion + as.numeric(qry$entity.serialIndex)/10)
		si2row <- split(1:nrow(qry), qry$entity.serialIndex)
#		qry$entity.serialIndex <- oldSerialIndices
	}else{
		si2row <- split(1:nrow(qry), qry$entity.serialIndex)
	}
	mostRecentEntities <- matrix(NA, nr=length(si2row), nc=ncol(qry))
	colnames(mostRecentEntities) <- colnames(qry)
	for(i in 1:length(si2row)){
		qry2 <- qry[si2row[[i]],]
		id <- max(qry2$entity.revisionNumber)
		qry2 <- qry2[which(qry2$entity.revisionNumber == id),]
		if(nrow(qry2) > 0){
			# Find one with largest version
			id <- which.max(qry2$entity.tcgaVersion)
			mostRecentEntities[i,] <- as.matrix(qry2[id,])
		}else{
			mostRecentEntities[i,] <- as.matrix(qry2)
		}
	}
	as.data.frame(mostRecentEntities, stringsAsFactors=FALSE)
}



.deleteThis <- function(x){
	qry <- synapseQuery('select id, name from folder where folder.parentId=="syn1450029"')
	sapply(qry$folder.id, function(x){ 
				qry2 <- synapseQuery(paste('select id, name from folder where folder.parentId=="', x, '"', sep=""))
				qry2$folder.name
			}) -> subFolders
	ourWorkflows <- c('pd.genomewidesnp.6', 'tissue_images', 'hgu133plus2', 'clinicalFiles', 'mageTabFiles', 'agilentg4502a_07', 'hthgu133a', 'hg-cgh-244a', 'hg-cgh-415k_g4124a', 'huex10stv2', 'h-mirna_8x15k', 'h-mirna_8x15kv2', 'cgh-1x1m_g4447a')
	needWorkflows <- setdiff(unique(unlist(subFolders)), ourWorkflows)
	
	qry <- qry[-1,]
	for(i in 2:nrow(qry)){
		cat(qry$folder.name[i],"\n")
		qry2 <- synapseQuery(paste('select id, name from folder where folder.parentId=="', qry$folder.id[i], '"', sep=""))
		folders <- setdiff(qry2$folder.name, ourWorkflows)
		qry2 <- qry2[match(folders, qry2$folder.name),]
		if(nrow(qry2)==0){ next }
		for(j in 1:nrow(qry2)){
			qry3 <- synapseQuery(paste('select id, name from folder where folder.parentId=="', qry2$folder.id[j],'"', sep=""))
			qry3 <- qry3[grep("level", qry3$folder.name),]
			rootFolder <- getEntity(qry2$folder.id[j])
			newFolder <- Folder(list(name="batchProcessed", parentId=qry2$folder.id[j]))
			annotValue(newFolder, 'platform') <- annotValue(rootFolder, 'platform')	
			annotValue(newFolder, 'tissueType') <- annotValue(rootFolder, 'tissueType')	
			annotValue(newFolder, 'disease') <- annotValue(rootFolder, 'disease')	
			annotValue(newFolder, 'acronym') <- annotValue(rootFolder, 'acronym')	
			annotValue(newFolder, 'repository') <- annotValue(rootFolder, 'repository')
			ids <- which.max(as.numeric(gsub("level_", "", qry3$folder.name)))
			annotValue(newFolder, 'fromLevel') <- qry3$folder.name[ids]
			cat(qry2$folder.name[j], "\t", qry3$folder.name[ids],"\n")
			newFolder <- createEntity(newFolder)
			cat(propertyValue(newFolder, 'id'), "\n")
		}
	}
}
