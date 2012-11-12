inputFolderId <- 'syn1461449'
inputFolderId <- 'syn1465260'

.findNewDataToProcess <- function(inputFolderId){
	
	inputFolder <- getEntity(inputFolderId)
	alreadyProcessedBatches <- synapseQuery(paste('select * from entity where entity.parentId=="', inputFolderId, '"',sep=""))
	fromLevel <- annotValue(inputFolder, 'fromLevel')
	fromFolderId <- synapseQuery(paste('select id, name from entity where entity.parentId=="', propertyValue(inputFolder, 'parentId'),
					'" and entity.name=="', fromLevel,'"',sep=""))
	mostRecentBatches <- .getMostRecentBatches(fromFolderId$entity.id)
	
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
		return("Nothing to process")
	}else{
		batchesToProcess <- c(newBatches, updatedBatches)
		batchesToProcess <- batchesToProcess[!is.na(batchesToProcess )]
		return(	batchesToProcess)
	}
}



.getMostRecentBatches <- function(folderId){ 
	qry <- synapseQuery(paste('select id,name, revisionNumber, serialIndex, tcgaVersion, numSamples from entity where entity.parentId=="', folderId,'"',sep=""))	
	si2row <- split(1:nrow(qry), qry$entity.serialIndex)
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


