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


