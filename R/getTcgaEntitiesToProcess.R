getTcgaEntitiesToProcess <- function(private=FALSE) {
	parentId1='syn301192'
	parentId2='syn150935'
	if(isTRUE(private)){ 
		parentId1='syn310038'
		parentId2='syn308810'
	}
	
	# Gets output from tcga crawler.  All of these have already been contributed to the SCR.  This is faster than crawling SCR to find all entities.
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="tcgaCrawlerOutput" and entity.parentId=="', parentId1, '"', sep=""))
	ent <- loadEntity(qry$entity.id)
	tcga <- ent$objects$tcga	
	tcga <- tcga[!duplicated(tcga$data.name),]
	rownames(tcga) <- tcga$data.name
	
	# Gets the list of studies already in the SCR
	tcgaQry <- synapseQuery(paste('select id, name from study where study.repository == "TCGA" and study.parentId == "',parentId2,'"',sep=""))
	newEntities <- oldEntities <- NA
	
	# Now we want to determine if the data identified by our crawler exists in the SCR.
	for(i in 1:nrow(tcgaQry)){
		# Announce which cancer is being processed
		x <- tcgaQry[i,]
		cat(paste(x,collapse="\t"),"\n")
		
		# Get the set of existing entities already in the SCR, then filter to remove any that are not level 1 data.
		qry <- paste('select id, name, platform, numSamples from entity where entity.parentId == "', x[2],'"',sep="")
		existingEntities <- synapseQuery(qry)
		existingEntities <- existingEntities[grep("Level_1",existingEntities$entity.name),]
		existingEntities <- existingEntities[!grepl('tissue_images',existingEntities$entity.name),]
		rownames(existingEntities) <- existingEntities$entity.name
		
		# Determine which of the ones identified by our crawler already exist in the SCR.
		# Used to remove SNP6.0 Level 1 files, which we actually process from the private-SCR project.
		# We might want to delete these entities at some point in the future.
		th <- intersect(existingEntities$entity.name, rownames(tcga))
		existingEntities <- existingEntities[th,]
		existingEntities <- existingEntities[which(existingEntities$entity.platform != "NULL"),]
		
		# If there aren't any new entities then move onto the next one.
		if(nrow(existingEntities) == 0){
			next
		}
		
		# For each platform, get the lastest versions
		plat2row <- split(1:nrow(existingEntities), as.character(unlist(existingEntities$entity.platform)))
		sapply(plat2row, function(p){
					existingEntities2 <- existingEntities[p,]
					batchInfo <- sapply(existingEntities2$entity.name, .getBatchInfo)
					serial2col <- split(1:ncol(batchInfo), batchInfo[1,])
					serial2rev <- split(as.numeric(batchInfo[2,]), batchInfo[1,])
					rev2col <- split(1:ncol(batchInfo), batchInfo[2,])		
					ids <- c(9)
					for(j in 1:length(serial2col)){
						ids[j] <- serial2col[[j]][which.max(serial2rev[[j]])]
					}
					p[ids]
				}) -> ids
		
		# Remove all but the latest version of each entity.
		existingEntities <- existingEntities[unlist(ids),]
		decisions <- matrix(NA, nr=nrow(existingEntities), nc=ncol(existingEntities)+1)
		colnames(decisions) <- c("entity.numSamples", "entity.platform", "entity.name","entity.id", "entity.status")
		# For each potential entity, determine if it needs to be processed.
#		apply(existingEntities, 1, function(x){ 
		for(j in 1:nrow(existingEntities)){
			# Build an entity and set name equal to its corresponding name in metaGenomics.
			newEnt <- .makeMetaGenomicsEntityName(existingEntities[j,"entity.id"])
			# Retrieve the revision number for this entity.
			newRevNumber <- as.numeric(annotValue(newEnt, 'revisionNumber'))
			# Get the parent entity, which should be the TCGA study in SCR
			qry <- synapseQuery(paste('select name, acronym from study where study.id == "',propertyValue(newEnt,"parentId"),'"',sep=""))
			parentName <- qry$study.name
			parentId <- qry$study.id
			# Get the corresponding study from metaGenomics and query for the entity id as it might already exist.
			mGStudy <- synapseQuery(paste('select id, name from study where study.parentId == "syn275039" and study.name=="',parentName,' - SNM Normalized"',sep=""))$study.id
			existingEntityId <- synapseQuery(paste('select id, name, revisionNumber from entity where entity.parentId=="', mGStudy, '" and entity.name=="', propertyValue(newEnt,'name'), '"', sep=""))$entity.id					
			if(is.null(existingEntityId)){
				# New entity, so it should be processed
				#cat(x, "\t", newRevNumber, "\t", 'new', "\n")
				decisions[j,1:4] <- as.character(unlist(existingEntities[j,]))
				decisions[j,5] <- 'new'
			}else{
				# Old entity, so it should only be processed if its an updated revision of the corresponding batch.
				ent <- getEntity(existingEntityId)
				existingRevNumber <- as.numeric(annotValue(ent, 'revisionNumber'))
				if(newRevNumber > existingRevNumber) {
					# Updated version of an existing entity, so it should be processed
					decisions[j,1:4] <- as.character(unlist(existingEntities[j,]))
					decisions[j,5] <- 'new'
				}else{
					# Old entity
					#cat(x, "\t", newRevNumber, "\t", 'old', "\n")
					decisions[j,1:4] <- as.character(unlist(existingEntities[j,]))
					decisions[j,5] <- 'old'
				}
			}
		}
		
		# Annotate the decisions matrix and then split into new and old matrices.
		newOnes <- as.matrix(decisions[which(decisions[,5] == 'new'),])
		oldOnes <- as.matrix(decisions[which(decisions[,5] == 'old'),])
		if(nrow(newOnes) > 0){ 	
			if(ncol(newOnes)==1){
				newOnes <- t(newOnes)
			}
			if(is.na(newEntities[1])){
				newEntities <- newOnes
			}else{
				newEntities <- rbind(newEntities, newOnes)
			}
		}
		if(nrow(oldOnes) > 0){ 	
			if(ncol(oldOnes) == 1){
				oldOnes <- t(oldOnes)
			}
			if(is.na(oldEntities[1])){
				oldEntities <- oldOnes
			}else{
				oldEntities <- rbind(oldEntities, oldOnes)
			}
		}
	}
	
	# Add file to Synapse containing new entities to process.  
	# Needs to add files for both public and private tier.
	if(nrow(newEntities) > 0){
		newEntities <- newEntities[order(as.numeric(newEntities[,'entity.numSamples'])),]
		fname <- 'tcgaNewDataEntitiesToProcess'
		if(isTRUE(private)){
			fname <- 'tcgaNewDataEntitiesToProcessPrivate'
		}
		write.table(newEntities ,file="tcgaNewDataEntitiesToProcess.txt",sep="\t",row.names=FALSE, quote=FALSE, append=FALSE)                  
		qry <- synapseQuery(paste('select id, name from entity where entity.name=="',fname,'" and entity.parentId=="',parentId1,'"',sep=""))
		if(is.null(qry)){
			myLayer <- Data(list(name = "tcgaNewDataEntitiesToProcess", parentId = parentId1, type="M"))
			if(isTRUE(private)){
				myLayer <- Data(list(name = "tcgaNewDataEntitiesToProcessPrivate", parentId = parentId1, type="M"))
			}
			myLayer <- addFile(myLayer, "tcgaNewDataEntitiesToProcess.txt")
			myLayer <- createEntity(myLayer)
			newEntities <- as.data.frame(newEntities)
			myLayer <- addObject(myLayer, newEntities)
		}else{
			myLayer <- loadEntity(qry$entity.id)
			myLayer <- deleteFile(myLayer, "tcgaNewDataEntitiesToProcess.txt")
			myLayer <- addFile(myLayer, "tcgaNewDataEntitiesToProcess.txt")
			newEntities <- as.data.frame(newEntities)
			myLayer <- addObject(myLayer, newEntities)
			myLayer <- updateEntity(myLayer)
		}
		myLayer <- storeEntity(myLayer)	
	}else{
		myLayer <- NULL
	}
	myLayer
}

.getBatchInfo <- function(name) {
	if(!grepl('Level',name)){
		return(name)
	}
	domain <- strsplit(name,"_")[[1]][1]
	m <- regexpr("Level_\\d\\.\\d+\\.\\d+", name,perl=TRUE)
	if(m[1] == -1){
		return(c(NA,NA))
	}
	mtch <- regmatches(name, m)
	mtch <-strsplit(mtch, '\\.')[[1]]
	serialIndex <- mtch[2]
	revisionNumber <- mtch[3]
	c(serialIndex, revisionNumber, domain)
}


.makeMetaGenomicsEntityName <- function(id,update=FALSE) {
	ent <- getEntity(id)
	name <- propertyValue(ent,'name')
	if(!grepl('Level',name)){
		return(ent)
	}
	acronym <- annotValue(ent,'acronym')
	domain <- strsplit(name,"_")[[1]][1]
	platform <- propertyValue(ent,'platform')
	m <- regexpr("Level_\\d\\.\\d+\\.\\d+", name,perl=TRUE)
	if(m[1] == -1){
		return(c(NA,NA))
	}
	mtch <- regmatches(name, m)
	mtch <-strsplit(mtch, '\\.')[[1]]
	serialIndex <- mtch[2]
	revisionNumber <- mtch[3]
	new.name <- paste(acronym, platform, serialIndex, sep="_")
	if(propertyValue(ent,'platform') == "pd.genomewidesnp.6") {			
		new.name <- paste(new.name, 'probeRatios', sep="_") 
	}
	if(!grepl('Level',name)){
		new.name <- name
	}
	annotValue(ent,'serialIndex') <- serialIndex  
	annotValue(ent,'institution') <- domain
	annotValue(ent,'revisionNumber') <- revisionNumber
	propertyValue(ent,'name') <- new.name
	if(isTRUE(update)){
		ent <- updateEntity(ent)
	}	
	ent
}
