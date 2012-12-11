#1) Get all the entity ids for the new folders in the SCR.
allFolders <- synapseQuery('select id, name from folder where folder.parentId=="syn1491485"')
#2) Read in all the cancer studies we currently have
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

#3) Get the folders for the previously curated breast and prostate data
folders <- synapseQuery('select id, name from folder where folder.parentId=="syn1353107"')
lapply(folders$folder.id, function(x){ 
		qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="', x,'"',sep=""))
		}) -> res

toMove <- matrix(NA, nrow=sum(sapply(res,nrow)), nc=3)
rowCounter <- 1
for(i in 1:length(res)){
	for(j in 1:nrow(res[[i]])){
		toMove[rowCounter,] <- c(res[[i]][j,1],res[[i]][j,2],folders$folder.name[i])
		rowCounter <- rowCounter + 1
	}
}
toMove[which(!grepl("^GSE", toMove[,1])),1] <- paste(toMove[which(!grepl("^GSE", toMove[,1])),3],'_',toMove[which(!grepl("^GSE", toMove[,1])),1],sep="")
toMove[,3] <- gsub("_tmp","",toMove[,3])
toMove <- cbind(toMove, allFolders$folder.id[match(toMove[,3], allFolders$folder.name)])
toMove <- toMove[-c(4:7),]
toMove <- toMove[-c(4:6,8:11),]
newName <- paste(propertyValue(ent, 'name'),'tsv',sep=".")

newEnts <- list()
# Move each one over
for(i in 9:nrow(toMove)){
	cat(i, "\t", toMove[i,], "\n")
	inputId <- toMove[i,2]; 
	parentId <- toMove[i,4]; 
	newName <- toMove[i,1]
	newEnts[[i]] <- .buildNewEntity(inputId, parentId, newName)
	newEnts[[i]] <- storeEntity(newEnts[[i]])
}

# Now get all the processed, but not moved over entities.
procQry <- synapseQuery('select id, name from entity where entity.parentId=="syn1353107"')


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


.moveMetadata <- function(inputId, study,parentId){
	parent <- getEntity(parentId)
	ent <- loadEntity(inputId)
	content <- read.table(file.path(ent$cacheDir, ent$files),sep="\t",stringsAsFactors=FALSE)
	newName <- paste(study,'_metadata.tsv', sep="")
	eme <- .writeEntity(content, newName)
	ent <- Data(list(name=newName,parentId=parentId))
	annotValue(ent, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	ent <- addFile(ent, newName)
	annotValue(ent, 'study') <- study
	annotValue(ent, 'repository') <- annotValue(parent,'repository')
	annotValue(ent, 'status') <- 'metadata'
	ent <- storeEntity(ent)
}

.writeEntity <- function(content, newName){
	write.table(content, sep="\t",quote=FALSE, file=newName)
}

.buildNewEntity <- function(inputId, parentId, newName){
	ent <- loadEntity(inputId)
	dat <- try(exprs(ent$objects$eset),silent=TRUE)
	if(class(dat) == "try-error"){
		dat <- try(exprs(ent$objects[[1]]$eset),silent=TRUE)
	}
	if(class(dat) == "try-error"){
		cat("ERROR!!!")
		return(ent)
	}
	rownames(dat) <- gsub("_mt", "_eg", rownames(dat))
	dat <- round(dat,3)
	write.table(dat,
			file=newName,
			sep="\t",
			quote=FALSE)
	newEnt <- Data(list(name=newName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- propertyValue(ent, 'platform')
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- propertyValue(ent, 'numSamples')
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,newName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	newEnt
}
