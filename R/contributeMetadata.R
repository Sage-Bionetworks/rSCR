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
for(i in 1:nrow(toMove)){
	cat(i, "\t", toMove[i,], "\n")
	inputId <- toMove[i,2]; 
	parentId <- toMove[i,4]; 
	newName <- toMove[i,1]
	newEnts[[i]] <- .buildNewEntity(inputId, parentId, newName)
	newEnts[[i]] <- storeEntity(newEnts[[i]])
}

# Now get all the processed, but not moved over entities.
procQry <- synapseQuery('select id, name from entity where entity.parentId=="syn1353107"')
procQry <- procQry[grepl("_processed",procQry$entity.name),]
procQry$entity.study <- gsub("_processed", "", procQry$entity.name)
toProcess <- procQry[match(allCancers[,1], procQry$entity.study),]
toProcess <- toProcess[which(!is.na(toProcess$entity.study)),]
toProcess$folder.id <- allFolders$folder.id[match(toProcess$entity.study, allFolders$folder.name)]
toProcess <- toProcess[!duplicated(toProcess$entity.id),]
toProcess$entity.platform <- allCancers$platform[match(toProcess$entity.study, allCancers$name)]
toProcess <- toProcess[which(grepl("h",toProcess$entity.platform)),]
toProcess <- toProcess[which(toProcess$entity.platform != 'hgu133plus2'),]

qry <- synapseQuery('select id, name, numSamples from entity where entity.benefactorId == "syn1450028" and entity.status=="processed"')



newEnts <- list()
# Move each one over
for(i in 32:nrow(toProcess)){
	cat("\n\n", i,"\n\n")
	inputId <- toProcess[i,2]; 
	parentId <- toProcess[i,4]; 
	newName <- paste(toProcess[i,3], toProcess[i,5], 'gene.tsv', sep="_")
	newEnts[[i]] <- .buildNewEntity2(inputId, parentId, newName)
	if(class(newEnts[[i]]) == "Data"){
		newEnts[[i]] <- storeEntity(newEnts[[i]])
	}
}
		

.buildNewEntity2 <- function(inputId, parentId, newName){
	ent <- loadEntity(inputId)
	dat <- try(exprs(ent$objects$eset),silent=TRUE)
	annot <- try(annotation(ent$objects$eset),silent=TRUE)
	if(class(dat) == "try-error"){
		dat <- try(exprs(ent$objects[[1]]$eset),silent=TRUE)
		annot <- try(annotation(ent$objects[[1]]$eset),silent=TRUE)
	}
	if(class(dat) == "try-error" | class(annot) == "try-error"){
		cat("ERROR!!!")
		return("Error")
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
	propertyValue(newEnt, 'platform') <- annot
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,newName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	newEnt
}


for(i in 6:length(ents)){
	newName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	newName <- paste(newName, '_pd.genomewidesnp.6_probe.tsv', sep="")
	parentId <- toProcess[i,4]
	ent <- ents[[i]]
	dat <- try(exprs(ent$objects$eset),silent=TRUE)
	annot <- try(annotation(ent$objects$eset),silent=TRUE)
	if(class(dat) == "try-error"){
		dat <- try(exprs(ent$objects[[1]]$eset),silent=TRUE)
		annot <- try(annotation(ent$objects[[1]]$eset),silent=TRUE)
	}
	if(class(dat) == "try-error" | class(annot) == "try-error"){
		cat("ERROR!!!")
		return("Error")
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
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,newName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'probe'
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	probeEnt <- storeEntity(newEnt)
	
	fits <- fit.pset(dat, all.gene2row)
	sapply(fits$probe.weights, function(x){
				nms <- names(x)
				if(is.null(nms)){
					rep("Unkown", length(x))
				}else{
					nms
				}
			}) -> probeNames
	
	# Probe weights
	probeWeights <- data.frame(symbols = rep(names(fits$probe.weights), 
					sapply(fits$probe.weights,length)),
			probeNames=unlist(probeNames),
			probeWeights = round(unlist(fits$probe.weights),3))
	rownames(probeWeights) <- NULL
	wtsName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	wtsName <- paste(wtsName, '_pd.genomewidesnp.6_probeWeights.tsv', sep="")	
	write.table(probeWeights, file=wtsName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=wtsName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,wtsName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(wtsName)$size)
	annotValue(newEnt, 'status') <- 'probeWeights'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	annotValue(newEnt, 'summaryType') <- 'gene'
	wtsEnt <- storeEntity(newEnt)
	
	# FIC
	fic <- fits$singular.values
	ficName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	ficName <- paste(ficName, '_pd.genomewidesnp.6_FIC.tsv', sep="")	
	fic <- round(fic,3)
	write.table(fic, file=ficName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=ficName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,ficName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(ficName)$size)
	annotValue(newEnt, 'status') <- 'fic'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	ficEnt <- storeEntity(newEnt)
	
	# Summarized Data
	dat <- fits$estimated.rna.concentration
	datName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	datName <- paste(datName, '_pd.genomewidesnp.6_gene.tsv', sep="")	
	dat <- round(dat,3)
	write.table(dat, file=datName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=datName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,datName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(datName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	dataEnt <- storeEntity(newEnt)
	
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

.moveMetadata <- function(inputId, study, parentId){
	parent <- getEntity(parentId)
	ent <- loadEntity(inputId)
	content <- read.table(file.path(ent$cacheDir, ent$files),sep="\t",stringsAsFactors=FALSE,header=TRUE)
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

