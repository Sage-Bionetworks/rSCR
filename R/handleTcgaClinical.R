#
#  Should take as input a synpase entity ID.  maybe a tcga cancer name
#  It will first get the clinical data layer from synapse.
#  Next, it will make the big master file from these files.
#  Then, it will add on the info for each of the SDRF files.
#

setMethod(
		f = "handleTcgaClinical",
		signature = signature("character"),
		definition = function(entity){
			if(entity=="all") {
				# Get all TCGA datasets
				studies <- synapseQuery('select * from study where study.repository == "TCGA" and study.parentId == "syn1124721"')
				studiesMG <- synapseQuery('select * from study where study.repository == "TCGA" and study.parentId == "syn275039"')
				studies <- studies[match(studiesMG$study.acronym, studies$study.acronym),]
				mergedData <- list()
				for( i in 1:nrow(studies)){
					cat("\n\n\n", i, ":", studies$study.acronym[i],"\n\n\n")
					studyId <- studies$study.id[i]
					mergedData[[i]] <- .oneTcgaCancer(studyId)
					if(class(mergedData[[i]]) == "try-error"){
						next;
					}
					metadata <- mergedData[[i]]	
					file.name <- paste(studies$study.acronym[i],"_mergedClinical.txt",sep="")
					write.table(mergedData[[i]], file.name, sep="\t",quote=FALSE, row.names=FALSE)
					qry <- synapseQuery(paste('select id, name from entity where entity.name=="',
									gsub(".txt",'',file.name),'" and entity.parentId=="', studiesMG$study.id[i],'"',sep=""))
					if(is.null(qry)){
						data <- Data(list(parentId=studiesMG$study.id[i],
										name=gsub(".txt",'',file.name),
										tissueType=studiesMG$study.tissueType[i], 
										disease=studiesMG$study.disease[i],
										species=studiesMG$study.species[i]))
						annotValue(data,'repository') <- 'TCGA'
						annotValue(data,'acronym') <- studiesMG$study.acronym[i]
						data <- addFile(data,file.name)
						data <- addObject(data, metadata)
						data <- createEntity(data)
						data <- storeEntity(data)							
					}else{
						data <- loadEntity(qry$entity.id)
						data <- deleteFile(data, file.name)
						data <- deleteObject(data, 'metadata')
						data <- updateEntity(data)
						data <- addFile(data,file.name)
						data <- addObject(data, metadata)
						data <- updateEntity(data)
						data <- storeEntity(data)						
					}
				}
				return(mergedData)
			}else{
				# Assume contribution is a specific layer ID
				mergedData <- .oneTcgaCancer(entity)
			}
		}
)

setMethod(
		f = "handleTcgaClinical",
		signature = signature("numeric"),
		definition = function(entity){
			handleTcgaClinical(  as.character(entity)  )
		}
)

setMethod(
		f = "handleTcgaClinical",
		signature = signature("Data"),
		definition = function(entity){
			.handleTcgaClinicalLayer(  entity )
		}
)

.oneTcgaCancer <- function(studyId){
	clinicalLayers <- synapseQuery(paste('select * from entity where entity.tcgaLevel == "clinical" and entity.parentId == "', studyId, '"', sep=''))						
	id2load <- grep("^clinical_public", clinicalLayers$entity.name,perl=TRUE)[1]
	if(is.na(id2load)){
		var <- paste('Cannot find file beginning with clinical_public for ',studyId) 
		class(var) <- 'try-error'
		return(var)
	}
	clinicalData <- loadEntity(clinicalLayers$entity.id[id2load])
	cacheDir <- clinicalData$cacheDir
	if(grepl("tar.gz",list.files(cacheDir)[1])) {
		bn <- basename(tempfile())
		tDir <- paste(tempdir(),bn,sep='/')
		dir.create(tDir)
		system(paste("tar -xzf ", cacheDir,'/',list.files(cacheDir)[1], ' -C ', tDir,sep=""))
		clinicalData@archOwn@fileCache$cacheDir <- list.dirs(tDir)[1]
	}
	clinicalMerge <- .handleTcgaClinicalLayer(clinicalData)	
	# Get information from mage files
	clinicalMerge <- .handleMageFiles(clinicalLayers, clinicalMerge)
	return(clinicalMerge)
}

.handleMageFiles <- function(clinicalLayers, clinicalMerged) {
	# This function merges information from each of the MAGE files
	mageLayers <- clinicalLayers[grep("mage",tolower(clinicalLayers$entity.name)),]
	if(nrow(mageLayers) == 0){
		return(clinicalMerged)
	}
	mageLayers <- mageLayers[which(!grepl('tissue_images',mageLayers$entity.name)),]
	mageLayers <- mageLayers[which(!grepl('RPPA',mageLayers$entity.name)),]
	mageLayers <- mageLayers[which(!grepl('IlluminaGA_DNASeq',mageLayers$entity.name)),]
	mageLayers <- .findLargestVersions(mageLayers)  
	for(i in 1:nrow(mageLayers)){
		cat(mageLayers[i,"entity.name"],"\n")
		ent <- loadEntity(mageLayers[i,"entity.id"])
		cacheDir <- ent$cacheDir
		if(grepl("tar.gz",list.files(cacheDir)[1])){
			cacheDir <- ent$cacheDir				
			tDir <- tempdir()
			system(paste("tar -xzf ", cacheDir,'/',list.files(cacheDir)[1], ' -C ', tDir,sep=""))
			ent@archOwn@fileCache$cacheDir <- paste(tDir,'/',gsub(".tar.gz","",mageLayers[i,"entity.name"]),sep="")
		}
		sdrf.file <- paste(ent$cacheDir,'/',list.files(ent$cacheDir,pattern="sdrf")[1],sep="")
		mage.tab <- .read(sdrf.file)
		id <- .findColumnToMergeOn(mage.tab, clinicalMerged$bcr_aliquot_barcode)
		if(id[2] == 0) {
			id <- .findColumnToMergeOn(mage.tab, clinicalMerged$bcr_aliquot_uuid)
			cns <- setdiff(1:ncol(mage.tab), id[2])
			colnames(mage.tab)[cns] <- paste(gsub('.mage-tab',"",mageLayers$types[i]),colnames(mage.tab)[cns],sep="-")
			clinicalMerged <- merge(clinicalMerged, mage.tab, by.y=colnames(mage.tab)[id[1]], by.x="bcr_aliquot_uuid", all.x=TRUE)
		}else{
			cns <- setdiff(1:ncol(mage.tab), id[2])
			colnames(mage.tab)[cns] <- paste(gsub('.mage-tab',"",mageLayers$types[i]),colnames(mage.tab)[cns],sep="-")
			clinicalMerged <- merge(clinicalMerged, mage.tab, by.y=colnames(mage.tab)[id[1]], by.x="bcr_aliquot_barcode", all.x=TRUE)
		}
		cat("Total:", nrow(mage.tab), " Found: ", id[2], " Column: ", colnames(mage.tab)[id[1]],"\n")
#		if(sum(colnames(mage.tab) == "Comment..Aliquot.UUID.") > 0 & sum(colnames(clinicalMerged) == "bcr_aliquot_uuid") > 0){  
#			mage.tab[,"Comment..Aliquot.UUID."] <- tolower(mage.tab[,"Comment..Aliquot.UUID."])
#			cns <- which(colnames(mage.tab) != "Comment..Aliquot.UUID.")
#			colnames(mage.tab)[cns] <- paste(gsub('.mage-tab',"",mageLayers$types[i]),colnames(mage.tab)[cns],sep="-")
#			clinicalMerged <- merge(clinicalMerged, mage.tab, by.y="Comment..Aliquot.UUID.", by.x="bcr_aliquot_uuid", all.x=TRUE)
#		}else if(sum(colnames(mage.tab) == "Extract.Name") > 0 ){
#			cns <- which(colnames(mage.tab) != "Extract.Name")
#			colnames(mage.tab)[cns] <- paste(gsub('.mage-tab',"",mageLayers$types[i]),colnames(mage.tab)[cns],sep="-")
#			clinicalMerged <- merge(clinicalMerged, mage.tab, by.y="Extract.Name", by.x="bcr_aliquot_barcode", all.x=TRUE)
#		}else if(sum(colnames(mage.tab) == "Source.Name") > 0){
#			cns <- which(colnames(mage.tab) != "Source.Name")
#			colnames(mage.tab)[cns] <- paste(gsub('.mage-tab',"",mageLayers$types[i]),colnames(mage.tab)[cns],sep="-")
#			clinicalMerged <- merge(clinicalMerged, mage.tab, by.y="Source.Name", by.x="bcr_aliquot_barcode", all.x=TRUE)			
#		}else{
#			cat("Found nothing for", mageLayers[i,"entity.name"],"\n")
#		}		
	}
	clinicalMerged
}

.findColumnToMergeOn <- function(mage.tab,y){
	cts <- apply(mage.tab,2,function(x){ sum(x %in% y)})
	c(which.max(cts), max(cts))
}

.findLargestVersions <- function(allLayers){ 
	# Code finds layer prefix and version suffix.  Turns version suffix into integer.
	# Sorts by prefix and integer within prefix, with version in descending order.
	# Then finds prefixes that are not duplicated.  These correspond to highest version.
	t(sapply(allLayers[,"entity.name"], function(x){
						x <- gsub(".tar.gz","",x)
						info <- strsplit(x,"\\.")[[1]]
						i <- length(info)
						type <- paste(info[c(-((i-2):(i)))],collapse=".")
						version <- as.numeric(info[(i-2)]) * 10000000000 + as.numeric(info[(i-1)]) * 100000 + as.numeric(info[i])
						c(type, version)
					})) -> types.and.versions
	layers <- data.frame(version=as.numeric(types.and.versions[,2]),
			types = types.and.versions[,1])
	version.order <- order(layers$types, layers$version,decreasing=TRUE)
	allLayers <- allLayers[version.order,]
	layers <- layers[version.order,]
	allLayers$types <- layers$types
	tmp2 <- allLayers[!duplicated(layers$types),]
	tmp2
}

.handleTcgaClinicalLayer <- function(layer){ 
	cacheDir <- layer$cacheDir
# Step 1: Load in aliquot file.
	files <- dir(cacheDir, pattern = ".txt", full.names=TRUE)
	if(!any(grepl("aliquot",files))) {
		stop("No Aliquot file found\n");
	}
	aliquot.file <- files[grep("aliquot",files)]
	cat(aliquot.file, "\n")
	master <- .read(aliquot.file)
	
# Step 2: add patient and sample identifiers
	t(apply(master, 1, function(x){
						obj <- strsplit(x[["bcr_aliquot_barcode"]],"-")
						patient <- paste(obj[[1]][1:3], collapse="-")
						sample <- paste(obj[[1]][1:4], collapse="-")
						c(patient, sample)
					})) -> pat
	colnames(pat) <- c("bcr_patient_barcode", "bcr_sample_barcode")
	master.start <- master <- cbind(master, pat)
#	cat("Master started with", nrow(master),"rows and", ncol(master),"columns\n");
	for(i in 1:length(files)){
		if(! grepl("aliquot",files[i])) {
			cat("Adding file", files[i],"\n")
			master <- .handleDuplicates(files[i], master)
			#	cat("Master now has", nrow(master),"rows and", ncol(master),"columns\n");
			#cat(i,"\n")
		}else{
			#cat("Skipping",files[i],"\n")
		}
	}
	master[,'bcr_aliquot_uuid'] <- tolower(master[,'bcr_aliquot_uuid'])
	return(master)
}

.loadSDRF <- function(sdrfName, master){
	cat("Trying to merge sdrf file:", sdrfName,"\n");
	sdrf <- .read(sdrfName)
	if(! "Extract.Name" %in% colnames(sdrf)){
		stop(cat("Cannot find column named Extract.Name in file",sdrfName,"\n"));
	}
	id <- which(colnames(sdrf) == "Extract.Name")
	colnames(sdrf) <- paste(colnames(sdrf), strsplit(sdrfName,"\\.")[[1]][4],sep="_")
	master.updated <- merge(master, sdrf, by.x="bcr_aliquot_barcode",by.y=colnames(sdrf)[id], all.x=TRUE)
	return(master.updated)
}

.handleDuplicates <- function(fileName, master){
	if(!file.exists(fileName)) {
		stop("Cannot find file\n");
	}
	file <- .read(fileName)
	string <- "bcr_sample_barcode"
	id <- which(colnames(file) == "bcr_sample_barcode")
	if(length(id) == 0) {
		id <- which(colnames(file) == "bcr_patient_barcode")
		string <- 'bcr_patient_barcode'
	}
	# Are there multiple rows per bcr_patient_barcode?
	num.dup <- sum(duplicated(file[,id]))
	
	if(num.dup > 0){ 
		otherIds <-setdiff(1:ncol(file), id)
		sapply(otherIds, function(x){
					sapply(split(file[,x], file[,id]), function(y){
								paste(y,collapse=",")
							})
				}) -> hmm
		ids <- names(split(file[,otherIds[1]], file[,id]))  
		if(length(ids) > 1){
			ret <- cbind(ids, hmm)
			colnames(ret) <- colnames(file)[c(id,otherIds)]
		}else{
			ret <- matrix(c(ids,hmm), nr=1)
			colnames(ret) <- colnames(file)[c(id,otherIds)]
		}
		localInfo <- as.data.frame(ret)
	}else{
		localInfo <- file
	}
	in.common <- intersect(colnames(localInfo),colnames(master))
	in.common <- setdiff(in.common, c("bcr_sample_barcode", "bcr_patient_barcode"))
	if(length(in.common) > 0){ 
		for(i in 1:length(in.common)){
			x <- in.common[i]
			thisCol <- which(colnames(localInfo) == x)
			colnames(localInfo)[thisCol] <- paste(colnames(localInfo)[thisCol], basename(fileName))
		}
	}
	if(sum(grepl(string, colnames(localInfo))) != 0 & sum(grepl(string, colnames(master))) != 0){
		master.updated <- merge(master, localInfo, by=string, all.x=TRUE)
		master.updated
	}else{
		master
	}
}

.read <- function(file) {
	fileContents <- read.delim(file, 
			header=TRUE, 
			stringsAsFactors=FALSE,
			quote="",strip.white=TRUE, 
			na.strings=c("[Not Available]", "[Not Reported]", "[Not Applicable]", "<-"))
	idx <- which(!( apply(fileContents, 2, function(x) sum(!is.na(x)) < 1) | apply(fileContents, 2, function(x) length(unique(x))==1) ))
	fileContents <- fileContents[,idx]
}

