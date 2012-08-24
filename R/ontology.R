# Build synonym tables
# Build controlled vocabulary tables

# Assume inputs come in two forms. 
# The first are tab-delimited files, with samples along the rows and metadata along the columns
# The second is the four column format, with study, sample, key, value
# So one function we need is a translator from first to second file type.

translateMatrixInto4Columns <- function(mat, studyName){
	sampleColumn <- which(colnames(mat) == 'sampleName')
	allKeys <- rep(colnames(mat[,-sampleColumn]), each=nrow(mat))
	allValues <- as.character(unlist(mat[,-sampleColumn]))
	allSamples <- rep(mat[,sampleColumn], times=ncol(mat)-1)
	allStudies <- rep(studyName, length(sampleNames))
	obj <- as.data.frame(cbind(allStudies, allSamples, allKeys, allValues))
	names(obj) <- c('studies', 'samples','keys','values')
	return(obj)
}

standardizeMatrix <- function(mat, key2value, syn2value){
	sampleColumn <- which(colnames(mat) == 'sampleName')
	colnames(mat) <- cleanTerms(colnames(mat))
	apply(mat, 2, function(x){ 
			cleanTerms(x)
		}) -> cleanedMat
	
	validKeys <- which(colnames(cleanedMat) %in% names(key2value))
	inValidKeys <- setdiff(1:ncol(mat), validKeys)
	unmappedValues <- matrix(,nc=2)
	for( i in 1:length(validKeys)){
		vk <- validKeys[i]
		validValues <- key2value[[colnames(mat)[vk]]]
		tmp <- match(cleanedMat[,vk], validValues)				
		inValidTerms <- which(is.na(tmp))
		validTerms <- which(!is.na(tmp))
		if(length(inValidTerms) == 0 | length(validTerms) == 0){
			next;
		}
		tmp <- match(cleanedMat[inValidTerms,vk], names(syn2value))
		still.inValidTerms <- inValidTerms[is.na(tmp)]
		cleanedMat[inValidTerms[!is.na(tmp)],vk] <- as.character(unlist(syn2value[tmp[!is.na(tmp)]]))
		key <- colnames(cleanedMat)[vk]
		tmp <- cbind(rep(key, length(still.inValidTerms)),
					cleanedMat[still.inValidTerms,vk])
		tmp <- tmp[!duplicated(tmp[,2]),]		
		cat(i,"\n")
		unmappedValues <- rbind(unmappedValues,tmp)
	}
	return(list(matrix=cleanedMat, unmapped=unmappedValues))
}

cleanTerms <- function(terms){
	terms <- tolower(gsub("[\\.\\s]+",'_', terms,perl=TRUE))
	terms <- gsub("[\\s]+$", "", terms, perl=TRUE)
	terms
}

loadKeyValues <- function(ontologyDir = '/Users/brig.mecham/Documents/workspace/metaGEO/trunk/geoAnnotation/Ontologies'){
	# These should be two column files.
	allKeys <- list.files(ontologyDir)
	allDirs <- list.dirs(ontologyDir)
	allDirs <- allDirs[!grepl("svn",allDirs)]
	allDirs <- allDirs[allDirs != ontologyDir]
	tmp <- read.table(list.files(allDirs[1], full.names=TRUE), stringsAsFactors=FALSE)[,1]
	keyValues <- cbind(rep(basename(allDirs[1]), length(tmp)), tmp)
	colnames(keyValues) <- c("key","value")
	for(i in 2:length(allDirs)){
		if(i == 36 | i == 60){ 
			next; 
		}
		x <- allDirs[i]
		cat(i,"\t",basename(x),"\n")
		fname <- paste(x, '/', basename(x),'_ont.txt',sep="")		
		tmp <- read.table(fname, stringsAsFactors=FALSE,sep="\t",quote="")[,1]
		L.keyValues <- cbind(rep(basename(x), length(tmp)), tmp)
		keyValues <- rbind(keyValues, L.keyValues)				
	}
	keyValues[,1] <- cleanTerms(keyValues[,1])
	keyValues[,2] <- cleanTerms(keyValues[,2])
	key2value <- split(keyValues[,2], keyValues[,1])
	return(key2value)
}

loadSynonyms <- function(synDir='/Users/brig.mecham/Documents/workspace/metaGEO/trunk/geoAnnotation/synonym_tables'){ 
  allSyns <- list.files(synDir, full.names=TRUE)
	allSyns <- allSyns[grepl('_syn.txt',allSyns)]
	tmp <- read.table(allSyns[1], stringsAsFactors=FALSE, sep="\t",quote="")
	synonyms <- cbind(rep(basename(allSyns[1]), nrow(tmp)), tmp)
	colnames(synonyms) <- c("key","synonym", "value")
	for(i in 2:length(allSyns)){
		if(i == 10){ 
			next; 
		}
		x <- allSyns[i]
		cat(i,"\t",basename(x),"\n")
		fname <- x		
		tmp <- read.table(fname, stringsAsFactors=FALSE,sep="\t",quote="")
		L.synonyms<- cbind(rep(basename(x), nrow(tmp)), tmp)
		colnames(L.synonyms) <- c("key","synonym", "value")
		synonyms <- rbind(synonyms, L.synonyms)				
	}
	synonyms[,2] <- cleanTerms(synonyms[,2])
	synonyms[,3] <- cleanTerms(synonyms[,3])
	syn2value <- split(synonyms[,3], synonyms[,2])
	return(syn2value)
}

#
key2value <- loadKeyValues()
syn2value <- loadSynonyms()
#
mat <- read.table( '~/metaGEO/GSE2990_metadata.xls',
		sep="\t",
		header=TRUE, 
		stringsAsFactors=FALSE, 
		quote="")
mat[sample(nrow(mat),10),2] <- 'adenocarcinoma-very_poorly_differentiated'
mat[sample(nrow(mat),10),2] <- 'adnocrcnma'
mat[sample(nrow(mat),10),4] <- 'madeUp'
res <- standardizeMatrix(mat,key2value, syn2value)

#
#createSQLiteOntology <- function(synapseId='syn310033') {
#	# This function downloads the synapse ontology and builds
#	# a SQLite database.  It returns a database connection.
#	m <- dbDriver("SQLite")
#	ontFile <- 'synapseOntology.sql'
#	con <- dbConnect(m, dbname = ontFile)
#	ent <- loadEntity(synapseId)
#	files <- paste(ent$cacheDir,'/',ent$files,sep="")
#	ontology <- read.table(files[1],sep="\t",header=FALSE, stringsAsFactors=FALSE)
#	names(ontology) <- c("key","value")
#	ontology$key <- tolower(ontology$key)
#	ontology$value <- tolower(ontology$value)
#	synonymKeys <- read.table(files[2],sep="\t",header=FALSE, stringsAsFactors=FALSE)
#	names(synonymKeys) <- c("synonym","key")
#	synonymKeys$synonym <- tolower(synonymKeys$synonym)
#	synonymKeys$key <- tolower(synonymKeys$key)
#	synonymValues <- read.table(files[3],sep="\t",header=FALSE, stringsAsFactors=FALSE)
#	names(synonymValues) <- c("key","synonym","value")
#	synonymValues$synonym <- tolower(synonymValues$synonym)
#	synonymValues$key <- tolower(synonymValues$key)
#	synonymValues$value <- tolower(synonymValues$value)
#	dbWriteTable(con, "ontology", ontology)
#	dbWriteTable(con, "synonymValues", synonymValues)
#	dbWriteTable(con, "synonymKeys", synonymKeys)
#	con	
#}
#
#con <- createSQLiteOntology()
## query
#rs <- dbSendQuery(con, 'select * from synonymValues where synonym=="gbm"')
#fetch(rs,-1)	
#rs <- dbSendQuery(con, 'select * from synonymKeys where synonym=="estrogen_receptor_status"')
#fetch(rs,-1)	
#
#setMethod(
#		f = "standardizeVocabulary",
#		signature = c("character","character","SQLiteConnection"),
#		definition = function(key, value, con){
#			rs <- dbSendQuery(con, paste('select * from ontology where key == "',key,'" and value=="',value,'"', sep=""))
#			allHits <- fetch(rs,-1)
#			if(nrow(allHits) == 1){ 
#			# found the match directly.. This key:value pair is already in the ontology.  
#				return(key,value)
#			}else if(nrow(allHits)> 1){
#				# found more than one match. Throw and error
#				var <- paste('found multiple matches for', key, value)
#				class(var) <- 'try-error'
#				return(var)
#			}else{
#				# Didn't find any matches.
#				
#			}
#				
#		})
#
#file.remove('synapseOntology.sql')
#
#d1 <- fetch(rs, n = 10) # extract data in chunks of 10 rows
#dbHasCompleted(rs)
#d2 <- fetch(rs, n = -1) # extract all remaining data
#dbHasCompleted(rs)
#dbClearResult(rs)
#dbListTables(con)
## clean up
#dbDisconnect(con)
#file.info(ontFile)
