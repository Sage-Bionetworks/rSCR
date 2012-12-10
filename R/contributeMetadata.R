




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

syn341502

parentId <- 'syn1533995'
parent <- getEntity(parentId)

ent <- Data(list(name="GSE14996_metadata.tsv",parentId=parentId))
annotValue(ent, 'fileSize') <- .prettifyFileSize(file.info('GSE14996_metadata.tsv')$size)
ent <- addFile(ent, 'GSE14996_metadata.tsv')
annotValue(ent, 'study') <- propertyValue(parent,'name')
annotValue(ent, 'repository') <- annotValue(parent,'repository')
annotValue(ent, 'status') <- 'metadata'


ent <- loadEntity('syn1125641')
dat <- exprs(ent$objects$eset)
rownames(dat) <- gsub("_mt", "_eg", rownames(dat))
dat <- round(dat,3)
write.table(dat,file=paste(propertyValue(ent, 'name'),".tsv",sep=""),sep="\t",quote=FALSE)

newEnt <- Data(list(name=paste(propertyValue(ent, 'name'),'tsv',sep="."), parentId="syn1533995"))
propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
propertyValue(newEnt, 'platform') <- propertyValue(ent, 'platform')
propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
propertyValue(newEnt, 'numSamples') <- propertyValue(ent, 'numSamples')
annotations(newEnt) <- annotations(ent)
newEnt <- addFile(newEnt,'GSE14996_pd.genomewidesnp.6_gene.tsv')
annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(paste(propertyValue(ent, 'name'),'tsv',sep="."))$size)
annotValue(newEnt, 'status') <- 'processed'

