updateScrWithTcgaCrawlerOutput <- function(tcgaCrawlerOutput, private=FALSE, update=TRUE) {
	parentId1='syn301192'
	parentId2='syn150935'
	if(isTRUE(private)){ 
		parentId1='syn310038'
		parentId2='syn308810'
	}
	# Adds new TCGA layers to the SCR
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="tcgaCrawlerOutput" and entity.parentId=="', parentId1, '"', sep=""))
	if(!is.null(qry)){
		tcgaCrawlerEntity <- getEntity(qry$entity.id)
	}else{
		tcgaCrawlerEntity <- Data(list(name = "tcgaCrawlerOutput", parentId = parentId1, type="M"))
	}
	
	fname <- paste(tempdir(), 'tcgaCrawlerOutput.txt',sep='/')
	write.table(tcgaCrawlerOutput, file=fname, sep="\t", quote=FALSE,row.names=FALSE)
	tcgaCrawlerEntity <- addFile(tcgaCrawlerEntity,fname)
	tcgaCrawlerEntity <- addObject(tcgaCrawlerEntity,tcgaCrawlerOutput)
	if(is.null(qry)){
		tcgaCrawlerEntity <- createEntity(tcgaCrawlerEntity)		
	}else{
		tcgaCrawlerEntity <- updateEntity(tcgaCrawlerEntity)
	}
	tcgaCrawlerEntity <- storeEntity(tcgaCrawlerEntity)
	tcgaQry <- synapseQuery(paste('select id, name, acronym from study where study.repository == "TCGA" and study.parentId == "',parentId2,'"',sep=""))
	apply(tcgaQry, 1, function(x){
				cat(paste(x,collapse="\t"),"\n")
				qry <- paste('select id,name from entity where entity.parentId == "', x[2],'"',sep="")
				cancer <- synapseQuery(qry)
				intersect(tcgaCrawlerOutput$data.name, cancer$entity.name)
			}) -> existingEntities
	newEntities <- setdiff(tcgaCrawlerOutput$data.name, unlist(existingEntities))
	tcgaCrawlerOutput <- tcgaCrawlerOutput[!duplicated(tcgaCrawlerOutput$data.name),]
	rownames(tcgaCrawlerOutput) <- tcgaCrawlerOutput$data.name
	ent <- "No new entities to process!"
	if(length(newEntities) > 0){
		toUpdate <- tcgaCrawlerOutput[newEntities,]
		if(isTRUE(private)){
			toUpdate <- toUpdate[which(toUpdate$data.tcgaLevel!="clinical"),]
			toUpdate <- toUpdate[which(toUpdate$data.type!="C"),]
		}
		if(nrow(toUpdate)> 0){
			toUpdate$parentId <- parentId2
			if(isTRUE(update)){
				ent <- contributeData(toUpdate)
			}else{
				ent <- toUpdate
			}
		}else{
			ent <- "No new entities to process!"
		}
	}else{
		ent <- "No new entities to process!"
	}	
	return(ent)
}
