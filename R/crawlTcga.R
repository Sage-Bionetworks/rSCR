

crawlTcga <- function(user='anonymous',pwd='anonymous',private=FALSE, debug=FALSE){
	# Define the header
	header <- c("study.name\tacronym\ttissueType\trepository\tdisease\tspecies\tstudy.numPatients\tdata.url\tdata.name\tdata.md5\tdata.type\tdata.platform\tdata.status\tdata.add\tdata.lastUpdate\tdata.numSamples\tdata.tcgaLevel")
  # Write out the header and define the base url for the public or private sites.
	if(isTRUE(private)){
		url <- "https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/"
		write.table(header, file="tcgaPrivate.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
	}else{
		url <- "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"
		write.table(header, file="tcga.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)
	}
	# Crawl the site
	links <- .getLinks(url,debug=debug,user=user,pwd=pwd,clinicalOnly=clinicalOnly,private=private)
	# Write out and store the crawled content.
	if(isTRUE(private)){
		tcga <- read.table("tcgaPrivate.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
		res <- .storetcgaInSynapse(tcga, private, fname='tcgaPrivate.txt')
	}else{
		tcga <- read.table("tcga.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
		res <- .storetcgaInSynapse(tcga, private, fname='tcga.txt')
	}
	return(tcga)
}

.storetcgaInSynapse <- function(tcga, private=FALSE, fname='tcga.txt'){ 
	# Adds new TCGA layers to the SCR
	entName <- 'tcgaPublicCrawlerOutput'
	if(private=="TRUE"){
		entName <- 'tcgaPrivateCrawlerOutput'
	}
	qry <- synapseQuery(paste('select id, name from entity where entity.name=="',entName,'" and entity.parentId=="syn1452692"', sep=""))
	if(!is.null(qry)){
		tcgaCrawlerEntity <- getEntity(qry$entity.id)
	}else{
		tcgaCrawlerEntity <- Data(list(name = entName, parentId = 'syn1452692'))
	}
	tcgaCrawlerEntity <- addFile(tcgaCrawlerEntity,fname)
	tcgaCrawlerEntity <- addObject(tcgaCrawlerEntity,tcga)
	if(is.null(qry)){
		tcgaCrawlerEntity <- createEntity(tcgaCrawlerEntity)		
	}else{
		tcgaCrawlerEntity <- updateEntity(tcgaCrawlerEntity)
	}
	tcgaCrawlerEntity <- storeEntity(tcgaCrawlerEntity)
}



.getURL <- function(url, curlHandle){
	# getURL robust to server side errors.
	response <- "error"
	class(response) <- 'try-error'
	i <- 1
	while(class(response) == "try-error"){
		if(i > 1){ cat("\rattempt",i) }
		response <- getURL(url, curl = curlHandle)
	}
	return(response)
}

.loadFiles <- function(){
	tcga.input.file <- paste(file.path(.path.package("rSCR"),"extdata"),"/tcgaAnnotations.txt",sep="")
	tcgaAnnotations<- read.table(tcga.input.file,stringsAsFactors=FALSE, sep="\t", header=TRUE,quote="")
	k <- ncol(tcgaAnnotations)
	tcgaAnnotations <- tcgaAnnotations[,-k]
	rownames(tcgaAnnotations) <- tolower(tcgaAnnotations[,2])
	sapply(rownames(tcgaAnnotations), function(x){ 
				paste(tcgaAnnotations[x,],collapse="\t")
			}) -> res
	return(res)
}

.getLinks <- function(url,debug=FALSE,user="anonymous",pwd="anonymous",clinicalOnly=FALSE,private=FALSE){
# Load TCGA Annotations
	tcgaAnnotations <- .loadFiles();
# Get the page
	cat(url,"\n")	
	if(isTRUE(debug)){
		if(user=='anonymous'){
			if(!grepl("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/gbm",url) & url != "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/"){
				return(NA)
			}
		}else{
			if(!grepl("https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/gbm",url) & url != "https://tcga-data-secure.nci.nih.gov/tcgafiles/tcga4yeo/tumor/"){
				return(NA)
			}
		}
	}
	h = getCurlHandle(header = FALSE, userpwd = paste(user,pwd,sep=":"), netrc = TRUE)
	text <- .getURL(url, curl = h)
	# Find all links
	links.and.dates <- regmatches(text,
			gregexpr('<a href=\\\"[^"]+\\\">[^<]+<\\/a>\\s+\\d+-[^-]+-\\d+',text))[[1]]	
	links <- unlist(regmatches(links.and.dates,
					gregexpr('<a href=\\\"[^"]+\\\">[^<]+<',links.and.dates)))
	# Drop the ones we don't care about
	drop <- grepl("abi",links) + grepl("pathology",links) + grepl('Name<',links) + grepl('microsat_i',links) + grepl('Last modified<',links) + grepl('Size<',links) + grepl('Parent',links) + grepl('diagnostic_images',links) + grepl("stanford.edu",links)
 # Get the ones we want to follow or push to Synapse
	links <- links[drop == 0]
	if(isTRUE(clinicalOnly)){
		drop <- grepl('cgcc',links) + grepl('gsc', links)
		links <- links[drop == 0]
	}
	if(length(links)==0){ return(url) }
	# Remove the carrots from the text
	dates <- as.Date(unlist(regmatches(links.and.dates,
							gregexpr('\\d+-[^-]+-\\d+',links.and.dates))),"%d-%b-%Y")
	dates <- dates[drop == 0]
	sapply(links, function(link){
				l <- gsub("[><]","",regmatches(link,gregexpr('>[^<]+<',link))[[1]])
				l
			}) -> cleanedLinks
	links <- as.character(cleanedLinks)
	if(!grepl('bcr', url) & sum(grepl('tar.gz',links) > 0) & isTRUE(!clinicalOnly)) {
		# There's some tar.gz files and we aren't only updating the clinical portion
		links <- .handleDataLinks(url,links,tcgaAnnotations,dates,user=user,pwd=pwd,private=private)
		return(links)
	}else if(any(grepl("clinical_", links) | grepl("slide_images",url) | grepl("\\.bio\\.L",links))){
		links <- .handleClinicalLinks(url,links,tcgaAnnotations,dates,private=private,user=user, pwd=pwd)
		return(links)
	}else{
		if(any(grepl("IlluminaGA_mRNA_DGE",links))){
			return(links)
		}
		for( i in 1:length(links) ){
			next.link <- paste(url, links[i], sep="")
			res <- .getLinks(next.link,debug=debug,user=user,pwd=pwd, clinicalOnly=clinicalOnly, private=private)
		}
	}
	return(res)
}

.handleDataLinks <- function(url, links,tcgaAnnotations, dates,user="anonymous",pwd="anonymous", private=FALSE){
	platforms <- .loadPlatformKey()
	types <- .loadTypeKey()
	split.urls <- strsplit(url,'/')[[1]]
	platform <- split.urls[length(split.urls)-1]
	type <- split.urls[length(split.urls)]
	platform <- platforms[platform]
	type <- types[type]
	drop <- grepl("Level",links) + grepl("mage-tab",links)
	links <- links[which(drop > 0)]
	dates <- dates[which(drop > 0)]
	cancer.type <- regmatches(url, gregexpr("tumor\\/[^\\/]+",url))[[1]]
	cancer.type <- gsub("tumor/","",cancer.type)
	if(any(grepl(".tar.gz$",links,perl=TRUE))) {
		ids <- grep("tar.gz$",links)
		tar.gzs <- links[grepl(".tar.gz$",links,perl=TRUE)]		
		sapply(ids, function(x){ 
					link <- links[x]
					url2 <- paste(url,gsub(".tar.gz","",link),"/",sep="")
					h = getCurlHandle(header = FALSE, userpwd = paste(user,pwd,sep=":"), netrc = TRUE)
					text <- .getURL(url2, curl = h)
					if(grepl("HTTP 404",text)){
						num.files <- "NA"
					}else{
						links.and.dates <- regmatches(text,gregexpr('<a href=\\\"[^"]+\\\">[^<]+<\\/a>\\s+\\d+-[^-]+-\\d+',text))[[1]]
						links <- unlist(regmatches(links.and.dates,gregexpr('<a href=\\\"[^"]+\\\">[^<]+<',links.and.dates)))
						drop <- grepl("DESCRIPTION",links)  + grepl("MANIFEST",links) + grepl("README",links)
						links <- links[which(drop==0)]
						num.files <- length(links)
					}
					date <- dates[x]
					md5Link <- paste(c(url,link,".md5"),collapse="")
					h = getCurlHandle(header = FALSE, userpwd = paste(user,pwd,sep=":"), netrc = TRUE)
					text <- .getURL(md5Link, curl = h )
					md5 <- strsplit(text,"\\s+",perl=TRUE)[[1]]
					level <- .getTcgaLevel(md5[2])
					if(grepl("mage",link)){
						type <- 'C'
						dataLink <- paste(tcgaAnnotations[[cancer.type]],paste(c(url,md5[2]),collapse=""), md5[2], 'NA',type, platform, "raw", "TRUE",date,num.files,level,sep="\t")
					}else{
						dataLink <- paste(tcgaAnnotations[[cancer.type]],paste(c(url,md5[2]),collapse=""), md5[2], md5[1],type, platform, "raw", "FALSE",date,num.files,level,sep="\t")
					}
					if(isTRUE(private)){
						write.table(dataLink, file="tcgaPrivate.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
					}else{
						write.table(dataLink, file="tcga.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
					}
					
					c(paste(c(url,md5[2]),collapse=""),md5[1])
				}) -> dataLinks
	}else{
		dataLinks <- NA
	}
	return(dataLinks)
}	

.getTcgaLevel <- function(dataName){
	if(grepl('Level_1',dataName)){
		return("level_1")
	}else if(grepl('Level_2',dataName)){
		return("level_2")
	}else if(grepl('Level_3',dataName)){
		return("level_3")
	}else{
		return("clinical")
	}
}

.handleClinicalLinks <- function(url,links,tcgaAnnotations,dates,private,user=user,pwd=pwd) {	
	platforms <- .loadPlatformKey()
	types <- .loadTypeKey()
	split.urls <- strsplit(url,'/')[[1]]
	platform <- split.urls[length(split.urls)-1]
	type <- split.urls[length(split.urls)]
	platform <- as.character(unlist(platforms[platform]))
	type <- as.character(unlist(types[type]))
	
	cancer.type <- regmatches(url, gregexpr("tumor\\/[^\\/]+",url))[[1]]
	cancer.type <- gsub("tumor/","",cancer.type)
	id <- c(grep(".tar.gz$",links,perl=TRUE),grep(".tar$",links,perl=TRUE))
	
	for(i in 1:length(id)){
		dataAdd <- 'TRUE'
		tar.gz <- links[id[i]]
		date <- as.character(unlist(dates[id[i]]))
		dataLink <- paste(c(url,tar.gz),collapse="")
		#download <- try(synapseClient:::.curlWriterDownload(dataLink,"eme.tar.gz"))
		md5 <- NA # as.character(tools::md5sum("eme.tar.gz"))
		if(grepl('tissue_images', dataLink)){
			md5Link <- paste(c(url,tar.gz,".md5"),collapse="")
			h = getCurlHandle(header = FALSE, userpwd = paste(user,pwd,sep=":"), netrc = TRUE)
			text <- .getURL(md5Link, curl = h )
			md5 <- strsplit(text,"\\s+",perl=TRUE)[[1]][1]
			dataAdd <- 'FALSE'
		}
		num.files <- "NA"
		level <- .getTcgaLevel(tar.gz)
		if(length(platform)==0){
			platform <- "NA"
		}
		
		output <- paste(c(tcgaAnnotations[[cancer.type]], dataLink, tar.gz, md5,type,platform, "raw",dataAdd,date,num.files,level),collapse="\t")	
		if(isTRUE(private)){
			write.table(output, file="tcgaPrivate.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}else{
			write.table(output, file="tcga.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
	return(tar.gz)
}

.loadTypeKey<- function(){
	types <- list(
			RNASeq = 'E',
			RNASeqV2 = 'E',
			rnaseq=	'E',
			miRNASeq=	'E',
			transcriptome	='E',
			fragment_analysis	=NA,
			protein_exp	='P',
			mutations ='G',
			mirnaseq ='E',
			mirna	='E',
			snp	='G',
			exon ='E',
			slide_images='M',
			clin ='C',
			methylation	= 'G',
			cna	= 'G')
	return(types)
}

.loadPlatformKey <- function(){
	platforms <- list(
			'agilentg4502a_07_3'= "agilentg4502a_07",
			'illuminahiseq_mirnaseq'= "illuminahiseq_mirnaseq",
			'humanmethylation450'= "humanmethylation450",
			'microsat_i'= "microsat_i",
			'mda_rppa_core'= "mda_rppa_core",
			'humanmethylation27'= "humanmethylation27",
			'illuminadnamethylation_oma003_cpi'= "illuminadnamethylation",
			'human1mduo'= "human1mduo",
			'hg-u133_plus_2'= "hgu133plus2",
			'illuminaga_mirnaseq'= "illuminaga_mirnaseq",
			'illuminaga_dnaseq'= "illuminaga_dnaseq",
			'agilentg4502a_07_1'= "agilentg4502a_07",
			'huex-1_0-st-v2'= "huex10stv2",
			'illuminadnamethylation_oma002_cpi'= "illuminadnamethylation",
			'bcgsc.ca'= "bcgsc.ca",
			'illuminaga_rnaseq'= "illuminaga_rnaseq",
			'h-mirna_8x15k'= "h-mirna_8x15k",
			'agilentg4502a_07_2'= "agilentg4502a_07",
			'h-mirna_8x15kv2'= "h-mirna_8x15kv2",
			'hg-cgh-244a'= "hg-cgh-244a",
			'minbiotab'= "minbiotab",
			"mda_rppa_core" = "mda_rppa_core",
			'genome_wide_snp_6'= "pd.genomewidesnp.6",
			'humanhap550'= "humanhap550",
			'hg-cgh-415k_g4124a'= "hg-cgh-415k_g4124a",
			'solid_dnaseq'= "solid_dnaseq",
			'cgh-1x1m_g4447a'= "cgh-1x1m_g4447a",
			'illuminahiseq_dnaseqc'= "illuminahiseq_dnaseqc",
			'ht_hg-u133a'= "hthgu133a",
			'tissue_images' = 'tissue_images',
			'illuminahiseq_rnaseq'= "illuminahiseq_rnaseq",
			'illuminahiseq_rnaseqv2'= "illuminahiseq_rnaseqv2")
	return(platforms)
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
