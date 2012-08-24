#
#

.setUp <-
		function()
{
	
	myQD <- 
			function(contribution)
	{
		if(contribution$qry==2){
			return(NULL)
		}
		qryResult <- list(study.id="12345", study.lastUpdate="2004-02-11")
		return(qryResult)
	}
	myGD <- 
			function(x)
	{
		load(file.path(.path.package("rSCR"), "extdata/testData/egStudy.Rda"))
		return(dataset)
	}
	myCE <- 
			function(x)
	{
		return(x)
	}
	
	unloadNamespace("rSCR")
	assignInNamespace(".queryDatasetExistence", myQD, "rSCR")
	unloadNamespace("synapseClient")
	assignInNamespace("getEntity",    myGD , "synapseClient")
	assignInNamespace("createEntity", myCE , "synapseClient")
	assignInNamespace("updateEntity", myCE , "synapseClient")
	attachNamespace("synapseClient")
	attachNamespace("rSCR")
}

.tearDown <-
		function()
{
	unloadNamespace("rSCR")
	unloadNamespace("synapseClient")
	attachNamespace("synapseClient")
	attachNamespace("rSCR")
}




unitTestcreateStudy <-
		function()
{
	
	contribution <- list(study.name = "GSE7765_eg",
			species = "Homo sapiens",
			description = "MCF7 cells were treated with DMSO or 100 nM Dioxin for 16 hr. Gene expression changes were quantified by microarray ...",
			numSamples = 12,
			platform = "hgu133a;hgu133b",
			cellLine = "MCF7",
			parentId = "syn299343",
			data.url = "ftp://anonymous:anonymous@ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE7765/GSE7765_RAW.tar",
			data.status = "raw",
			data.type = "E",
			data.lastUpdate = "2007/05/12",                     
			data.name = "GSE7765 Raw Data Layer from NCBI GEO",
			data.compound = "DMSO,Dixoin")
	contribution$qry = 2

	res <- createStudy(contribution,"eg.txt")		
	checkEquals(class(res)[1], "Study")
	
	contribution$qry = 3	
	res <- createStudy(contribution, "eg.txt")		
	checkEquals(class(res)[1], "Study")
	
	contribution$study.lastUpdate <- "2004-02-12"
	res <- createStudy(contribution, "eg.txt")		
	checkEquals(class(res)[1], "Study")	
	
}
