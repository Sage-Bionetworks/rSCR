#
# This tests the createLayer function.  We reset the functions that make calls out to 
# Synapse.  We play with the example to test the logic around updating an existing data
# set depending upon when it was most recently updated externally and in Synapse.
#

.setUp <-
		function()
{
	
	myQL <- 
			function(contribution)
	{
		if(contribution$qry==2){
			return(NULL)
		}
		qryResult <- list(layer.lastUpdate = "2004-02-11", layer.id = "1234")
		return(qryResult)
	}
	myGL <- 
			function(x)
	{
		load(file.path(.path.package("rSCR"), "extdata/testData/egLayer.Rda"))
		return(layer)
	}
	myCE <- 
			function(x)
	{
		return(x)
	}	
	myHL <- 
			function(x,y)
	{
		return(x)
	}
	
	unloadNamespace("rSCR")
	assignInNamespace(".checkLayerExistence", myQL, "rSCR")
	assignInNamespace("handleLocation", myHL, "rSCR")
	unloadNamespace("synapseClient")
	assignInNamespace("getEntity", myGL , "synapseClient")
	assignInNamespace("createEntity", myCE , "synapseClient")
	assignInNamespace("updateEntity", myCE , "synapseClient")
	attachNamespace("synapseClient")
	attachNamespace("rSCR")
}

.tearDown <-
		function()
{
	unloadNamespace("rSCR")
	attachNamespace("rSCR")
}


unitTestcreateLayer <-
		function()
{
	
	contribution <-
			list(
					dataset.name="GSE30884_test2",
					dataset.url = 'http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE30884',
					dataset.contact='someone@somewhere',
					layer.name="GSE30884_rawData",
					layer.status="raw",
					layer.type="E",
					parentId = "102610",
					Species="Homo sapiens",
					rawDataAvailable=TRUE,
					layer.lastUpdate="2004-02-12",
					layer.url="ftp://anonymous:anonymous@ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE30884/GSE30884_RAW.tar");
	
	contribution$qry = 3
	contribution$checksum = "59d01989613bb2883c1d19ac6caa9a3e"	
	contribution$layerLocation<- list(layerLocation=list(path="ftp://anonymous:anonymous@ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/series/GSE30884/GSE30884_RAW.tar"),
			type="external")

	# Make certain we get through function call and return object of class ExpressionLayer
	res <- createLayer(contribution, FALSE)
	checkEquals(class(res),"ExpressionLayer")

	# Should throw error as lastUpdate date is not valid date 
	contribution$layer.lastUpdate <- "hello"
	res <- createLayer(contribution, FALSE)
	checkEquals(class(res),"try-error")
	
	# Should throw error as lastUpdate date suggests data was last updated externally prior to most recent synapse update 
	contribution$layer.lastUpdate <- "2004-02-10"
	res <- createLayer(contribution, FALSE)
	checkEquals(class(res),"try-error")

	# Should throw error as lastUpdate date suggests data is up to date 
	contribution$layer.lastUpdate <- "2004-02-13"
	res <- createLayer(contribution, FALSE)
	checkEquals(class(res),"ExpressionLayer")

	#Now we try to build a layer when it doesn't already exist
	contribution$qry = 2
	res <- createLayer(contribution, FALSE)
	checkEquals(class(res),"ExpressionLayer")
			
}

