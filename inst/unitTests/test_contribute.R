#
# This test the array express crawler.  We reset the functions that hits the 
# array express api.  Here they simply load in a small version of the downloaded objects.
# the test then runs the crawlArrayExpress() function and checks to see 
#

.setUp <-
		function()
{
	
	myCD <- 
			function(contribution)
	{
		load(file.path(.path.package("rSCR"), "extdata/testData/egDataset.Rda"))
		return(dataset)
	}
	myQR <- 
			function()
	{
		load(file.path(.path.package("rSCR"), "extdata/testData/qryResult.Rda"))
		return(qryResult)
	}
	myCL <- 
			function(contribution)
	{
		load(file.path(.path.package("rSCR"), "extdata/testData/egLayer.Rda"))
		return(layer)
	}
	
	unloadNamespace("rSCR")
	assignInNamespace("createStudy", myCD, "rSCR")
	assignInNamespace(".checkLayerExistence", myQR , "rSCR")
	assignInNamespace("createLayer", myCL , "rSCR")
	attachNamespace("rSCR")
}

.tearDown <-
		function()
{
	unloadNamespace("rSCR")
	attachNamespace("rSCR")
}

unitTestcontribute <-
		function()
{

#	load(file.path(.path.package("rSCR"), "extdata/testData/arrayExpress.Rda"))
#	contribution <-  arrayExpress[[1]]
#	contribution$parentid <- "115399"
#	contribution$rawDataAvailable <- FALSE
	# Test situation where raw data is not available
#	checkTrue(is.null(contribute(contribution)))

}

