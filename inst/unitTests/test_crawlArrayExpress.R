#
# This tests the array express crawler.  We reset the functions that makes calls out to the 
# array express api.  Here they simply load in a small version of the downloaded objects.
# The test then runs the crawlArrayExpress() function and checks to see if we get a specific
# answer.
#



.setUp <-
		function()
{
	
	myGAEJ <- 
			function()
	{
		load("../extdata/testData/arrayExpressJSON.Rda")
		return(arrayExpressJSON)
	}
	
	myPAEJ <- 
			function(arrayExpressJSON, aE.platformMap, project)
	{
		load("../extdata/testData/parsedArrayExpressJSON.Rda")
	}
	
	myGetURL <-
			function(parsedArrayExpressJSON, ...)
	{
		load("../extdata/testData/arrayExpressURLs.Rda")
		return(arrayExpressURLs)		
	}

	
	
	unloadNamespace("rSCR")
	assignInNamespace(".getArrayExpressJSON", myGAEJ, "rSCR")
	#assignInNamespace(".parseAEJsonObject", myPAEJ, "rSCR")
	assignInNamespace(".getArrayExpressRawDataURLs", myGetURL, "rSCR")
	attachNamespace("rSCR")
}

.tearDown <-
		function()
{
	unloadNamespace("rSCR")
	unloadNamespace("synapseClient")
	unloadNamespace("RCurl")
	attachNamespace("RCurl")
	attachNamespace("synapseClient")
	attachNamespace("rSCR")
}



unitTestCrawlArrayExpress <-
		function()
{
#	testProject <-	new(Class="Project")
#	propertyValue(testProject, "id") <- "12345"
#	ans <- crawlArrayExpress(testProject)
	## replace these with actual tests to make sure crawlArrayExpress works
#	load("../extdata/testData/arrayExpress.Rda")	
#	checkEquals(ans, arrayExpress)
}

