#
# This test the array express crawler.  We reset the functions that make calls out the 
# array express api.  Here they simply load in a small version of the downloaded objects.
# the test then runs the crawlArrayExpress() function and checks to see 
#

.setUp <-
		function()
{
	
}

.tearDown <-
		function()
{
	unloadNamespace("rSCR")
	attachNamespace("rSCR")
}



unitTesthandleLocation <-
		function()
{
	
#	load(file.path(.path.package("rSCR"), "extdata/testData/contributionToLayer.Rda"))
#	checkTrue(class(handleLocation(contribution,FALSE))=="list")
	
	# Test situation where user isn't providing a layer name
#	load(file.path(.path.package("rSCR"), "extdata/testData/contributionToLayer.Rda"))
#	contribution <- contribution[-which(names(contribution)=="url")]
#	checkTrue(class(try(handleLocation(contribution,FALSE),silent=TRUE)) == "try-error")
		
}

