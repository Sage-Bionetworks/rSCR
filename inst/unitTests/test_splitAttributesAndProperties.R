
#.setUp <-
	#	function()#
#{
	
#}

#.tearDown <-
#		function()
#{
#	unloadNamespace("rSCR")
#	attachNamespace("rSCR")
#}



unitTestsplitAttributesAndPropeties <-
		function()
{
	
#	load(file.path(.path.package("rSCR"), "extdata/testData/contributionToLayer.Rda"))
#	names(contribution)[which(names(contribution)=="parentid")] <- "parentId"
	
#	checkTrue(class(splitAttributesAndPropeties(contribution)) == "list")
#	contribution <- contribution[-which(names(contribution)=="name")]
#	checkTrue(class(try(splitAttributesAndPropeties(contribution))) == "try-error")
	
}

