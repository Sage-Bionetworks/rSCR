

unitTesthandleContribution <-
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
	
	load(file.path(.path.package("rSCR"), "extdata/testData/egProject.Rda"))
		
	res <- handleContribution(project, contribution, "me")
	checkEquals(class(res),"list")
	# Test that contributions without any of url, name, and type are invalid
	tmp <- contribution[-which(names(contribution) == "layer.url")]
	res <- try(handleContribution(project, tmp, "me"),silent=TRUE)
	checkEquals(class(res),"try-error")
	
	tmp <- contribution[-which(names(contribution) == "layer.name")]
	res <- try(handleContribution(project, tmp, "me"),silent=TRUE)
	checkEquals(class(res),"try-error")
	
	tmp <- contribution[-which(names(contribution) == "layer.type")]
	res <- try(handleContribution(project, tmp, "me"),silent=TRUE)
	checkEquals(class(res),"try-error")
	
}

