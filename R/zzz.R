# Lifted from zzz.R file in synapseClient package 1-9-12
.SCRcache <- new.env(parent=emptyenv())

## package-local 'getter'
.getCacheSCR <-
		function(key)
{
	.SCRcache[[key]]
}

## package-local 'setter'
.setCacheSCR <-
		function(key, value)
{
	.SCRcache[[key]] <- value
}
## delete cache
.deleteCacheSCR <-
		function(keys)
{
	indx <- which(keys %in% ls(.SCRcache))
	if(length(indx) > 0)
		rm(list=keys[indx], envir=.SCRcache)
}

#on load
.onLoad <-
		function(libname, pkgname)
{
	.setCacheSCR("arrayExpressJSONUrl","http://www.ebi.ac.uk/arrayexpress/json/v2/experiments?directsub=true&array=A-AFFY*")
}

#on unload
.onUnload <- function(libpath) .Last.lib()

#Last lib
.Last.lib <- function(...) {
	try(stoppedStep <- stopStep(), silent=TRUE)
}
