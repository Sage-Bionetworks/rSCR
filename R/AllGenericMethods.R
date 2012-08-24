setGeneric(
		name = "handleLocation",
		def = function(contribution,keepLocal){
			standardGeneric("handleLocation")
		}
)

setGeneric(
		name = "contribute",
		def = function(contribution, keepLocal,logFile){
			standardGeneric("contribute")
		}
)

setGeneric(
		name = "handleContribution",
		def = function(project, contribution, logFile){
			standardGeneric("handleContribution")
		}		
)

setGeneric(
		name = "handleProject",
		def = function(project, build,update){
			standardGeneric("handleProject")
		}		
)

setGeneric(
		name = "createData",
		def = function(contribution,keepLocal){
			standardGeneric("createData")
		}		
)

setGeneric(
		name = "createStudy",
		def = function(contribution,logFile){
			standardGeneric("createStudy")
		}		
)

setGeneric(
		name = "handleTcgaClinical",
		def = function(entity){
			standardGeneric("handleTcgaClinical")
		}		
)

setGeneric(
		name = "standardizeVocabulary",
		def = function(key, value, SQLiteConnection){
			standardGeneric("standardizeVocabulary")
		}		
)




