

contributeData <- function(contribution, project=NULL, keepLocal=FALSE, logFile=TRUE){

	if(logFile == "TRUE") {
		logFile <- .setUpLogFile(logFile)
		cat("log file for this workflow: ", logFile, "\n")
	}

	project <- handleProject(project)

	translatedContribution <- handleContribution(project, contribution, logFile)

	if(class(translatedContribution[[1]]) == "character"){
	  data <- contribute(translatedContribution, keepLocal,logFile)
	}else{ 
	  data <- sapply(translatedContribution, function(x){ contribute(x,keepLocal,logFile) })
	}
	return(data)
}


