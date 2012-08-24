tryCatch.W.E <- function(expr, logFile, expectedClass=NULL, msg=NULL)
{
	msg <- paste("metaGenomics purpose:",msg)
	W <- NULL
	w.handler <- function(w){ # warning handler
		W <<- w
		invokeRestart("muffleWarning")
	}
	ret <- list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), 
					warning = w.handler),
			warning = W)
	if( !is.null(ret$warning) ){
		# value contains the error message, warning contains the warning message.
		if(class(ret$value) != expectedClass){
			# If we get a warning, don't get an error, and don't return the correct class type, then write out to logFile and
			# return full object
			class(ret) <- "try-error"
			if( !is.null(msg)){ write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)}
			write.table(as.character(unlist(ret)),file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
			ret
		}else{
			# If we get a warning, but get the correct class type, just return the value
			if(is.list(ret$warning)) {
				if( !is.null(msg)){ write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)}
				write.table(as.character(unlist(ret$warning)),file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
			}else{
				if( !is.null(msg)){ write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)}
				write.table(as.character(ret$warning),file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
			}
			ret$value
		}
	}else if(!is.null(ret$value$message)){
		# An error occurred
		# value contains the error message, warning contains the warning message.
		class(ret) <- "try-error"
		if( !is.null(msg)){ write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)}
		write.table(as.character(unlist(ret)),file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE) 
		ret
	}else{
		if( !is.null(msg)){ write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)}
		ret$value
	}	
}

setUpLogFile <-  function(logFile) { 
	if(logFile=="TRUE"){ 
		logFile <- paste("mG_", paste(sample(seq(0,9),10,replace=TRUE),collapse=""),".log",sep="")
	}
	logFile
}

sendLogMessage <- function(msg,logFile){ 
	write.table(msg,file=logFile,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

