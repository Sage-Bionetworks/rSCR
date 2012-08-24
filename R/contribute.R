####################################################################################
#
#  These methods control how various types of input are contributed to synapse.
#
####################################################################################


setMethod(
		f = "contribute",
		signature = c("list","logical"),
		definition = function(contribution, keepLocal, logFile){
			
			###############################################################
			# Get or build dataset from contribution
			###############################################################
			cat("Creating study...")
			study <- createStudy(contribution,logFile)
			if(class(study) == "try-error"){
				cat(study,"\n");
				return(study)
			}else{
				cat("... success!\n")
			}
			
			url <- contribution[which(names(contribution) == "data.url")]
			if(!grepl("\\.tar$|\\.gz$|\\.zip$",url)) {
				# Check that user provided url points to a .tar, .gz, or .zip file
				var <- "Raw data cannot be built as data.url does not describe a .tar, .gz, or .zip file\n"
				class(var) <- "try-error"
				return(var)
			}

			# Switch parentId to the datasets synapse Id so we can correctly
			# place the data in synapse.
			contribution$parentId <- propertyValue(study,"id")
			
			###############################################################
			# Try to contribute data to Synapse.
			###############################################################
			cat("Creating data...")
			data <- createData(contribution, keepLocal)
			if(class(data) != "try-error"){
				cat("... success!\n")
			}
			return(data)
		}
)


