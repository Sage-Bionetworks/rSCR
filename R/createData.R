
setMethod(
		f = "createData",
		signature = c("list","logical"),
		definition = function(contribution,keepLocal){
			
			# Remove any dataset annotations or properties
			if(any(grepl("study",names(contribution)))){	
				contribution <- contribution[-grep("study",names(contribution))]
			}
			names(contribution) <- gsub("data\\.","",names(contribution))
			if(! all(c("type","status","url","name","parentId") %in% names(contribution))) {				
				val <- paste("All of data.name, parentId, data.url, data.status, and data.type are required.\n")
				class(val) <- "try-error"
				return(val)
			}
			
			###############################################################
			# Determine if data already exists... 
			###############################################################
			qryResult <- .checkDataExistence(contribution)
			
			##################################################################
			# If it does, compare the last update dates. If date is equal to #
			# current data in synapse then stop and send appropriate message#
			##################################################################
			if(!is.null(qryResult)){
				res <- .compareDataDates(contribution, qryResult)
				# For now we just skip anything that already exists.
				if(class(res) == "try-error") {
					return(res)
				}
			}
			#################################################################
			# The data either does not exist or it needs to be udpated.
			# First we handle the location information.
			#################################################################
			res <- handleLocation(contribution,keepLocal)
			if(class(res) == "try-error"){
				return(res)
			}
			contribution <- res$contribution
			destfile <- res$destfile
			
			#Create the data
			cat("Creating data", contribution$name,"...")
			tmp <- splitAttributesAndProperties(contribution)
			#Either add the downloaded file to the data or add the location object
			if('add' %in% names(contribution) & isTRUE(contribution$add)){
				data <- Data(tmp$properties)
				data <- addFile(data, destfile)
			}else{
				data <- Data(tmp$properties)
			}
			annotationValues(data) <- tmp$annotations
			
			# If data does not exist then create or store it.  Otherwise update the md5 sum and move on.
			if (is.null(qryResult)) {
				if('add' %in% names(contribution) & isTRUE(contribution$add)){					
					data <- storeEntity(data)
				}else{
					data <- createEntity(data)
				}
				cat("Data created!\n")
			}else{
				data <- getEntity(qryResult$entity.id)
				if('add' %in% names(contribution) & isTRUE(contribution$add)){
					data <- loadEntity(qryResult$entity.id)
					data <- deleteFile(data, data$files)
					annotValue(data, 'lastUpdate') <- contribution$lastUpdate
					data <- addFile(data, destfile)
					data <- updateEntity(data)
					data <- storeEntity(data)
				}else{
					propertyValue(data, "md5") <- contribution$checksum
					data <- updateEntity(data)
				}
				cat("Data updated!\n")
			}
			# If the file exists, then we want to remove it unless the user set keppLocal to TRUE.
			if(!isTRUE(keepLocal) & file.exists(destfile)){
				file.remove(destfile)
			}
			return(data) 
		}
)

.checkDataExistence <- function(contribution){
	qryString <- sprintf("select id, lastUpdate from entity where entity.parentId == \"%s\" and entity.name == \"%s\"", 
			contribution$parentId, contribution$name)
	qryResult <- synapseQuery(qryString)
}

.compareDataDates <- 
		function(contribution, qryResult){
	names(qryResult) <- tolower(names(qryResult))
	time.since=1
	# Determine if dataset has been updated in ArrayExpress since we stored it in Synapse:
	if (!is.na(qryResult$entity.lastupdate) & class(try(as.Date(contribution$lastUpdate),silent=TRUE)) == "Date" & class(try(as.Date(qryResult$entity.lastupdate),silent=TRUE)) == "Date" ) {
		# some datasets don't have a lastUpdate date.  In this case we
		# set the last update date to the date at which the data was initially contributed to Synapse.  
		time.since <- julian(as.Date(contribution$lastUpdate))[1] - julian(as.Date(qryResult$entity.lastupdate))[1]
	}else if ( class(try(as.Date(contribution$lastUpdate),silent=TRUE)) == "Date" ){
		# set time.since equal to a positive number so we know to update or create the data
		time.since = 1
	}else{
		# These contributions have an invalid last update date.
		var <- "Error: Last update date not valid."
		class(var) <- "try-error"
		cat(var)
		return(var)
	}
	if(time.since < 0){ 
		# This situation should not arise and if it does we want to be made aware of the issue.
		var <- "Error: Data last update date in Synapse is later than dataset last update date."
		class(var) <- "try-error"
		cat(var)
		return(var)
	}else if (time.since==0){
		# data is up to date
		var <- "Data up to date...\n"
		class(var) <- "try-error"
		cat(var)
		return(var)
	}
	# otherwise we need to update the data in synapse.
	return(1)
}

