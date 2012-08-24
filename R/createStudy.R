################################################################################
# Start by checking if the study of the same name already exists.  
# If it does, then we get it from Synapse and check its last update date
# to determine if it has been updated in its external location since the
# last time it was contributed to Synapse. If the update date is the same, 
# then we return the study.  If the last time it was
# updated extenally is more recent, then we update the study.  
#
#  If it does not yet exist in Synapse then we create it.
#################################################################################

setMethod(
		f = "createStudy",
		signature = signature("list","character"),
		definition = function(contribution, logFile){
			
			#
			# Remove any data annotations / properties from contribution
			#
			if(any(grepl("data",names(contribution)))){	
				contribution <- contribution[-grep("data",names(contribution))]
			}
			
			# Build the Study locally, then put the createEntity into a try statement.
			# if it fails, 
			names(contribution) <- gsub("study\\.","",names(contribution))
			
			# Create a new study
			splitBy <- splitAttributesAndProperties(contribution)
			if(class(splitBy) == "try-error") {
				return(splitBy)
			}
			
			# make sure there's a lastupdate date. If not, set it to today.  Then build
			# the study, add the annotations, and create it.
			if( ! ("lastUpdate" %in% names(contribution)) ) {
				splitBy$annotations$lastUpdate = format(Sys.time(), "%Y-%m-%d")
			}
			study <- Study(splitBy$properties)
			annotationValues(study) <- splitBy$annotations			
			study <- try(createEntity(study),silent=TRUE)
			
			if(class(study) == "try-error"){
				# study already exists
				qryResult <- .queryStudyExistence(contribution)
				# snippet below useful for running in workflow mode.
				if(!("study.id" %in% names(qryResult))){
					for(index in 1:10){
						cat("Waiting for annotations to be set properly!\n")
						qryResult <- .queryStudyExistence(contribution)
						if("study.id" %in% names(qryResult)){
							break;
						}
						Sys.sleep(0.2)
					}
				}
				# We couldn't get the query to successfully return after 11 tries.
				# throw an error and return
				if(!("study.id" %in% names(qryResult))){
					res <- "Cannot access study through Synapse API\n"
					f <- basename(tempfile())
					save(study, contribution, qryResult, file=f)
					class(res) <- 'try-error'
					return(res)
				}
				# We got the query result. Get the entity and return
				names(qryResult) <- tolower(names(qryResult))
				study <- getEntity(qryResult$study.id[1])
				return(study)
			}else if(class(study) == "Study"){
				# we created the study 
				cat("Created study:",contribution$name,"\n")
				return(study)
			}
			return(study)
		}
)

.queryStudyExistence <- function(contribution){
	qryString <- sprintf("select id, name, lastUpdate from study where study.parentId == \"%s\" and study.name == \"%s\"", 
			contribution$parentId, contribution$name)
	qryResult <- synapseQuery(qryString)
	return(qryResult)
}

.checkStudyInSynapse <- 
		function(contribution, qryResult) {
	# study exists.  Get entity from synapse
	names(qryResult) <- tolower(names(qryResult))
	study <- getEntity(qryResult$study.id[1])
	if("study.lastupdate" %in% names(qryResult)){
		# Determine if study has been updated since we stored it in Synapse:
		if (class(try(as.Date(contribution$lastUpdate),silent=TRUE)) == "Date" & qryResult$study.lastupdate != "Not Available") {
			time.since <- julian(as.Date(contribution$lastUpdate))[1] - julian(as.Date(qryResult$study.lastupdate))[1]
		}else{
			# Some contributions don't have a lastupdate value.  In this case we
			# set the last update date to the date at which the data was initially contributed.
			# When we run into such a study here we assume no change has been made since contribution.
			time.since=0
		}
	}else{
		# Set time.since to a positive number.  This forces us to update the existing study to add the lastUpdate annotation.
		time.since=1
	}
	if(time.since < 0){ 
		# This situation should not arise and if it does we want to be made aware of the issue.
		var <- "Error: study last update date in Synapse is later than study last update date from contribution."
		class(var) <- "try-error"
		return(var)
	}else if(time.since > 0){ 
		# Data set has been updated since it was last contributed so we want to update the study in synapse.
		# Note this replaces ALL existing annotations. 
		return(study)
	}else{
		# study is up to date.  Return study and continue on
		cat("study Exists...\n")
		return(study)
	}
}

