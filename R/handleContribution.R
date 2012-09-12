#
#  These functions handle the work associated with 
#  gathering the information about the contributed
#  dataset(s).
#

setMethod(
		f = "handleContribution",
		signature = signature("Project", "list","character"),
		definition = function(project, contribution, logFile) {
			# First remove any annotations set to NA
			if(sum(is.na(contribution))){
				contribution <- contribution[-which(is.na(contribution))]
			}
			# Next, make sure the user has provided the required information.
			if(! all(c("data.url","study.name","data.name", "data.type","data.status") %in% names(contribution))){				
				cat(setdiff(c("data.url","study.name","data.name", "data.type","data.status"), names(contribution)),"\n")
				stop("Mising one of study.name, data.url, data.name, data.type, or data.status.  Please try again.\n")
			}
			# Set the data.lastUpdate to NA if it wasn't provided
			if(! "data.lastUpdate" %in% names(contribution)){
				# Last update not included.  Set to NA
				contribution$data.lastUpdate = format(Sys.time(), "%d-%m-%Y")
			}
			contribution$parentId <- propertyValue(project,"id")
			if(any(names(contribution) == "tissue")){
				names(contribution)[which(names(contribution) == "tissue")] <- 'tissueType'
			}
			if(any(names(contribution) == "platformName")){
				names(contribution)[which(names(contribution) == "platformName")] <- 'platform'
			}
			return(contribution)
		}
)

setMethod(
		f = "handleContribution",
		signature = signature("Folder", "list","character"),
		definition = function(project, contribution, logFile) {
			# First remove any annotations set to NA
			if(sum(is.na(contribution))){
				contribution <- contribution[-which(is.na(contribution))]
			}
			# Next, make sure the user has provided the required information.
			if(! all(c("data.url","study.name","data.name", "data.type","data.status") %in% names(contribution))){				
				cat(setdiff(c("data.url","study.name","data.name", "data.type","data.status"), names(contribution)),"\n")
				stop("Mising one of study.name, data.url, data.name, data.type, or data.status.  Please try again.\n")
			}
			# Set the data.lastUpdate to NA if it wasn't provided
			if(! "data.lastUpdate" %in% names(contribution)){
				# Last update not included.  Set to NA
				contribution$data.lastUpdate = format(Sys.time(), "%d-%m-%Y")
			}
			contribution$parentId <- propertyValue(project,"id")
			if(any(names(contribution) == "tissue")){
				names(contribution)[which(names(contribution) == "tissue")] <- 'tissueType'
			}
			if(any(names(contribution) == "platformName")){
				names(contribution)[which(names(contribution) == "platformName")] <- 'platform'
			}
			return(contribution)
		}
)
setMethod(
		f = "handleContribution",
		signature = signature("NULL", "data.frame", "character"),
		definition = function(project, contribution, logFile){
			# here we assume user has supplied a dataframe, with
			# each row containing a single contribution
			if(! all(c("data.url","study.name","data.name", "data.type","data.status", "parentId") %in% names(contribution))){
				stop("Mising one of study.name, data.url, data.name, data.type, or data.status.  Please try again.\n")
			}
			proj <- contribution$parentId
			projects <- sapply(unique(proj), handleProject)
			names(projects) <- sapply(projects, function(project){ 
						propertyValue(project,"id")		
					})
			contribution <- lapply(1:nrow(contribution),function(x){ as.list(contribution[x,])})
			sapply(contribution, function(x){
						project <- projects[[as.character(x$parentId)]]
						handleContribution(project, x, logFile)
					}) -> contribution
			tmp <- list()
			if(class(contribution) == "matrix"){				
				tmp <- apply(contribution,2,as.list)
				contribution <- tmp
			}
			contribution
		}
)


setMethod(
		f = "handleContribution",
		signature = signature("NULL", "character", "character"),
		definition = function(project, contribution, logFile){
			# here we assume user has supplied a file.  we assume
			# each row contains a single contribution
			if(!file.exists(contribution)) {
				stop("Cannot find file", contribution,"\n.  Please try again.")
			}
			contribution <- read.delim(contribution, 
					sep="\t",
					quote="",
					strip.white=TRUE,
					stringsAsFactors=FALSE)
			contribution <- handleContribution(project, contribution, logFile)
			return(contribution)
		}
)

