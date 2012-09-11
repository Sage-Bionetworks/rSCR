#
#  These functions handle user input for a project.  If its a character or numeric
#  then it assumes its a synapse entity id, downlods it and returns it.
#  If its not a valid id then it throws an error.
#  If its a project, then it tries to build the project in synapse if need be.
#  Both return the project.  There's no logic to test here that doesn't require
#  downloading something from synapse, so I won't build a unit test.
#

setMethod(
		f = "handleProject",
		signature = signature("character"),
		definition = function(project) {
			project <- try(getEntity(project),silent=TRUE)
			if(class(project) != "Project" | class(project) != "Folder"){
				# We could create the project for the user instead of stopping the session.  Requires more input
        # like a name, etc.
				stop("User provided project id ", project, " does not exist.\nPlease provide a valid Synapse ID or object of class project.")
			}
			return(project)
		}
)

setMethod(
		f = "handleProject",
		signature = signature("numeric"),
		definition = function(project) {
			handleProject( as.character( project ))
		}
)



setMethod(
		f = "handleProject",
		signature = signature("Project"),
		definition = function(project) {
			userProjectID <- propertyValue(project,"id")
			qryString <- sprintf("select * from project ")
			qryResult <- synapseQuery(qryString)		
			qry <- which(qryResult$project.id == userProjectID)
			if(length(qry) > 0) { 
				# Project exists
				project <- updateEntity(project)
			}else{
				# Project does not exist, so create it.
				project <- createEntity(project)
			}
			return(project)
		}
)


setMethod(
		f = "handleProject",
		signature = signature("NULL"),
		definition = function(project) {
			# We assume the user has supplied a file.
			return(NULL)
		}
)
