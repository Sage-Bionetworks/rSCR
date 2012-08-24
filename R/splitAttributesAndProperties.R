
splitAttributesAndProperties <-function(a) {
	dataSetPropertyLabels<-c("name", "description", "platform", "createdBy", "parentId", "disease", "tissueType", "numSamples", "species", "md5","locations") 
	if("Description" %in% names(a)){
		names(a)[which(names(a) == "Description")] <- "description"
	}
	if( !("name" %in% names(a) & "parentId" %in% names(a))) {
		val <- "Cannot create dataset.  Both name and parentId are required inputs.";
		class(val) <- "try-error"
		return(val)
	}
	properties<-list()
	annotations<-list()
	for (i in 1:length(a)) {
		fieldName<-names(a[i])
		if(is.na(a[i]) | a[i] == "NA"){
			a[i] <- "Not Available"
		}
		if (any(dataSetPropertyLabels==fieldName)) {
			properties[fieldName]<-a[i]
		} else {
			annotations[fieldName]<-a[i]
		}
	}
	return(list(properties=properties, annotations=annotations))
}

