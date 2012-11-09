#
# Function to handle location provided by contributor.  Will download and 
# check md5 sum if required.
#

setMethod(
		f = "handleLocation",
		signature = c("list","logical"),
		definition = function(contribution, keepLocal){
			# Set all characters as lowercase for each name.  Check that url is provided
			# If not, then return an error.
			if(! any(names(contribution) %in% c("url"))){
				stop("Error in handleLocation. Contribution must have an element named data.url.\n")
			}
			url <- contribution[[which(names(contribution) == "url")]]
			
			# Use the URL and synapseCache variable to build the file location.  If it does not
			# exist and the user has not provided the md5, then download the file and calculate the md5.
			tmp <- .parseURL(url)
			destfile <- tmp[1]
			destdir <- tmp[2]
			if(! (file.exists(destfile) | 'md5' %in% names(contribution)) ){
				# If the file doesn't already exist and md5 sum not provided by user
				# then download it.
				fileSize <- .getFileSize(url)
				if(fileSize > (20 * 1073741824)){
					# File is more than 20GB large.
					var <- "File too large.  Will not download files larger than 20GB.\n"
					class(var) <- "try-error"
					return(var)
				}
				if(!file.exists(destdir)){
					dir.create(destdir, recursive = TRUE)
				}
				destfile2 =  tempfile()
				for(i in 1:5) {  
					# Try 5 times in case our download breaks part of the way 
					# through.  If it does break, remove whatever file was there.
					download <- try(synapseClient:::.curlWriterDownload(url,destfile))
					if(class(download) != "try-error"){
						break;
					}else{
						file.remove(destfile)
					}
				}
				if(i == 5){
					var <- "Could not download data\n"
					class(var) <- "try-error"
					return(var)
				}
			}
			# Take md5 sum and build location object.
			if(file.exists(destfile)){
				checksum <- as.character(tools::md5sum(destfile))
			}else if('md5' %in% names(contribution)){
				checksum <- contribution$md5
			}
			dataLocation = list(path=url,type="external")
			# Add objects to contribution and return
			contribution$md5 = checksum
			contribution$locations = list(dataLocation)
#			contribution <- contribution[[-which(names(contribution)=="dataLocation")]]
#			contribution <- contribution[[-which(names(contribution)=="checksum")]]
			return(list(contribution=contribution, destfile=destfile))
		}
)

.parseURL <- function(url){ 
	parsedUrl <- synapseClient:::.ParsedUrl(url)
	destfile <- file.path(synapseCacheDir(), gsub("^/", "", parsedUrl@path))
	destfile <- path.expand(destfile)
	destdir <- gsub(parsedUrl@file, "", destfile, fixed = TRUE)
	destdir <- gsub("[\\/]+$", "", destdir)			
	c(destfile, destdir)
}


.getFileSize <- function(url){ 
	curlHandle = getCurlHandle() 
	h = basicTextGatherer()
	response <- curlPerform(URL=url, 
			.opts=list(header=TRUE, nobody=TRUE),  
			writefunction=h$update, 
			curl=curlHandle) 
	curlInfo <- getCurlInfo(curlHandle)
	return(curlInfo$content.length.download)
}

