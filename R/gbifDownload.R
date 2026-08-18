### 1.1 ==== Utility functions used in GBIF Data Download and Visualisation ====

### 1.1.1 ---- Sanity check GBIF Username ----
#' @title Find GBIF Username
#'
#' @description A function that looks in the common places for a GBIF username
#' if the username is not passed as an argument to the function:
#' \code{'GBIF_USER'} system variable or from a call to
#' \code{getOption('gbif_user')}.
#'
#' @param val A character scalar containing the GBIF username
#'
#' @return A user name for GBIF
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
check_gbifUser <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF username: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF username: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_USER", "")
    if(inVal == "") {
      inVal <- getOption("gbif_user", stop("you must supply a GBIF username"))
    }
  }
  inVal
}

### 1.1.2 ---- Sanity check GBIF Password ----
#' @title Find GBIF Password
#'
#' @description A function that looks in the common places for a GBIF password
#' if the password is not passed as an argument to the function:
#' \code{'GBIF_PWD'} system variable or from a call to
#' \code{getOption('gbif_pwd')}.
#'
#' @param val A character scalar containing the GBIF password
#'
#' @return A password for GBIF
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
check_gbifPwd <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF password: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF password: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_PWD", "")
    if(inVal == "") {
      inVal <- getOption("gbif_pwd", stop("you must supply a GBIF password"))
    }
  }
  inVal
}

### 1.1.3 ---- Sanity check GBIF Email ----
#' @title Find GBIF Email
#'
#' @description A function that looks in the common places for a GBIF email
#' address if the address is not passed as an argument to the function:
#' \code{'GBIF_EMAIL'} system variable or from a call to
#' \code{getOption('gbif_email')}.
#'
#' @param val A character scalar containing the GBIF username
#'
#' @return A user name for GBIF
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
check_gbifEmail <- function(val) {
  inVal <- tryCatch(as.character(val), error = function(err) {
    stop("invalid entry for GBIF email: ", err)
  })
  inVal <- inVal[!is.na(inVal)]
  if(length(inVal) > 1) {
    warning("more than one value given as a GBIF email: only the first entry will be used")
    inVal <- inVal[1]
  } else if(length(inVal) <= 0) {
    inVal <- Sys.getenv("GBIF_EMAIL", "")
    if(inVal == "") {
      inVal <- getOption("gbif_email", stop("you must supply a GBIF email"))
    }
  }
  inVal
}

### 1.1.4 ---- Function to Retrieve Default CRS for GBIF ----
#' @title Retrieve GBIF's Default CRS
#'
#' @description A function that returns the default coordinate reference system
#' information that is used by GBIF's API
#'
#' @return A \code{crs} object containing GBIF's coordinate reference system
#' information
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
gbifDefaultCRS <- function() {
  sf::st_crs("+proj=longlat +datum=WGS84 +no_defs +type=crs")
}

### 1.1.5 ---- Process Dataset Titles ----
#' @title Process Dataset Titles
#'
#' @description Utility function to take the titles of datsets and wrap them
#' over multiple lines if they are long
#'
#' @param inTitles A character vector containing the titles of datasets
#' @param targetWidth An integer scalar giving the target width (in characters)
#' to wrap the titles into
#'
#' @return A character vector containing the processed dataset title names
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @noRd
processDatasetTitle <- function(inTitles, targetWidth) {
  sapply(X = ifelse(is.na(inTitles) | gsub("\\s+", "", inTitles, perl = TRUE) == "", "Unknown", inTitles), FUN = function(curTitle, targetWidth) {
    paste(strwrap(curTitle, width = targetWidth), collapse = "\n")
  }, targetWidth = targetWidth)
}

### 1.2 ==== Function to Retrieve Taxon Information from GBIF ====
#' @title Lookup GBIF Taxonomic Information
#'
#' @description A function that iteratively calls the
#' \code{\link[rgbif]{name_lookup}} function to retrieve information relating to
#' a taxon (or selection of taxa) contained at [GBIF](https://www.gbif.org/).
#'
#' @param kingdom A character vector containing the kingdoms that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the kingdoms specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param phylum A character vector containing the phyla that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the phyla specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param class A character vector containing the classes that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the classes specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param order A character vector containing the orders that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the orders specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param family A character vector containing the families that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the families specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param genus A character vector containing the genera that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the genera specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param species A character vector containing the species that the taxa
#' retrieved from GBIF's taxonomies database must belong to. An entry with a
#' length greater than one will ensure that all retrieved entries must belong
#' to one of the species specified in the input. A value of \code{NULL} ensures
#' that entries will not be filtered according to any criteria at this taxonomic
#' level.
#' @param other A character vector containing the values for the taxonomic level
#' specified in the \code{rank} argument that the taxa retrieved from GBIF's
#' taxonomies database must belong to. An entry with a length greater than one
#' will ensure that all retrieved entries must belong to one of the values
#' specified in the input. A value of \code{NULL} ensures that entries will not
#' be filtered according to any criteria at this taxonomic level.
#' @param limit An integer scalar containing the maximum number of entries to
#' return from the query. A value of \code{NULL} indicates that all valid
#' entries will be returned.
#' @param rank If the \code{other} parameter is non-\code{NULL} then a scalar
#' character must be provided here defining the taxonomic level to conduct the
#' query at. See \code{\link[rgbif]{name_lookup}} for acceptable values to use
#' for the \code{rank} argument.
#' @param ... Other parameters passed to the \code{\link[rgbif]{name_lookup}}
#' function when querying GBIF's taxonomies database.
#'
#' @return A \code{data.frame} containing the taxonomic information for all
#' matching taxa according to the criteria specified in the function arguments.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[rgbif]{name_lookup}}
#' @export
taxonLookup <- function(
    kingdom = NULL,
    phylum = NULL,
    class = NULL,
    order = NULL,
    family = NULL,
    genus = NULL,
    species = NULL,
    other = NULL,
    limit = NULL,
    rank = NULL,
    ...
) {
  # Set the hard limit set by the GBIF API
  gbifLimit <- 99999
  # Sanity check the limit
  inLimit <- tryCatch(as.integer(limit), error = function(err) {
    stop("error encountered processing the limit argument: ", err)
  })
  if(length(inLimit) <= 0) {
    inLimit <- Inf
  } else if(length(inLimit) > 1) {
    warning("limit argument has length greater than one: only the first element will be used")
    inLimit <- inLimit
  }
  if(is.na(inLimit)) {
    inLimit <- Inf
  }
  # Function to process the taxa specification
  taxaTextProcess <- function(inVec) {
    outTaxa <- tryCatch(as.character(inVec), error = function(err) {
      stop("error encountered processing taxonomic information: ", err)
    })
    outTaxa[!is.na(outTaxa) & outTaxa != ""]
  }
  # Populate a taxonomic search list based on the user entries
  taxaSpec <- list(
    other = taxaTextProcess(other),
    species = taxaTextProcess(species),
    genus = taxaTextProcess(genus),
    family = taxaTextProcess(family),
    order = taxaTextProcess(order),
    class = taxaTextProcess(class),
    phylum = taxaTextProcess(phylum),
    kingdom = taxaTextProcess(kingdom)
  )
  # Set the rank entry (only used if 'other' parameter is specified)
  rankInfo <- taxaTextProcess(rank)
  if(length(rankInfo) > 1) {
    warning("rank argument has length greater than one: only the first element will be used")
    rankInfo <- rankInfo[1]
  }
  # Set the taxonomic level of the query
  qryLvlIndx <- which(sapply(X = taxaSpec, FUN = length) > 0)
  if(length(qryLvlIndx) <= 0) {
    stop("no taxonomic information has been provided to query")
  }
  qryLvlIndx <- qryLvlIndx[1]
  qryLvl <- names(taxaSpec)[qryLvlIndx[1]]
  # Retrieve the list of ... parameters
  extraArgs <- eval(substitute(list(...)))
  if(!is.null(extraArgs$verbose)) {
    warning("entry for 'verbose' argument provided but this value will be over-ridden")
  }
  extraArgs$verbose <- FALSE
  startIndex <- 0
  # Prcoess the 'start' argument if it is provided in
  if(!is.null(extraArgs$start)) {
    startIndex <- tryCatch(as.integer(extraArgs$start), error = function(err) {
      stop("invalid entry for the start argument: ", err)
    })
    if(length(startIndex) <= 0) {
      startIndex <- 0
    } else if(length(startIndex) > 1) {
      warning("start argument has length greater than one: only the first element will be used")
      startIndex <- startIndex[1]
    }
    if(is.na(startIndex) || startIndex < 0) {
      stop("invalid entry for the start argument: value is NA or less than zero")
    }
  }
  outList <- lapply(X = taxaSpec[[qryLvl]], FUN = function(curQueryVal, taxaSpec, qryLvl, extraArgs, rankInfo, inLimit, gbifLimit, startIndex) {
    # Retrieve the rank of the current query
    curRank <- qryLvl
    if(curRank == "other") {
      if(length(rankInfo) <= 0) {
        stop("rank information not set when the 'other' taxonomic level is specified")
      }
      curRank <- rankInfo
    }
    # Break up the query into smaller chunks (with maximum size set by the GBIF limit)
    # and then stitch them all together
    outFrame <- NULL
    curStart <- startIndex
    endOfRecords <- FALSE
    message("retrieving entries from GBIF names database that have \"", curRank, "\" values similar to ", curQueryVal, "...")
    while(!endOfRecords && curStart < inLimit) {
      queryOut <- do.call(rgbif::name_lookup, c(list(
        query = curQueryVal, rank = curRank,
        limit = min(inLimit - curStart, gbifLimit),
        start = curStart), extraArgs))
      curStart <- curStart + nrow(queryOut$data)
      endOfRecords <- all(queryOut$meta[, "endOfRecords"], na.rm = TRUE)
      outFrame <- rbind(outFrame, as.data.frame(queryOut$data))
    }
    # Filter the data frame based on the other rank restrictions
    curQryInd <- which(names(taxaSpec) == qryLvl)[1] + 1
    rowsToUse <- rep(TRUE, nrow(outFrame))
    if(curQryInd < length(taxaSpec)) {
      rowsToUse <- apply(X = sapply(X = curQryInd:length(taxaSpec), FUN = function(curQryIndVal, taxaSpec, outFrame) {
        curAllowableValues <- taxaSpec[[curQryIndVal]]
        curCol <- names(taxaSpec)[curQryIndVal]
        rowsToUse <- rep(TRUE, nrow(outFrame))
        if(length(curAllowableValues) > 0) {
          rowsToUse <- tolower(outFrame[, curCol]) %in% tolower(curAllowableValues)
        }
        rowsToUse
      }, taxaSpec = taxaSpec, outFrame = outFrame), MARGIN = 1, FUN = all)
    }
    outFrame[rowsToUse, ]
  }, taxaSpec = taxaSpec, qryLvl = qryLvl, extraArgs = extraArgs, rankInfo = rankInfo, inLimit = inLimit, gbifLimit = gbifLimit, startIndex = startIndex)
  # Helper function to bind table when columns differ (occasionally happens due to weird ways the taxon lookup happens)
  rbindPad <- function(inList) {
    allColNames <- unique(unlist(lapply(X = inList, FUN = colnames)))
    outList <- lapply(X = inList, FUN = function(curFrame, allColNames) {
      outFrame <- as.data.frame(curFrame)
      isInFrame <- allColNames %in% colnames(outFrame)
      if(!all(isInFrame)) {
        missingCols <- allColNames[!isInFrame]
        outFrame <- cbind(outFrame, as.data.frame(stats::setNames(replicate(length(missingCols), rep(NA, nrow(outFrame)), simplify = FALSE), missingCols)))
      }
      outFrame[, allColNames]
    }, allColNames = allColNames)
    do.call(rbind, outList)
  }
  rbindPad(outList)
}

### 1.3 ==== Function to Produce a JSON Query for GBIF's API ====
#' @title Create a JSON Query for GBIF's Occurrence API
#'
#' @description A function to create a JSON query for the
#' [occurrence download API](https://techdocs.gbif.org/en/data-use/api-downloads.html)
#' of [GBIF](https://www.gbif.org/).
#'
#' @param ... A set of queries to combine to create for GBIF's occurrence
#' download API (see details below).
#' @param user A character string containing a valid username associated with a
#' GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_USER')} and
#' \code{getOptions('gbif_user')} are searched for a value instead.
#' @param email A character string containing a valid email associated with a
#' GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_EMAIL')} and
#' \code{getOptions('gbif_email')} are searched for a value instead.
#' @param sendNotification A scalar logical that, if \code{TRUE}, adds a portion
#' to the query that requests that an email be sent to the account specified in
#' \code{email} when the download is ready.
#' @param geodistance.values Argument is only used when using a 'geodistance'
#' query (see details below).
#'
#' @details Queries to be converted to JSON are specified using named
#' parameters. Any named parameter that shares a name with one of the parameters
#' listed in GBIF's
#' [occurrence search API](https://techdocs.gbif.org/en/openapi/v1/occurrence#/Searching%20occurrences/searchOccurrence)
#' will be interpreted as a filter on the relevant records. For example, having
#' an argument entitled \code{COUNTRY} with value \code{"NO"} will result a
#' query that will only download records that have an ISO-3166-1 country code of
#' NO (that is all record registered in Norway). If the vector provided for this
#' argument has a length greater than one, then all records that match any of
#' the composite values will be returned. These query types correspond to the
#' 'equals' and 'in' predicates described in the
#' [occurrence download API document](https://techdocs.gbif.org/en/data-use/api-downloads.html)
#' respectively. Other predicates can be specified by adding a suffix to the
#' named argument. For example:
#' \describe{
#'  \item{\code{key.min}}{Will set a query that will only return records that
#'  have a value for the \code{key} that are greater than or equal to the value
#'  provided in the argument. For example, supplying \code{YEAR.min = 1960} as
#'  an argument will result in a query that only returns records that have been
#'  collected from the year 1960 onward. There is also the suffix \code{.xmin}
#'  that can be used to create a query that will retrieve records only if they
#'  have a relevant value for \code{key} that are greater than the value
#'  provided in the argument.}
#'  \item{\code{key.max}}{Will set a query that will only return records that
#'  have a value for the \code{key} that are less than or equal to the value
#'  provided in the argument. For example, supplying \code{YEAR.max = 1990} as
#'  an argument will result in a query that only returns records that have been
#'  collected in the year 1990 or earlier onward. There is also the suffix
#'  \code{.xmax} that can be used to create a query that will retrieve records
#'  only if they have a relevant value for \code{key} that are less than the
#'  value provided in the argument.}
#'  \item{\code{key.range}}{Will set a query that will only return records that
#'  have a value for the \code{key} that are greater than or equal to the
#'  minimum value and less than or equal to the maximum value provided in the
#'  vector of values provided for this argument. For example, supplying
#'  \code{ELEVATION.range = c(400, 600)} as an argument will result in a query
#'  that only returns records that have been collected at an elevation of
#'  between 400 and 600 metres above sea level.}
#'  \item{\code{key.like}}{Will set a query that will only return records that
#'  have a value for the \code{key} that match a pattern set in the value given
#'  for the argument. Details on how the string provided for the pattern
#'  matching is interpreted in the specification for the 'like' predicate in the
#'  [occurrence download API document](https://techdocs.gbif.org/en/data-use/api-downloads.html).}
#'  \item{\code{key.notnull}}{If the value of this argument is \code{TRUE} then
#'  it will result in a query that will only include records that have a
#'  non-empty value for \code{key}. Similarly a value of \code{FALSE} will
#'  result in a query that will only include records that have an empty value
#'  for \code{key}.}
#'  \item{\code{key.not}}{This presence of this suffix represents the logical
#'  negation for all components of the suffix that comes before it. For example,
#'  supplying \code{DATASET_KEY.not = "4fa7b334-ce0d-4e88-aaae-2e0c138d049e"}
#'  will result in a query that includes records from all datasets except that
#'  with the provided dataset key. Similarly, supplying
#'  \code{ELEVATION.range.not = c(400, 600)} as an argument will result in a
#'  query that retrieves all records that do not fall within the elevational
#'  band that lies between 400 and 600 metres above sea level.}
#' }
#' This function will accept arguments that are named either using the uppercase
#' notation used in the examples provided here or using the camel case notation
#' used in the \code{\link[rgbif]{pred}} function.
#'
#' By default the arguments \code{HAS_GEOSPATIAL_ISSUE = FALSE},
#' \code{HAS_COORDINATE = TRUE}, \code{OCCURRENCE_STATUS = "PRESENT"}, and
#' \code{BASIS_OF_RECORD.not = c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")} are
#' added to the list of arguments provided by the user. These are commonly
#' applied defaults and are added to avoid accidental requests for erroneous
#' records. If the user does not want these default to be applied then they can
#' override this by providing \code{NULL} as the relevant argument. For example,
#' supplying \code{OCCURRENCE_STATUS = NULL} will result in a query that
#' requests records that represent both presences and absences.
#'
#' If an argument with the name '\code{geometry}' appears in the list of query
#' arguments with a value represented by an object of class \code{sfc} (or that
#' inherits from it), then a query will be produced that only includes records
#' that fall within the polygon specified in the object. Similarly, if an
#' argument with the name '\code{geodistance}' appears in the list of query
#' arguments with a value with class that inherits from \code{sfc} than a query
#' will be produced that only included records that fall within a distance
#' (specified by the \code{geodistance.values} argument) of the coordinates (as
#' resulting from the \code{\link[sf]{st_coordinates}} function) that make up
#' the relevant spatial feature.
#'
#' If any arguments are passed to the function with names that match those in
#' \code{\link{taxonLookup}} and/or \code{\link[rgbif]{name_lookup}} functions,
#' then a call will be made to \code{taxonLookup} using those arguments. The
#' resulting data frame will be added as a \code{"taxonFrame"} attribute to
#' the output. The keys associated with the finest taxonomic resolution in the
#' taxonomic lookup query will then be used as keys to refine the occurrence
#' record search. For example, providing the arguments
#' \code{species = "Aeshna grandis"} and \code{order = "Odonata"}, will result
#' in a data frame returned in the \code{"taxonFrame"} attribute of the output
#' that will include the column \code{"speciesKey"} representing the keys in
#' GBIF's taxonomic database that are possible entries for the species
#' Aeshna grandis in the order Odonata. These keys values will be used in a
#' query as though they were passed to the function using the argument notation
#' \code{SPECIES_KEY = vectorOfSpeciesKeys} and therefore ensuring that only
#' records corresponding to this species are requested.
#'
#' @return A \code{jsonlite} object containing a query to use in GBIF's
#' occurrence download API. This object may have the \code{"taxonFrame"}
#' attribute defined if a call to \code{\link{taxonLookup}} was made. If so,
#' this attribute will be set to the data frame resulting from the taxonomic
#' query.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[rgbif]{name_lookup}}, \code{\link{taxonLookup}},
#'  \code{\link{gbifOccDownload}}, \code{\link[sf]{st_sfc}},
#'  \code{\link[sf]{st_coordinates}}
#' @export
jsonOccFormulation <- function(
    ...,
    user = NULL, email = NULL,
    sendNotification = TRUE,
    geodistance.values = NULL
) {
  # Function to check units support
  unitSupport <- function(inVec) {
    outVal <- tryCatch(units(inVec), error = function(err) {
      NULL
    })
    if(is.null(outVal)) {
      outVal <- ""
    }
    outVal
  }
  # Sanity check the arguments
  inSendNotification <- tryCatch(as.logical(sendNotification), error = function(err) {
    stop("error encountered processing the \"sendNotification\" argument: ", err)
  })
  if(length(inSendNotification) <= 0) {
    inSendNotification <- FALSE
  } else if(length(inSendNotification) > 1) {
    warning("\"sendNotification\" argument has length greater than one: only the first element will be used")
    inSendNotification <- inSendNotification[1]
  }
  if(is.na(inSendNotification)) {
    inSendNotification <- FALSE
  }
  geodistText <- NULL
  if(!is.null(geodistance.values)) {
    if(unitSupport(geodistance.values) == "") {
      geodistText <- tryCatch(paste0(as.character(geodistance.values), "m"), error = function(err) {
        stop("unable to process the geodistance values: ", err)
      })
    } else {
      geodistText <- tryCatch(paste0(as.character(geodistance.values), units(geodistance.values)), error = function(err) {
        stop("unable to process the geodistance values: ", err)
      })
    }
  }
  if(length(geodistText) <= 0) {
    geodistText <- NULL
  }
  # Retrieve the ... arguments
  predArgs <- eval(substitute(list(...)))
  # If any arguments are the same as the non ... arguments used in the taxon lookup
  # function then call that function first to get a list of the taxon names
  taxonLookupArgs <- methods::formalArgs(taxonLookup)
  taxonLookupArgs <- taxonLookupArgs[taxonLookupArgs != "..."]
  nameLookupArgs <- methods::formalArgs(rgbif::name_lookup)
  nameLookupArgs <- nameLookupArgs[!(nameLookupArgs %in% taxonLookupArgs)]
  taxonFrame <- NULL
  if(any(names(predArgs) %in% taxonLookupArgs)) {
    # Check to see the level that the taxon look up is being made
    lookupLevel <- NA
    taxaLevels <- c("species", "genus", "family", "order", "class", "phylum", "kingdom")
    taxaIsInPreds <- sapply(X = taxaLevels, FUN = function(curLevel, predArgs) { !is.null(predArgs[[curLevel]]) }, predArgs = predArgs)
    if(any(taxaIsInPreds)) {
      lookupLevel <- taxaLevels[which(taxaIsInPreds)[1]]
    }
    taxonFrame <- do.call(taxonLookup, predArgs[c(taxonLookupArgs, nameLookupArgs)[c(taxonLookupArgs, nameLookupArgs) %in% names(predArgs)]])
    # Remove the taxon-lookup related arguments from the processing list and also
    # add a search for the taxa found in the lookup table
    taxonLookupIDs <- unique(taxonFrame[, ifelse(is.na(lookupLevel), "key", paste0(lookupLevel, "Key"))])
    predArgs <- append(
      predArgs[!(names(predArgs) %in% c(taxonLookupArgs, nameLookupArgs))],
      stats::setNames(list(taxonLookupIDs[!is.na(taxonLookupIDs)]), paste(ifelse(is.na(lookupLevel), "TAXON", toupper(lookupLevel)), "KEY", sep = "_"))
    )
  } else if(any(names(predArgs) %in% nameLookupArgs)) {
    # Remove those arguments that are only arguments to the names lookup function of rgbif
    predArgs <- predArgs[!(names(predArgs) %in% nameLookupArgs)]
  }
  # Go through and create the list of predicates
  allPreds <- lapply(X = names(predArgs), FUN = function(curArgName, predArgs, geodistText) {
    curArgs <- predArgs[[curArgName]]
    if(is.logical(curArgs)) {
      # Logical vectors are converted into lowercase true/false values
      curArgs <- tolower(as.character(curArgs))
    } else if(inherits(curArgs, "sfg")) {
      # Input is a simple feature geometry so take its WKT string form and
      # use it directly
      curArgs <- sf::st_as_text(curArgs)
    } else if(!is.null(curArgs) && !inherits(curArgs, "sf") && !inherits(curArgs, "sfc")) {
      # Otherwise, if it is not an sf or sfc object, convert it into a character
      curArgs <- paste0(as.character(curArgs), unitSupport(curArgs))
    }
    # Function to process the predicate
    processPredicate <- function(curArgName, curArgs, geodistText) {
      outPred <- character()
      tempPred <- NULL
      attr(outPred, "key") <- NA
      if(is.null(curArgs)) {
        attr(outPred, "key") <- rgbif::pred_notnull(curArgName)$parameter
      } else if(tolower(curArgName) == "geometry") {
        # Current argument is a special geometry predicate
        if(!is.character(curArgs)) {
          # Argument has an sf or sfc class
          if(inherits(curArgs, "sf")) {
            # Just retrieve the geometry object from the sf class
            curArgs <- sf::st_geometry(curArgs)
          }
          if(inherits(curArgs, "sfc_MULTIPOLYGON") || inherits(curArgs, "sfc_POLYGON")) {
            # Firstly ensure that the spatial object has the same coordinate
            # system then the one the GBIF API uses
            curArgs <- sf::st_transform(curArgs, gbifDefaultCRS())
            # Now union all the elements in the feature and then retrieve the
            # WKT representation of the object
            curArgs <- sf::st_as_text(sf::st_union(curArgs)[[1]])
            curArgs <- gsub("\\s+", " ", gsub("^MULTIPOLYGON\\s+\\(", "MULTIPOLYGON(", gsub("^POLYGON\\s+\\(", "POLYGON(", toupper(curArgs), perl = TRUE), perl = TRUE), perl = TRUE)
          } else {
            stop("geometry key must have a value that is a POLYGON or MULTIPOLYGON simple feature")
          }
        }
        outPred <- c(
          "{",
          "\t\"type\": \"within\",",
          paste0("\t\"geometry\": \"", curArgs, "\""),
          "}")
        attr(outPred, "key") <- "geometry"
      } else if(tolower(curArgName) == "geodistance") {
        # Current argument is a special geospatial distance predicate
        if(!is.null(geodistText)) {
          stop("the \"geodistance.values\" argument must be non-NULL when the geodistance predictate is used")
        }
        if(is.character(curArgs)) {
          # If the input is a WKT text string then use that directly
          curArgs <- sf::st_as_sf(data.frame(wkttext = curArgs), wkt = "wkttext", crs = gbifDefaultCRS())
        }
        # Retrieve the coordinates associated with the
        curArgs <- sf::st_coordinates(curArgs)
        if(nrow(curArgs) == 1) {
          outPred <- c(
            "{",
            "\t\"type\": \"geoDistance\",",
            paste0("\t\"latitude\": \"", curArgs$Y, "\","),
            paste0("\t\"longitude\": \"", curArgs$X, "\","),
            paste0("\t\"distance\": \"", geodistText[1], "\","),
            "}")
        } else if(nrow(curArgs) > 1) {
          outPred <- c(
            "{",
            "\t\"type\": \"or\",",
            "\t\"predicates\": [",
            sapply(X = 1:nrow(curArgs), FUN = function(curIndex, xVals, yVals, distVals) {
              paste0("\t\t{\"type\": \"geoDistance\", \"latitude\": \"", yVals[curIndex], "\", \"longitude\": \"", xVals[curIndex], "\", \"distance\": \"", distVals[curIndex], "\"}")
            }, xVals = curArgs$X, yVals = curArgs$Y, distVals = geodistText[(1:nrow(curArgs) - 1) %% nrow(curArgs) + 1]),
            "\t]",
            "}")
        }
        attr(outPred, "key") <- "GEODISTANCE"
      } else {
        if(grepl("\\.min$", curArgName, perl = TRUE)) {
          # Argument is a "greaterThanOrEquals" predicate
          curArgName <- gsub("\\.min", "", curArgName, perl = TRUE)
          lowVal <- min(as.numeric(curArgs), na.rm = TRUE)
          tempPred <- rgbif::pred_gte(curArgName, lowVal)
          outPred <- c(
            "{",
            paste0("\t\"type\": \"", tempPred$type, "\","),
            paste0("\t\"key\": \"", tempPred$key, "\","),
            paste0("\t\"value\": \"", as.character(lowVal), "\""),
            "}")
        } else if(grepl("\\.minx$", curArgName, perl = TRUE)) {
          # Argument is a "greaterThan" predicate
          curArgName <- gsub("\\.minx", "", curArgName, perl = TRUE)
          lowVal <- min(as.numeric(curArgs), na.rm = TRUE)
          tempPred <- rgbif::pred_gt(curArgName, lowVal)
          outPred <- c(
            "{",
            paste0("\t\"type\": \"", tempPred$type, "\","),
            paste0("\t\"key\": \"", tempPred$key, "\","),
            paste0("\t\"value\": \"", as.character(lowVal), "\""),
            "}")
        } else if(grepl("\\.max$", curArgName, perl = TRUE)) {
          # Argument is a "lessThanOrEquals" predicate
          curArgName <- gsub("\\.max$", "", curArgName, perl = TRUE)
          highVal <- max(as.numeric(curArgs), na.rm = TRUE)
          tempPred <- rgbif::pred_lte(curArgName, highVal)
          outPred <- c(
            "{",
            paste0("\t\"type\": \"", tempPred$type, "\","),
            paste0("\t\"key\": \"", tempPred$key, "\","),
            paste0("\t\"value\": \"", as.character(highVal), "\""),
            "}")
        } else if(grepl("\\.maxx$", curArgName, perl = TRUE)) {
          # Argument is a "lessThan" predicate
          curArgName <- gsub("\\.maxx$", "", curArgName, perl = TRUE)
          highVal <- max(as.numeric(curArgs), na.rm = TRUE)
          tempPred <- rgbif::pred_lt(curArgName, highVal)
          outPred <- c(
            "{",
            paste0("\t\"type\": \"", tempPred$type, "\","),
            paste0("\t\"key\": \"", tempPred$key, "\","),
            paste0("\t\"value\": \"", as.character(highVal), "\""),
            "}")
        } else if(grepl("\\.range$", curArgName, perl = TRUE)) {
          # Argument is a range of values determined by two predicates
          curArgName <- gsub("\\.range$", "", curArgName, perl = TRUE)
          lowVal <- min(as.numeric(curArgs), na.rm = TRUE)
          highVal <- max(as.numeric(curArgs), na.rm = TRUE)
          tempPred <- list(
            rgbif::pred_gte(curArgName, lowVal),
            rgbif::pred_lte(curArgName, highVal))
          outPred <- c(
            "{",
            "\t\"type\": \"and\",",
            "\t\"predicates\": [",
            paste0("\t\t{\"type\": \"", tempPred[[1]]$type, "\", \"key\": \"", tempPred[[1]]$key, "\", \"value\": \"", as.character(lowVal), "\"},"),
            paste0("\t\t{\"type\": \"", tempPred[[2]]$type, "\", \"key\": \"", tempPred[[2]]$key, "\", \"value\": \"", as.character(highVal), "\"}"),
            "\t]",
            "}")
        } else if(grepl("\\.like$", curArgName, perl = TRUE)) {
          # Argument is a "like" predicate
          curArgName <- gsub("\\.like$", "", curArgName, perl = TRUE)
          if(length(curArgs) == 1) {
            # Argument is of length 1 so only one search done
            tempPred <- rgbif::pred_like(curArgName, curArgs)
            outPred <- c(
              "{",
              paste0("\t\"type\": \"", tempPred$type, "\", "),
              paste0("\t\"key\": \"", tempPred$key, "\", "),
              paste0("\t\"value\": \"", as.character(curArgs), "\""),
              "}")
          } else if(length(curArgs) > 1) {
            # Argument is of length greater than one so wrap the multiple 'like'
            # predicates in a big 'or' predicate
            tempPred <- lapply(X = as.character(curArgs), FUN = function(curArg, curArgName) {
              rgbif::pred_like(curArgName, curArg)
            }, curArgName = curArgName)
            outPred <- c(
              "{",
              "\t\"type\": \"or\",",
              "\t\"predicates\": [",
              paste0(sapply(X = tempPred, FUN = function(curPred) {
                paste0("\t\t{\"type\": \"", curPred$type, "\", \"key\": \"", curPred$key, "\", \"value\": \"", curPred$value, "\"}")
              }), c(rep(",", length(tempPred) - 1), "")),
              "\t]",
              "}")
          }
        } else if(grepl("\\.notnull$", curArgName, perl = TRUE)) {
          # Argument is a "notnull" predicate
          curArgName <- gsub("\\.notnull$", "", curArgName, perl = TRUE)
          tempPred <- rgbif::pred_notnull(curArgName)
          curArgs <- as.logical(curArgs)
          if(length(curArgs) > 0) {
            if(length(curArgs) > 1) {
              warning("vector of length greater than one provided as an argument for a \"notnull\" predicate: only the first element will be used")
              curArgs <- curArgs[1]
            }
          }
          if(curArgs) {
            # If arg.notnull = TRUE then only return those entries with non-null
            # entries for arg
            outPred <- c(
              "{",
              paste0("\t\"type\": \"", tempPred$type, "\","),
              paste0("\t\"parameter\": \"", tempPred$parameter, "\""),
              "}")
          } else {
            # If arg.notnull = FALSE then only return those entries with null
            # entries for arg
            outPred <- c(
              "{",
              "\t\"type\": \"isNull\",",
              paste0("\t\"parameter\": \"", tempPred$parameter, "\""),
              "}")
          }
        } else {
          if(length(curArgs) == 1) {
            # Argument has length of one so predicate is "equals"
            tempPred <- rgbif::pred(curArgName, curArgs)
            outPred <- c(
              "{",
              paste0("\t\"type\": \"", tempPred$type, "\","),
              paste0("\t\"key\": \"", tempPred$key, "\","),
              paste0("\t\"value\": \"", as.character(curArgs), "\""),
              "}")
          } else if(length(curArgs) > 1) {
            # Argument has length of one so predicate is "in"
            tempPred <- rgbif::pred_in(curArgName, curArgs)
            outPred <- c(
              "{",
              paste0("\t\"type\": \"", tempPred$type, "\","),
              paste0("\t\"key\": \"", tempPred$key, "\","),
              paste0("\t\"values\": [", paste0("\"", as.character(curArgs), "\"", collapse = ", "), "]"),
              "}")
          }
        }
        # Set the key attribute
        if(is.null(tempPred$key)) {
          if(is.null(tempPred$parameter)) {
            attr(outPred, "key") <- tempPred[[1]]$key
          } else {
            attr(outPred, "key") <- tempPred$parameter
          }
        } else {
          attr(outPred, "key") <- tempPred$key
        }
      }
      outPred
    }
    outPred <- character()
    if(grepl("\\.not$", curArgName, perl = TRUE)) {
      procPred <- processPredicate(gsub("\\.not$", "", curArgName, perl = TRUE), curArgs, geodistText)
      if(length(procPred) > 0) {
        procPred[1] <- paste("\"predicate\":", procPred[1])
        # Argument is a negation predicate so set the appropriate conditions
        outPred <- c(
          "{",
          "\t\"type\": \"not\",",
          paste0("\t", procPred),
          "}"
        )
      }
    } else {
      # Otherwise do the standard processing of the predicate
      outPred <- processPredicate(curArgName, curArgs, geodistText)
    }
    outPred
  }, predArgs = predArgs, geodistText = geodistText)
  # If some particular predicates have not been set then add some sensible defaults
  allKeys <- sapply(X = allPreds, FUN = function(curPred) { attr(curPred, "key") })
  allKeys <- allKeys[!is.na(allKeys)]
  if(!("HAS_GEOSPATIAL_ISSUE" %in% allKeys)) {
    allPreds <- append(allPreds, list(c(
      "{",
      "\t\"type\": \"equals\",",
      "\t\"key\": \"HAS_GEOSPATIAL_ISSUE\",",
      "\t\"value\": \"FALSE\"",
      "}"
    )))
  }
  if(!("HAS_COORDINATE" %in% allKeys)) {
    allPreds <- append(allPreds, list(c(
      "{",
      "\t\"type\": \"equals\",",
      "\t\"key\": \"HAS_COORDINATE\",",
      "\t\"value\": \"TRUE\"",
      "}"
    )))
  }
  if(!("OCCURRENCE_STATUS" %in% allKeys)) {
    allPreds <- append(allPreds, list(c(
      "{",
      "\t\"type\": \"equals\",",
      "\t\"key\": \"OCCURRENCE_STATUS\",",
      "\t\"value\": \"PRESENT\"",
      "}"
    )))
  }
  if(!("BASIS_OF_RECORD" %in% allKeys)) {
    allPreds <- append(allPreds, list(c(
      "{",
      "\t\"type\": \"not\",",
      "\t\"predicate\": {",
      "\t\t\"type\": \"in\",",
      "\t\t\"key\": \"BASIS_OF_RECORD\",",
      "\t\t\"values\": [\"FOSSIL_SPECIMEN\", \"LIVING_SPECIMEN\"]",
      "\t}",
      "}"
    )))
  }
  allPreds <- allPreds[sapply(X = allPreds, FUN = length) > 0]
  outQueryText <- NA
  predText <- ""
  if(length(allPreds) > 0) {
    if(length(allPreds) == 1) {
      # There is only one predicate so only use that
      allPreds[[1]][1] <- paste0("\"predicate\": ", allPreds[[1]][1])
      predText <- paste0("\t", allPreds[[1]], collapse = "\n")
    } else {
      # Otherwise wrap all the predicates in a big 'and' predicate
      predText <- paste0(c(
        "\t\"predicate\": {",
        "\t\t\"type\": \"and\",",
        "\t\t\"predicates\": [",
        paste0(sapply(X = allPreds, FUN = function(curPreds) {
          paste0("\t\t\t", curPreds, collapse = "\n")
        }), collapse = ",\n"),
        "\t\t]",
        "\t}"
      ), collapse = "\n")
    }
    outQueryText <- paste("{",
                          paste0("\t\"creator\":\"", check_gbifUser(user), "\","),
                          paste0("\t\"notification_address\":[\"", check_gbifEmail(email), "\"],"),
                          paste0("\t\"sendNotification\": ", tolower(as.character(inSendNotification)), ","),
                          predText,
                          "}\n", sep = "\n")
    outQueryText <- jsonlite::prettify(outQueryText, indent = 1)
  }
  # If a taxon lookup has been performed then add an attribute that returns the
  # outcome of that search
  if(!is.null(taxonFrame)) {
    attr(outQueryText, "taxonFrame") <- taxonFrame
  }
  outQueryText
}

### 1.4 ==== Download Occurrence Records from GBIF ====
#' @title Download Occurrence Records from GBIF
#'
#' @description Function to query the
#' [occurrence download API](https://techdocs.gbif.org/en/data-use/api-downloads.html)
#' of [GBIF](https://www.gbif.org/) and retrieve matching records.
#'
#' @param user A character string containing a valid username associated with a
#' GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_USER')} and
#' \code{getOptions('gbif_user')} are searched for a value instead.
#' @param pwd A character string containing a valid password associated with a
#' GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_PWD')} and
#' \code{getOptions('gbif_pwd')} are searched for a value instead.
#' @param email A character string containing a valid email associated with a
#' GBIF account. If \code{NULL} then \code{Sys.getenv('GBIF_EMAIL')} and
#' \code{getOptions('gbif_email')} are searched for a value instead.
#' @param curlopts List of named curl options passed on to
#' \code{\link[rgbif]{occ_download}}.
#' @param body A character scalar containing a syntactically valid JSON query
#' to pass to the GBIF occurrence download API. If this is \code{NULL} then the
#' arguments in \code{...} are instead passed to the
#' \code{\link{jsonOccFormulation}} to generate a query.
#' @param pingtime An integer scalar giving the number of seconds to wait
#' between pings to GBIF's servers to query the status of the data download.
#' @param tmploc A character scalar providing the location to store temporary
#' files created in the download (that are then deleted on completion).
#' @param assf A logical scalar that if \code{TRUE}, converts the output of the
#' occurrence record download to an \code{sf} object.
#' @param ... If \code{body} is \code{NULL} then these parameters are passed to
#' the \code{\link{jsonOccFormulation}} to generate a query.
#'
#' @return An object of type \code{gbifOcc} that inherits from an \code{sf}
#' object if \code{assf} is \code{TRUE} but otherwise directly inherits from
#' a \code{data.frame}. In addition the object will have the following
#' attributes:
#' \describe{
#'  \item{\code{DOI}}{A character scalar giving the DOI allocated to the data
#'  download.}
#'  \item{\code{citation}}{A character scalar giving the citation to use when
#'  referring to the dataset in published literature.}
#'  \item{\code{query}}{A \code{jsonlite} object containing the full query used
#'  to generate the download.}
#'  \item{\code{dataset}}{A data frame containing information on the constituent
#'  datasets used in the occurrence download.}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[rgbif]{name_lookup}}, \code{\link{taxonLookup}},
#'  \code{\link{jsonOccFormulation}}
#' @export
gbifOccDownload <- function(user, pwd, email, curlopts = list(), body = NULL, pingtime = 5, tmploc = tempdir(), assf = TRUE, ...) {
  # Sanity check the arguments
  inPingTime <- tryCatch(as.integer(pingtime), error = function(err) {
    stop("error encountered processing the \"pingtime\" argument: ", err)
  })
  if(length(inPingTime) <= 0) {
    inPingTime <- formals(gbifOccDownload)$pingtime
  } else if(length(inPingTime) > 1) {
    warning("length of \"pingtime\" argument greater than one: only the first element will be used")
    inPingTime <- inPingTime[1]
  }
  if(is.na(inPingTime)) {
    inPingTime <- formals(gbifOccDownload)$pingtime
  }
  inTmpDir <- tryCatch(as.character(tmploc), error = function(err) {
    stop("error encountered processing the \"tmploc\" argument: ", err)
  })
  if(length(inTmpDir) <= 0) {
    inTmpDir <- formals(gbifOccDownload)$tmploc
  } else if(length(inTmpDir) > 1) {
    warning("length of \"tmploc\" argument greater than one: only the first element will be used")
    inTmpDir <- inTmpDir[1]
  }
  if(is.na(inTmpDir)) {
    inTmpDir <- formals(gbifOccDownload)$tmploc
  }
  inQuery <- body
  if(is.null(inQuery)) {
    # If the body is not provided then generate it from the ... parameters
    inQuery <- do.call(jsonOccFormulation, append(list(user = user, email = email, curlopts = curlopts), eval(substitute(list(...)))))
  }
  isAsSF <- tryCatch(as.logical(assf), error = function(err) {
    stop("error encountered processing the \"assf\" argument: ", err)
  })
  if(length(isAsSF) <= 0) {
    isAsSF <- formals(gbifOccDownload)$assf
  } else if(length(isAsSF) > 1) {
    warning("length of \"assf\" argument greater than one: only the first element will be used"=)
    isAsSF <- isAsSF[1]
  }
  if(is.na(isAsSF)) {
    isAsSF <- formals(gbifOccDownload)$assf
  }
  # Generate a download ID
  message("sending download request to GBIF servers...")
  downloadID <- rgbif::occ_download(body = inQuery, user = check_gbifUser(user), pwd = check_gbifUser(pwd), email = check_gbifEmail(email), format = "DWCA", curlopts = curlopts)
  # Delete any existing copy of the data
  occDataDir <- file.path(inTmpDir, paste0("occData", as.character(downloadID)))
  on.exit({
    # Tidy up temporary directory (if it exists)
    if(dir.exists(occDataDir)) {
      unlink(occDataDir, recursive = TRUE)
    }
    # Cancel a download if it is still being prepared/run
    if(rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "RUNNING" ||
       rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "PREPARING") {
      rgbif::occ_download_cancel(downloadID, user = check_gbifUser(user), pwd = check_gbifUser(pwd), curlopts = curlopts)
    }
  })
  # Wait a minimum of 20 seconds before starting to query the download status.  For large download requests
  # immediate querying of the download can result in "RemoteException: File does not exist" errors from
  # GBIF's servers
  Sys.sleep(max(20, inPingTime))
  message("GBIF download being prepared...")
  while(rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "RUNNING" ||
        rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "PREPARING") {
    Sys.sleep(inPingTime)
  }
  if(rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "FAILED") {
    stop("error encountered processing the download request at the GBIF servers")
  } else if(rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "KILLED" ||
            rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "CANCELLED" ||
            rgbif::occ_download_meta(downloadID, curlopts = curlopts)$status == "FILE_ERASED") {
    stop("download request killed or cancelled (or file is no longer stored see https://www.gbif.org/faq?question=for-how-long-will-does-gbif-store-downloads)")
  }
  if(dir.exists(occDataDir)) {
    unlink(occDataDir, recursive = TRUE)
  }
  if(!dir.create(occDataDir)) {
    stop("unable to create temporary directory")
  }
  # Import the data and add the spatial information to it
  message("retrieving download from the GBIF servers...")
  occDataGet <- do.call(rgbif::occ_download_get, append(list(
    key = as.character(downloadID), path = occDataDir, overwrite = FALSE
  ), curlopts))
  # occDataGet <- rgbif::occ_download_get(as.character(downloadID), path = occDataDir)
  message("unpacking downloaded data locally...")
  occDataRaw <- as.data.frame(rgbif::occ_download_import(occDataGet, downloadID, path = occDataDir))
  if(isAsSF && nrow(occDataRaw) > 0) {
    occDataRaw <- sf::st_as_sf(occDataRaw, coords = c("decimalLongitude", "decimalLatitude"), crs = gbifDefaultCRS())
  }
  # Add citation attributes to the dataset
  attr(occDataRaw, "DOI") <- attr(downloadID, "doi")
  attr(occDataRaw, "citation") <- attr(downloadID, "citation")
  message("Download fully processed. Please cite this dataset download using the following style in any publications resulting from this analysis:\n",
          paste(attr(occDataRaw, "citation"), collapse = "\n"), ".\nSee the GBIF citation guidelines for more information: https://www.gbif.org/citation-guidelines")
  attr(occDataRaw, "query") <- inQuery
  # Function to retrieve extra information on the constituent datasets in the download
  retrieveDatasetInfo <- function(downloadID, curlopts) {
    endOfRecords <- FALSE
    curIndex <- 0
    gbifLimit <- 1000
    outFrame <- NULL
    while(!endOfRecords) {
      queryOut <- rgbif::occ_download_datasets(downloadID, limit = gbifLimit, start = curIndex, curlopts = curlopts)
      endOfRecords <- queryOut$meta[, "endofrecords"]
      outFrame <- rbind(outFrame, queryOut$results)
      curIndex <- curIndex + nrow(queryOut$results)
    }
    outFrame <- as.data.frame(outFrame)
    rownames(outFrame) <- outFrame$datasetKey
    outFrame
  }
  attr(occDataRaw, "datasets") <- retrieveDatasetInfo(as.character(downloadID), curlopts)
  class(occDataRaw) <- c("gbifOcc", class(occDataRaw))
  occDataRaw
}

### 2.1 ==== Autoplotting Function for 'gbifOcc' objects ====
#' @title Plot GBIF Occurrence Data
#'
#' @description Function to visualise the data associated with a GBIF occurrence
#' download.
#'
#' @param object An object that inherits from class \code{gbifOcc} (such as that
#' produced by the \code{\link{gbifOccDownload}} function).
#' @param ... Named arguments passed to constituent plotting functions.
#' Arguments that are of the form \code{points.parameter} will be passed to the
#' \code{\link[ggplot2]{geom_sf}} function that determines the plotting
#' behaviour of the occurrence records (without the \code{points.} prefix). If
#' arguments of the form \code{boundary.parameter} are provided then a boundary
#' layer is added to the plot (with data drawn from the \code{boundary.data}
#' argument) and these arguments are similarly passed to the corresponding
#' \code{\link[ggplot2]{geom_sf}} function (without the \code{boundary.}
#' prefix). Optionally an argument titled \code{titleLabelWidth} can be provided
#' which is an integer scalar giving the number of characters to aim for when
#' wrapping text associated with the titles of dataset names in legends and axis
#' labels.
#'
#' @return A \code{ggplot} layer.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link{gbifOccDownload}}, \code{\link[ggplot2]{autoplot}},
#'  \code{\link[ggplot2]{geom_sf}}
#' @export
autoplot.gbifOcc <- function(object, ...) {
  # Convert the data frame to an sf object if it isn't already
  inObject <- object
  if(!inherits(inObject, "sf")) {
    if(inherits(inObject, "data.frame")) {
      inObject <- tryCatch(sf::st_as_sf(inObject, coords = c("decimalLongitude", "decimalLatitude"), crs = gbifDefaultCRS()), error = function(err) {
        stop("unable to create sf object: ", err)
      })
      attr(inObject, "datasets") <- attr(object, "datasets")
    } else {
      stop("input is an invalid object type")
    }
  }
  # Create default parameters for the plot
  defaultPlotParams <- alist(
  )
  if(is.null(attr(inObject, "datasets"))) {
    defaultPlotParams$points.mapping = ggplot2::aes(shape = ifelse(get("occurrenceStatus") == "PRESENT", "Presence", "Absence"))
  } else {
    inObject$datasetTitle <- processDatasetTitle(attr(inObject, "datasets")[inObject$datasetKey, "datasetTitle"], retrieveOtherArgs("titleLabelWidth", ..., defaultVal = 40))
    defaultPlotParams$points.mapping = ggplot2::aes(shape = ifelse(get("occurrenceStatus") == "PRESENT", "Presence", "Absence"), colour = get("datasetTitle"))
  }
  pointPlot <- ggplot2::ggplot(data = inObject)
  # Plot the border data if any is provided
  if("boundary.data" %in% names(eval(substitute(alist(...))))) {
    pointPlot <- pointPlot + do.call(ggplot2::geom_sf, retrievePrefixArgs("boundary", ..., defaultArgs = defaultPlotParams))
  }
  # Create a plot of the presence (and absence) information
  pointPlot <- pointPlot + do.call(ggplot2::geom_sf, retrievePrefixArgs("points", ..., defaultArgs = defaultPlotParams)) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(20, 4), c("Presence", "Absence"))) +
    ggplot2::theme_classic() + ggplot2::theme(legend.title = ggplot2::element_blank())
  pointPlot
}

### 2.2 ==== Expanded Dataset Plotting Function for 'gbifOcc' objects ====

### 2.2.1 ---- Datasetplot Function ----
#' @title Plot GBIF Occurrence Data
#'
#' @description Function to visualise the data associated with a GBIF occurrence
#' download.
#'
#' @param object An object that inherits from class \code{gbifOcc} (such as that
#' produced by the \code{\link{gbifOccDownload}} function).
#' @param ... Named arguments passed to constituent plotting functions.
#' Arguments that are of the form \code{points.parameter} will be passed to the
#' \code{\link[ggplot2]{geom_sf}} function that determines the plotting
#' behaviour of the occurrence records (without the \code{points.} prefix). If
#' arguments of the form \code{boundary.parameter} are provided then a boundary
#' layer is added to the plot (with data drawn from the \code{boundary.data}
#' argument) and these arguments are similarly passed to the corresponding
#' \code{\link[ggplot2]{geom_sf}} function (without the \code{boundary.}
#' prefix). Argument with the prefixes \code{yearbar.} and \code{datasetbar.}
#' are passed to the \code{\link[ggplot2]{geom_bar}} functions that plot the
#' dataset occurrences counted by year and source dataset respectively.
#' Arguments with the \code{colour.} prefix are passed to the
#' \code{\link[ggplot2]{scale_colour_discrete}} and
#' \code{\link[ggplot2]{scale_fill_discrete}} functions that define the colour
#' palettes in the constituent plots. Finally, arguments with \code{arrange.}
#' prefix are passed to the \code{\link[gridExtra]{arrangeGrob}} function that
#' arranges the constituent plots.
#' @param titleLabelWidth An integer scalar giving the number of characters to
#' aim for when wrapping text associated with the titles of dataset names in
#' legends and axis labels.
#'
#' @return A \code{datasetplotGrob} object.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link{gbifOccDownload}}, \code{\link[ggplot2]{autoplot}},
#'  \code{\link[ggplot2]{geom_sf}}, \code{\link[gridExtra]{arrangeGrob}},
#'  \code{\link[ggplot2]{scale_colour_discrete}},
#'  \code{\link[ggplot2]{scale_fill_discrete}}
#' @export
datasetplot <- function(object, ..., titleLabelWidth = 40) {
  inObject <- object
  grobList <- list()
  # Create default parameters for the plot
  defaultPlotParams <- alist(
    arrange.name = paste("datasetplot", deparse(substitute(object)), sep = "_"),
    arrange.heights = ggplot2::unit(c(0.7, 0.3), "npc"),
    arrange.layout_matrix = matrix(1:2, 2, 1)
  )
  if(!is.null(attr(inObject, "datasets"))) {
    defaultPlotParams$yearbar.mapping <- ggplot2::aes(fill = get("datasetTitle"))
    defaultPlotParams$yearbar.colour <- "white"
    defaultPlotParams$datasetbar.colour <- NA
    defaultPlotParams$arrange.widths <- ggplot2::unit(c(0.4, 0.6), "npc")
    defaultPlotParams$arrange.layout_matrix <- matrix(c(1, 1, 2, 3), 2, 2)
    inObject$datasetTitle <- processDatasetTitle(attr(inObject, "datasets")[inObject$datasetKey, "datasetTitle"], titleLabelWidth)
  }
  # Create the dataset point plot
  occPoints <- do.call(autoplot.gbifOcc, append(list(object = inObject, titleLabelWidth = titleLabelWidth), eval(substitute(alist(...)))))
  # Create bar plot of samples over each year
  yearBar <- ggplot2::ggplot(data = inObject, ggplot2::aes(x = get("year"))) + do.call(ggplot2::geom_bar, retrievePrefixArgs("yearbar", ..., defaultArgs = defaultPlotParams)) +
    ggplot2::ylab("Count") + ggplot2::xlab("Year") + ggplot2::theme_classic()
  if(!is.null(attr(inObject, "datasets"))) {
    # Remove the colour legend from the points plot (and use any provided colour parameter to overwrite the plotting colour)
    occPoints <- occPoints + ggplot2::guides(colour = "none") + do.call(ggplot2::scale_colour_discrete, retrievePrefixArgs("colour", ..., defaultArgs = defaultPlotParams))
    # Create a bar plot of the datasets
    datasetBar <- ggplot2::ggplot(data = inObject, ggplot2::aes(y = get("datasetTitle"), fill = get("datasetTitle"), colour = get("datasetTitle"))) +
      do.call(ggplot2::geom_bar, retrievePrefixArgs("datasetbar", ..., defaultArgs = defaultPlotParams)) +
      do.call(ggplot2::scale_fill_discrete, retrievePrefixArgs("colour", ..., defaultArgs = defaultPlotParams)) +
      ggplot2::guides(fill = "none", colour = "none") + ggplot2::ylab("") + ggplot2::xlab("Count") + ggplot2::theme_classic() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7)))
    # Adapt the year bar plot to include information about datasets
    yearBar <- yearBar + ggplot2::guides(fill = "none") + do.call(ggplot2::scale_fill_discrete, retrievePrefixArgs("colour", ..., defaultArgs = defaultPlotParams))
    # Add the dataset bar to the beginning of the grob list
    grobList <- list(datasetBar)
  }
  grobList <- lapply(X = append(grobList, list(occPoints, yearBar)), FUN = ggplot2::ggplotGrob)
  # Call the grob arrange function to organise all the plots
  outGrob <- do.call(gridExtra::arrangeGrob, append(grobList, retrievePrefixArgs("arrange", ..., defaultArgs = defaultPlotParams)))
  class(outGrob) <- c("datasetplotGrob", class(outGrob))
  outGrob
}

### 2.2.2 ---- Overloaded print method ----
#' @title Print Function for Dataset Plot Objects
#'
#' @description Defined method for the automatic plotting of dataset
#' plot objects created by the \code{\link{datasetplot}} function
#'
#' @param x A \code{datasetplotGrob} object to plot
#' @param ... A set of parameter that can include: \code{newpage}, a logical
#' scalar that, if \code{TRUE}, plots the object on a new page
#'
#' @return A \code{datasetplotGrob} object
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
print.datasetplotGrob <- function(x, ...) {
  if(retrieveOtherArgs("newpage", ..., defaultVal = TRUE)) {
    grid::grid.newpage()
  }
  grid::grid.draw(x)
  invisible(x)
}
