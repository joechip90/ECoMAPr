# 1.1 ====  runSDM Function ====
#' Run a distribution model
#'
#' @description
#' This function runs a single species distribution model (or series of species
#' distribution models if the \code{group} parameter is not \code{NULL}).
#' The function is also responsible for collating the relevant metadata for the
#' analysis.
#'
#' @param sdmFunc A function that runs the species distribution model.
#' @param occData An \code{\link[sf]{sf}} object containing the occurrence data
#' of the taxa to be modelled.
#' @param group A character scalar that is used identify the column name in
#' \code{occData} that is used to assign the different batches that the species
#' distribution model is going to be applied.  If \code{NULL} (the default),
#' then only a single SDM is run on the entire occurrence data. If a valid
#' column name is given then a seperate SDM is run for each subset of the
#' occurrence data based on unique values in that column.
#' @param outLoc A character scalar containing the location to store the output
#' model objects. If \code{group} is \code{NULL} then the model outputs are
#' stored directly in the directory specified. Otherwise, a subfolder is
#' created for each SDM.
#' @param numCores Integer vector containing the number of processes to be run
#' simultaneously. If \code{group} is not \code{NULL} then the first element of
#' the integer vector is used to divide resources for processing each of the
#' SDMs for the groups.  All unused elements of the integer vector are passed to
#' \code{sdmFunc} (and may be used there to control parallelisation within the
#' SDM function). If \code{NULL}, then defaults to the number of cores
#' available.
#' @param ... Named arguments to be passed to \code{sdmFunc}.
#'
#' @return A logical vector giving the exit status of the model functions.
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @seealso \code{\link[terra]{sds}}, \code{\link[terra]{rast}},
#' \code{\link[base]{strptime}}
#' @export
#'
runSDM <- function(sdmFunc, occData, group = NULL, outLoc = tempdir(), numCores = NULL, ...) {
  # 1.1.1 ---- Process the input arguments ----
  # Process the output location parameter
  parOutLoc <- tryCatch(as.character(outLoc), error = function(err) {
    stop("error encountered processing the ouput location argument: ", err)
  })
  if(is.null(parOutLoc) || length(parOutLoc) <= 0) {
    parOutLoc <- tempdir()
  } else if(length(parOutLoc) > 1) {
    warning("output location has a length greater than one: only the first element will be used")
    parOutLoc <- parOutLoc[1]
  }
  if(is.na(parOutLoc)) {
    parOutLoc <- tempdir()
  }
  # Process the occurrence data
  parOccData <- tryCatch(sf::st_as_sf(occData), error = function(err) {
    stop("error encountered processing the occurrence data argument: ", err)
  })
  # Process the grouping parameter
  parGroup <- tryCatch(as.character(group), error = function(err) {
    stop("error encountered processing the grouping argument: ", err)
  })
  if(is.null(parGroup) || length(parGroup) <= 0) {
    groupVals <- rep(basename(parOutLoc), nrow(parOccData))
    parOutLoc <- dirname(parOutLoc)
  } else if(length(parGroup) == 1 && parGroup %in% names(parOccData)) {
    groupVals <- tryCatch(as.character(parOccData[, parGroup]), error = function(err) {
      stop("error encountered processing the grouping argument: unable to convert selected column to character vector(", err, ")")
    })
  } else {
    groupVals <- parGroup[0:(nrow(parOccData) - 1) %% length(parGroup) + 1]
  }
  if(anyNA(groupVals)) {
    warning("NA values found in grouping variable: associated occurrence rows will be ignored")
    parOccData <- parOccData[!is.na(groupVals), ]
    groupVals <- groupVals[!is.na(groupVals)]
  }
  uniqueGroups <- unique(groupVals)
  # Process the number of cores parameter
  if(is.null(numCores) || length(numCores) <= 0) {
    parNumCores <- parallelly::availableCores()
  } else {
    parNumCores <- tryCatch(as.integer(numCores), error = function(err) {
      stop("error encountered processing the number of cores argument: ", err)
    })
    if(length(parNumCores) <= 0 || is.na(parNumCores[1]) || parNumCores[1] <= 0) {
      warning("invalid number of cores provided as argument: setting cores to 1")
      parNumCores[1] <- 1
    }
  }
  # Collect together the other parameters
  otherParams <- eval(substitute(alist(...)))
  # 1.1.2 ---- Initialise parallelisation ----
  parNumCores[1] <- min(parNumCores[1], length(uniqueGroups))
  # Distribute the groups to process between the cores
  groupsInCores <- lapply(1:parNumCores[1], FUN = function(coreNum, uniqueGroups, totCores) {
    uniqueGroups[0:(length(uniqueGroups) - 1) %% totCores + 1 == coreNum]
  }, uniqueGroups = uniqueGroups, totCores = parNumCores[1])
  # 1.1.3 ---- Function to process groups ----
  runGroups <- function(groupsToProc, parOutLoc, parOccData, groupVals, otherParams) {
    sapply(X = groupsToProc, FUN = function(curGroup, parOutLoc, parOccData, groupVals, otherParams) {
      # Filter out the occurrence data entries that belong in the current group
      curOccData <- parOccData[groupVals == curGroup, ]
      # Create an output folder to store the model outputs
      curOutLoc <- file.path(parOutLoc, curGroup)
      if(dir.exists(curOutLoc)) {
        unlist(curOutLoc, recursive = TRUE)
      }
      if(!dir.create(curOutLoc)) {
        stop("unable to create output folder at ", curOutLoc)
      }
    }, parOutLoc = parOutLoc, parOccData = parOccData, groupVals = groupVals, otherParams = otherParams)
  }
}
