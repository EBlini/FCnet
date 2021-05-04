
#' Read and load Functional Connectivity matrices or volumes
#'
#' This function is a helper function to quickly read and import in R functional
#' connectivity matrices stored as standard .csv files or lesion/disconnection
#' data stored as nifti files.
#' Files (.csv) must adhere to a strict format:
#' .csv files with no header nor row names.
#' It is of PARAMOUNT importance that FC matrices/volumes are correctly ordered
#' in their folder, according to participants' IDs (and therefore the behavioral
#' score(s) to predict). Note that characters are ordered alphabetically by
#' default in R, thus you may want to ensure that participants with low IDs
#' (e.g. # 5) are identified by file names preserving this order (e.g. FC 05.csv).
#'

#' @param path Path to the input files. If empty, it prompts the user
#' to select it manually via a selection menu. Only FC matrices or
#' nifti files must be included in the folder. The function performs
#' a cheap guess to understand which input is passed.
#' Files other than those are disregarded.
#' @param files An optional list of file names, e.g. as created by
#' `list.files()`. If this is supplied, path is overwritten as
#' the file addressed are assumed to be relative.
#' @param quick For csv files: the default uses base R, which may be slow if the FC matrices
#' are too many. Option "data.table" uses `data.table::fread()`, which is sensibly
#' faster when FC matrices are many; `data.table` must be installed. Option "parallel" uses
#' `future.apply::future_lapply()`: `future.apply` must be installed, your
#' machine should have multiple cores available for use, and threads should be defined explicitly
#' beforehand by the user (e.g. by calling `plan(multisession)`).

#' @return A list of lists including the FC matrices/volumes,
#' stored as data.frames or arrays.

#' @export

loadFC= function(path= choose.dir(),
                 files= NULL,
                 quick= c("default", "data.table", "parallel")){

  #match arguments
  quick= match.arg(quick)

  #read all files
  if(is.null(files)){
    files= list.files(path)
  } else {path= ".\\"}


  #if first in line is csv then read csvs
  if(grepl(".csv", files[1], fixed = T)){

    #disregard non-csv
    files= files[grepl(".csv", files, fixed = T)]

    #read all files
    if(quick== "default"){

      FC= lapply(files, function(x){

        #get proper path regardless of OS
        np= normalizePath(file.path(path, x),
                          winslash = "\\")

        matrix= read.csv(np, header= F)

        return(matrix)

      }) #end lapply
    }# end if "default"

    if(quick== "data.table"){

      FC= lapply(files, function(x){

        #get proper path regardless of OS
        np= normalizePath(file.path(path, x),
                          winslash = "\\")

        matrix= data.table::fread(np, header= F)
        matrix= as.data.frame(matrix)

        return(matrix)

      }) #end lapply
    }# end if "data.table"

    if(quick== "parallel"){

      FC= future.apply::future_lapply(files, function(x){

        #get proper path regardless of OS
        np= normalizePath(file.path(path, x),
                          winslash = "\\")

        matrix= read.csv(np, header= F)

        return(matrix)

      }, future.seed = T) #end future_lapply
    }# end if "parallel"

  } else {# end if .csv

    FC= lapply(files, function(x){

      #get proper path regardless of OS
      np= normalizePath(file.path(path, x),
                        winslash = "\\")

      matrix= oro.nifti::readNIfTI(np)

      matrix= oro.nifti::img_data(matrix)

      return(matrix)

    }) #end lapply
  } #end if not csv


  #return FC lists
  return(FC)

}
