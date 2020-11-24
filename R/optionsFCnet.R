# Variable, global to package's namespace.
# This function is not exported to user space and does not need to be documented.
opt= settings::options_manager(
  colorPaletteDefault= c("darkslateblue", "royalblue4",
                       "royalblue1", "cyan2",
                       "white",
                       "darkgoldenrod1", "darkgoldenrod3",
                       "red3", "firebrick4")

  )

# User function that gets exported:

#' Set or get options for FCnet
#'
#' @param colorPaletteDefault The default color palette for plotting matrices (a vector of colors).

#' @export

optionsFCnet <- function(...){

  # protect against the use of reserved words.
  settings::stop_if_reserved(...)
  # return/change
  opt(...)

}
