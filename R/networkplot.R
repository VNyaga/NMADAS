#' Network plot
#' @description  Draw a network plot.
#' @param data A data.frame with atleast two variables i.e study and test identifier.
#' @param SID A text indicating the name of the variable that identifies the different studies.
#' @param TID A text indicating the name of the variable that identifies the different tests.
#' @param ... additional options. See \link[pcnetmeta]{nma.networkplot} for more details.
#' @examples
#' data(demodata)
#' networkplot(data = demodata,
#'             SID = "study",
#'             TID = "Test")
#' @export
#' @author Victoria N Nyaga

networkplot.nmadas <- function(data,
  S.ID,
  T.ID,
  alphabetic = FALSE,
  edge.col = "grey",
  node.col = "orange",
  ...){

  df <- prepdata.nmadas(data,
     S.ID = S.ID,
     T.ID = T.ID)

  #if (grDevices::dev.interactive()) grDevices::dev.new()

  pcnetmeta::nma.networkplot(
     s.id = SID,
     t.id = TID,
     data = df@data,
     trtname = df@labels,
     alphabetic = alphabetic,
     edge.col = edge.col,
     node.col = node.col,
      ...)
  }




