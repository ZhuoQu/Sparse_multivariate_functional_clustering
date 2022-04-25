
subsetlist <- function(myList, elementNames) {
  lapply(elementNames, FUN = function(x) { myList[[x]] })
}
