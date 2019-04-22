#' Get input df for UpSetR
#'
#' @details Very much the same as \code{UpSetR::fromList()}, the only difference
#' being that I am adding rownames to the output.
#'
#' @param input list
#' @return data.frame where the column names correspond to the names of the vectors
#' of the list
#' @examples
#' listInput <- list(one = c(1, 2, 3, 5, 7, 8, 11, 12, 13), 
#'                  two = c(1, 2, 4, 5, 10),
#'                  three = c(1, 5, 6, 7, 8, 9, 10, 12, 13))
#' ups_input <- make_upsetr_input(listInput)
make_upsetr_input <- function(input){
  # get a vector of all entries
  universe <- unique(unlist(input))
  
  #
  data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(universe, x)) # NA will be introduced for every no-match
    }))
  
  data[is.na(data)] <- as.integer(0) # mark every non-match with a zero
  data[data != 0] <- as.integer(1) # every match gets a one
  # get the correct shape of the data.frame
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  row.names(data) <- universe
  return(data)
}

