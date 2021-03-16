  x <- c(1.2, 0.5, 0.103, 7, 0.1501)

  Percent <- function(Data, Type, digits = 2,
                      format = "f", ...) {
    # Create user-defined function
    if (Type == "Value") {
      paste0(formatC(Data/(1)), "%")
    } else {
      paste0(formatC((Data/sum(Data))*100, format = format, 
                     digits = digits, ...), "%")
    }
  }

  Percent(Data = x, Type = "Frame")  #Value, Frame
  