# Histogram function >>>
panel.hist.add <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE) # set breaks here >
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

# Histogram function >>>
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = TRUE) # set breaks here >
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}

# Smooth Scatterplot function, add TRUE >>>
panel.smooth.add <- function(x, y){
  smoothScatter(x,y, add=TRUE)
  r <- round(cor(x, y, use = "complete.obs"), digits=4)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.4, 0.85, txt)
  ylim=c(20,40)
  xlim=c(20,40)
  abline(0,1, col = "grey")
}

# Smooth Scatterplot function, add FALSE >>>
panel.smooth <- function(x, y){
  smoothScatter(x,y, add=FALSE)
  r <- round(cor(x, y, use = "complete.obs"), digits=4)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.4, 0.85, txt)
  ylim=c(20,40)
  xlim=c(20,40)
  abline(0,1, col = "grey")
}
