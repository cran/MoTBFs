#' Plots for \code{'motbf'} objects
#'
#' Draws an \code{'motbf'} function.
#' 
#' @param x An object of class \code{'motbf'} or one of its subclasses (\code{'univmotbf'}, \code{'piecewisemop'}, \code{'motbf.fit.node'}, or \code{'jointmotbf'}).
#' @param xlim,ylim Numeric vectors of length 2 specifying x and y axis limits. By default \code{0:1}. Used in \code{'univmotbf'} and \code{'piecewisemop'}.
#' @param type Character string specifying the plot type, as for \link{plot}. For \code{'univmotbf'} defaults to \code{'l'} (line); for \code{'jointmotbf'} defaults to \code{"contour"} (another option is \code{"perspective"}).
#' @param add Logical; if \code{TRUE}, adds the plot to the existing active graphics device. Used only in \code{'univmotbf'}.
#' @param panels Logical; if \code{TRUE}, displays plots in separate panels. Used only in \code{'motbf.fit.node'}.
#' @param ranges A \code{"numeric"} matrix containing the domain of the variables, by columns, which is used to specify the plotting range. Used only in \code{'jointmotbf'}.
#' @param orientation A \code{"numeric"} vector indicating the perpective of the plot in degrees. By default, it is set to \code{(5,-30)}. Used only in \code{'jointmotbf'}.
#' @param data An object of class \code{"data.frame"} containing two columns only.
#' This argument is used to draw the points over the main plot. By default, it is set to \code{NULL}. Used only in \code{'jointmotbf'}.
#' @param filled A logical argument; it is only used if \code{type = "contour"}.
#' is active. By default, it is \code{TRUE}, so filled contours are plotted. Used only in \code{'jointmotbf'}.
#' @param ticktype A \code{"character"} string, either \emph{simple} or \emph{detailed}. By default, it is set to \code{"simple"},
#'  which draws just an arrow parallel to the axis to indicate direction of increase. 
#'  In contrast, \code{"detailed"} draws normal ticks. This argument is only used in the \code{"perspective"} plot. Used only in \code{'jointmotbf'}.
#' @param main Title string for the plot. Used in \code{'jointmotbf'}.
#' @param \dots Further arguments to be passed as for \link{plot}.
#' @return A plot of the specificated function.
#' @name plot.motbf
#' @rdname plot.motbf
#' @exportS3Method graphics::plot motbf
#' @examples
#'
#'## 1. univmotbf - Example for univariate distributions 
#'## Data
#'X <- rexp(2000)
#'
#'f1 <- univMoTBF(X, POTENTIAL_TYPE = "MOP")
#'f2 <- univMoTBF(X, POTENTIAL_TYPE = "MTE", maxParam = 10)
#'
#'## Plots
#'plot(NULL, xlim = range(X), ylim = c(0,0.8), xlab="X", ylab="density")
#'plot(f1, xlim = range(X), col = 1, add = TRUE)
#'plot(f2, xlim = range(X), col = 2, add = TRUE)
#' 
#' 
#'## 2. jointmotbf - Example for join distributions 
#' ## Data
#' X <- data.frame(rnorm(500), rnorm(500))
#' 
#' ## Joint function
#' dim <- c(3,3) 
#' P <- jointmotbf.fit(X, dimensions = dim)
#' 
#' ## Plots
#' plot(P)
#' plot(P, type = "perspective", orientation = c(90,0))
#' 
plot.motbf <- function(x, ...){
  NextMethod("plot")
}

#' @rdname plot.motbf
#' @exportS3Method graphics::plot univmotbf
plot.univmotbf <-  function(x, xlim=NULL, ylim=NULL, type='l', add = FALSE,...){
  if(is.null(xlim)){
    xlim = unlist(x$Domain)
  }
  var = colnames(x$Domain)
  
  plot(as.function(x), xlim=xlim, ylim=ylim, type=type, add = add, xlab = var, ylab = 'Density',...)
}



#' @rdname plot.motbf
#' @exportS3Method graphics::plot piecewisemop
plot.piecewisemop <- function(x, xlim=NULL, ylim=NULL, ...){
  class(x) =  "piecewisemop" 
  dom = sapply(x, '[[', 'Domain')

  if(is.null(xlim)){
    xlim = c(min(dom), max(dom))
  }

  if(is.null(ylim)){
    ymax = c()
    var = getMotbfVar(x[[1]])
    for(i in 1:length(x)){
      f = eval(parse(text = paste("f <- function(",var,")",x[[i]]$Function)))
      domi = dom[,i]
      ymax[i] = max(f(seq(domi[1], domi[2], 0.01)))
    }

    ylim = c(0, round(max(ymax)+0.1,1))
  }

  plot(1, type = 'n', ylab = 'Density', xlab = getMotbfVar(x[[1]]), xlim = xlim, ylim = ylim)


  for(i in 1:length(x)){
    plot(x[[i]], add = TRUE, ...)
  }
}

#' @rdname plot.motbf
#' @exportS3Method graphics::plot motbf.fit.node
plot.motbf.fit.node <- function(x, panels = TRUE, ...){
  opar <- par(no.readonly =TRUE)       
  on.exit(par(opar)) 
  
  if(x$type == 'Continuous'){
    foo = x$functions
  }else{
    foo = collapseDiscreteCPD(x)    
  }
  
  n = ncol(foo)
  r = nrow(foo)
  
  var = colnames(foo)[n-1]
  # browser()
  # Get parent values to create a main title for the plot
  if(!is.null(x$parents)){
    if(panels){
      # par(mfrow=c(ceiling(r/3),r))
      columnas = ceiling(sqrt(r))
      filas =  ceiling(r/columnas)
      par(mfrow = c(filas,columnas))
    }
    
    parents = foo[,1:(n-2), drop = FALSE]
    parent_names = list()
    
    for(i in 1:(n-2)){
      p = colnames(parents)[i]
      if(is.list(parents[,i])){
        
        parent_names[[i]]=sapply(1:r, function(k) paste(round(parents[,i][[k]][1],4),'<', p, '<',round(parents[,i][[k]][2],4)))
      }else{
        parent_names[[i]]=paste(p, '=',parents[,i])
      }
    }
    
    main = apply(do.call(cbind, parent_names), 1, paste, collapse = '\n ')
    
  }else{
    main = paste('Fitted distribution of', var)
  }
  
  
  if(x$type == 'Continuous'){
    
    if(!exists("xlim")){
      xlim = as.vector(foo[,(n-1)][[1]])
    }
    
    if(!exists("ylim")){
      ymax = c()
      
      mops = foo[,n]
      for(i in 1:length(mops)){
        f = eval(parse(text = paste("f <- function(",var,")",mops[[i]]$Function)))
        domi = xlim
        ymax[i] = max(f(seq(domi[1], domi[2], 0.01)))
      }
      
      ylim = c(0, round(max(ymax)+0.1,1))
    }
    # browser()
    for(i in 1:r){
      plot(foo[,n][[i]], main = main[i])
    }
  }else{
    for(i in 1:r){
      graphics::barplot(unlist(foo[i,n]), names = unlist(foo[i,(n-1)]), xlab = var, 
              main = main[i], ylim = c(0, 1))
    }


  }

  



}


#' @rdname plot.motbf
#' @exportS3Method graphics::plot jointmotbf
plot.jointmotbf <- function(x, type="contour", ranges=NULL, orientation=c(5,-30), 
                            data=NULL, filled=TRUE, ticktype="simple", main = NULL, ...)
{
  opar <- par(no.readonly =TRUE)       
  on.exit(par(opar)) 
  
  varname <- attr(x$Domain, "dimnames")[[2]]
  
  if(length(getMotbfVar(x))!=2) stop("It is not possible plotting a joint function with more than 2 variables.")
  if(is.null(ranges)) ranges <- x$Domain
  #{ranges <- x$Domain; ranges[1,] <- ranges[1,]+2E-1; ranges[2,] <- ranges[2,]-2E-1}
  
  X <- seq(min(ranges[,varname[1]]), max(ranges[,varname[1]]), length = 20)
  Y <- seq(min(ranges[,varname[2]]), max(ranges[,varname[2]]), length = 25)
  Z <- outer(X, Y, as.function(x))
  Z[Z<=0]=1.0E-10
  
  if(is.null(main)){
    main = paste('Joint MoTBF of', paste(varname, collapse = ' and '))
  }
  if(!is.null(data)){
    if(!all(varname %in% colnames(data))){
      stop('Cannot find variable ', varname[which(!(varname %in% colnames(data)))], ' in dataset')
    }
    data = data[,varname]
  }
  if(type == "perspective"){    
    nrz <- nrow(Z)
    ncz <- ncol(Z)
    ## Create a function interpolating colors in the range of specified colors
    jet.colors <- colorRampPalette( c("blue", "green") )
    
    ## Generate the desired number of colors from this palette
    nbcol <- 100; color <- jet.colors(nbcol)
    
    ## Compute the z-value at the facet centres
    zfacet <- Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz]
    
    ## Recode facet z-values into color indices
    facetcol <- cut(zfacet, nbcol)
    persp(X, Y, Z, col = color[facetcol], phi = orientation[1], theta = orientation[2], ticktype=ticktype,
          xlab = varname[1], ylab = varname[2], zlab = '', main = main, ...) 
  }
  if(type == "contour"){
    if(!filled){
      if(!is.null(data)) plot(data, xlab = varname[1], ylab = varname[2], xlim=ranges[,1], ylim=ranges[,2], main = main, ...)
      if(is.null(data)) plot(NULL, xlim=ranges[,1], ylim=ranges[,2], xlab="", ylab="")
      contour(X,Y,Z, col=terrain.colors(15), xlab = varname[1], ylab = varname[2], add=TRUE)
      abline(h =round(ranges[,1])[1]:round(ranges[,1])[2] , v =round(ranges[,2])[1]:round(ranges[,2])[2],
             col = "gray", lty = 2, lwd = 0.1)   
    }
    if(filled){
      zlim <- range(Z, finite=TRUE)
      zlim[1] <- 0
      nlevels <- 20
      levels <- pretty(zlim, nlevels)
      nlevels <- length(levels)
      color <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"))
      filled.contour(X,Y,Z,nlevels=nlevels,levels=levels,
                     las=1,
                     col=color(nlevels),                           
                     plot.title = {title(xlab = varname[1], ylab = varname[2], cex.lab = 1.5, main = main)},...
      )
      if(!is.null(data)){
        mar.orig <- par("mar")
        w <- (3 + mar.orig[2]) * par("csi") * 2.54
        layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
        points(data)
        par(mfrow=c(1,1))
      }
    }
  }
}




#' Plot Conditional Functions
#' 
#' Plot conditional MoTBF densities.
#' 
#' @param conditionalFunction the output of function \link{conditionalMethod}. 
#' A list containing the the interval of the parent and the final conditional density (MTE or MOP).
#' @param data An object of class \code{data.frame}, corresponding to the dataset used to fit the conditional density.
#' @param nameChild A \code{character} string, corresponding to the name of the child variable in the conditional density. By default, it is \code{NULL}.
#' @param points A logical value. If \code{TRUE}, the sample points are overlaid.
#' @param color If not specified, a default palette is used. 
#' @param ... Additional graphical parameters passed to filled.contour().
#' @details If the number of parents is greater than one, then the error message 
#' "It is not possible to plot the conditional function." is reported.
#' @return A plot of the conditional density function.
#' @seealso \link{conditionalMethod}
#' @export
#' @examples
#' ## Data
#' X <- rnorm(1000)
#' Y <- rnorm(1000, mean=X)
#' data <- data.frame(X=X,Y=Y)
#' cov(data)
#' 
#' ## Conditional Learning
#' parent <- "X"
#' child <- "Y"
#' intervals <- 5
#' potential <- "MTE"
#' P <- conditionalMethod(data, nameParents=parent, nameChild=child, 
#' numIntervals=intervals, POTENTIAL_TYPE=potential)
#' plotConditional(conditionalFunction=P, data=data)
#' plotConditional(conditionalFunction=P, data=data, points=TRUE)
#' 
plotConditional <- function(conditionalFunction, data, nameChild=NULL, points=FALSE, color=NULL,...)
{
  opar <- par(no.readonly =TRUE)       
  on.exit(par(opar)) 
  
  ## Define the color
  if(is.null(color)) color <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"))
  
  nameParent <- conditionalFunction[[1]]$parent
  if(length(nameParent)==1){ 
    if(is.null(nameChild)){
      if(ncol(data)==2) nameChild <- colnames(data)[which(colnames(data)!=nameParent)]
      else (stop("The name of the child variable is needed because the number of columns in the dataset is bigger than two."))
    }else{
      if(ncol(data)==2)
        if(colnames(data)[which(colnames(data)!=nameParent)]!=nameChild)
          (stop("The name of the child variable has not been found in the dataset."))
    }
    X <- data[, nameParent]
    Y <- data[, nameChild]
    
    xgrid <- seq(min(Y),max(Y),length.out=50) 
    ygrid <- seq(min(X),max(X),length.out=50)
    griddata <- as.data.frame(expand.grid(xgrid,ygrid))
    D <- c()
    for(i in 1:length(conditionalFunction)){
      j <- which((griddata[,2]>=conditionalFunction[[i]]$interval[1])&(griddata[,2]<=conditionalFunction[[i]]$interval[2]))
      gdi <- griddata[j,]
      D <- c(D,as.function(conditionalFunction[[i]]$Px)(gdi[,1]))
    }
    d <- matrix(D, nrow = length(xgrid), ncol = length(ygrid), byrow=FALSE)
    zlim <- range(d, finite=TRUE)
    zlim[1] <- 0
    nlevels <- 20
    levels <- pretty(zlim, nlevels)
    nlevels <- length(levels)
    
    filled.contour(x = ygrid, y = xgrid, z = t(d), nlevels = nlevels,
                   levels = levels, col = color(nlevels),                           
                   plot.title = {title(xlab = nameParent, ylab = nameChild, cex.lab = 2)},
                   ...
    )
    if(points){
      mar.orig <- par("mar")
      w <- (3 + mar.orig[2]) * par("csi") * 2.54
      layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
      points(data[,c(nameParent, nameChild)])
      par(mfrow=c(1,1))
    }
  }else{
    (stop("It is not possible to plot the conditional function."))
  }
}


# definir colores
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
