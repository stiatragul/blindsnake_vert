# func_plotTreeboxplot.R
# This is a working version from: http://blog.phytools.org/2016/07/new-function-to-print-box-plot-next-to.html

plotTree.boxplot<-function(tree,x,args.plotTree=list(),
                           args.boxplot=list()){
  cw<-reorder(tree)
  if(!is.list(x)&&class(x)!="formula"){
    obj<-setNames(
      lapply(cw$tip.label,function(x,y) y[which(names(y)==x)],
             y=x),cw$tip.label)
  } else obj<-x
  if(class(x)=="formula") 
    args.boxplot$formula<-obj else args.boxplot$x<-obj
    args.boxplot$horizontal<-TRUE
    args.boxplot$axes<-FALSE
    args.boxplot$names.arg<-""
    args.boxplot$xlim<-c(1,Ntip(cw))
    args.boxplot$ylim<-NULL
    if(is.null(args.boxplot$space)) args.boxplot$space<-0.7
    if(is.null(args.boxplot$mar)) 
      args.boxplot$mar<-c(5.1,0,2.1,1.1)
    else args.boxplot$mar[2]<-0.1
    args.plotTree$tree<-cw
    if(is.null(args.plotTree$mar)) 
      args.plotTree$mar<-c(5.1,1.1,2.1,0)
    else {
      args.plotTree$mar[4]<-0
    }
    if(args.plotTree$mar[1]!=args.boxplot$mar[1])
      args.plotTree$mar[1]<-args.boxplot$mar[1]
    if(args.plotTree$mar[3]!=args.boxplot$mar[3])
      args.plotTree$mar[3]<-args.boxplot$mar[3]
    if(is.null(args.plotTree$ftype)) args.plotTree$ftype<-"i"
    if(is.null(args.plotTree$lwd)) args.plotTree$lwd<-1
    par(mfrow=c(1,2))
    do.call(plotTree,args.plotTree)
    par(mar=args.boxplot$mar)
    ii<-which(names(args.boxplot)%in%c("formula","x"))
    args.boxplot<-c(args.boxplot[ii],args.boxplot[-ii])
    obj<-do.call(boxplot,args.boxplot)
    axis(1)
    if(!is.null(args.boxplot$xlab)) title(xlab=args.boxplot$xlab)
    else title(xlab="x")
    invisible(obj)
}

# plotTree.boxplot <- function(tree, x, args.plotTree = list(), args.boxplot = list()) {
#   cw <- reorder(tree)
#   if (!is.list(x) && class(x) != "formula") {
#     obj <- setNames(
#       lapply(cw$tip.label, function(x, y) y[which(names(y) == x)],
#              y = x), cw$tip.label)
#   } else obj <- x
#   if (class(x) == "formula") 
#     args.boxplot$formula <- obj 
#   else args.boxplot$x <- obj
#   
#   args.boxplot$horizontal <- TRUE
#   args.boxplot$axes <- FALSE
#   args.boxplot$names.arg <- ""
#   args.boxplot$xlim <- c(1, Ntip(cw))
#   
#   if (is.null(args.boxplot$ylim)) 
#     args.boxplot$ylim <- NULL  # Default to NULL
#   
#   if (is.null(args.boxplot$space)) 
#     args.boxplot$space <- 0.7
#   
#   if (is.null(args.boxplot$mar)) 
#     args.boxplot$mar <- c(5.1, 0, 2.1, 1.1)
#   else 
#     args.boxplot$mar[2] <- 0.1
#   
#   args.plotTree$tree <- cw
#   
#   if (is.null(args.plotTree$mar)) 
#     args.plotTree$mar <- c(5.1, 1.1, 2.1, 0)
#   else 
#     args.plotTree$mar[4] <- 0
#   
#   if (args.plotTree$mar[1] != args.boxplot$mar[1]) 
#     args.plotTree$mar[1] <- args.boxplot$mar[1]
#   
#   if (args.plotTree$mar[3] != args.boxplot$mar[3]) 
#     args.plotTree$mar[3] <- args.boxplot$mar[3]
#   
#   if (is.null(args.plotTree$ftype)) 
#     args.plotTree$ftype <- "i"
#   
#   if (is.null(args.plotTree$lwd)) 
#     args.plotTree$lwd <- 1
#   
#   par(mfrow = c(1, 2))
#   do.call(plotTree, args.plotTree)
#   par(mar = args.boxplot$mar)
#   ii <- which(names(args.boxplot) %in% c("formula", "x"))
#   args.boxplot <- c(args.boxplot[ii], args.boxplot[-ii])
#   
#   # Check if ylim is specified and update if not NULL
#   if (!is.null(args.boxplot$ylim)) {
#     obj <- do.call(boxplot, c(list(x), args.boxplot))
#   } else {
#     # Apply default boxplot with default ylim
#     obj <- do.call(boxplot, args.boxplot)
#   }
#   
#   axis(1)
#   if (!is.null(args.boxplot$xlab)) 
#     title(xlab = args.boxplot$xlab)
#   else 
#     title(xlab = "x")
#   
#   invisible(obj)
# }
