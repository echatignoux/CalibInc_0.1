##' Internal to deal hatching of NS zones
##'
##'
##' @param poly : polygon
##' @param l : hatching lines spacement
##' @param shap : hatching lines orientation in degree
##' @return
##' @author Edouard Chatignoux
##' @export
##' @keywords internal
##' @examples
hach.poly<-function(poly,l=200,shap=-pi/4){
    require(rgeos)
    rg<-apply(summary(poly)$bbox,1,diff)
    x<-c(2*tan(-shap),2*tan(shap))*rg[1]+summary(poly)$bbox[1,1]
    shad.lines<-sp::SpatialLines(lapply(seq(x[1],x[2],l=l),function(x0)
      sp::Lines(sp::Line(cbind(x,
                               tan(shap)*(x-x0)+summary(poly)$bbox[2,1])
                         ),ID=x0)))
    do.call("rbind",
            lapply(poly@polygons,
                   function(dep){
                       a<-gIntersection(shad.lines,sp::SpatialPolygons(list(dep)))
                       a@lines[[1]]@ID<-dep@ID
                       do.call("rbind",
                               lapply(1:length(slot(a,"lines")[[1]]@Lines),
                                      function(l)
                                      data.frame(id=slot(a,"lines")[[1]]@ID,
                                                 l=paste(slot(a,"lines")[[1]]@ID,l,sep="."),
                                                 slot(a,"lines")[[1]]@Lines[[l]]@coords)))
                   }))
    }


##' Internal to add legends to maps
##'
##'
##' @param legDens set to  \code{TRUE} to represent distribution of variables as the legend
##' @param legBreak Legend breaks
##' @param p Map in a ggplot format
##' @param color Color vector
##' @param limits Limits to the maps values
##' @return
##' @author Edouard Chatignoux
##' @export
##' @keywords internal
##' @examples
ggLegMap<-function(legend=list(breaks=NULL,density=F,title = NULL),p,color,limits){
    p<-p+coord_cartesian()
    legBreak<-legend$breaks
    legDens<-legend$density
    legTitle<-legend$title
    data<-p$data%>%select(id,fill)%>%unique%>%na.omit
    if (is.null(limits)) {
        breaks <- (scales:::pretty_breaks(5))(range(na.omit(data$fill)))
        breaks <- breaks[-c(1, length(breaks))]
        if (!is.null(legBreak)) breaks<-legBreak
        dt=data_frame(x=1,y=1:length(color))
        bLeg<-(breaks-min(data$fill))/diff(range(data$fill))*(length(color))+.5
        bLeg<-data_frame(x=0.5,xend=0.75,id=1)%>%
            left_join(data_frame(y=bLeg,yend=bLeg,id=1))
    }
    else {
        breaks <- (scales:::pretty_breaks(5))(limits)
        if (!legDens) breaks <- breaks[-c(1, length(breaks))]
        if (!is.null(legBreak)) breaks<-c(legBreak,limits)%>%unique%>%sort
        dt=data_frame(x=1,y=1:(length(color)))
        bLeg<-(breaks-min(limits))/diff(limits)*(length(color))+.5
        bLeg<-data_frame(x=0.5,xend=0.75,id=1)%>%
            left_join(data_frame(y=bLeg,yend=bLeg,id=1))
    }
    bLeg$breaks<-breaks
    if (!legDens){
        theme_leg <- theme_bw()+
            theme(legend.position="none",
                  axis.text.x = element_blank(), axis.line = element_blank(),
                  axis.title = element_blank(), axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  plot.background = element_rect(fill = "white",colour="black"),
                  panel.border = element_blank(),
                  axis.text.y = element_text(size = 7))

        ggLegMap<-qplot(data=dt,x,y) +
            geom_raster(aes(fill=y),interpolate=T)+
            geom_segment(data=bLeg,aes(x=x,xend=xend,y=y,yend=yend),colour=color[1])+
            geom_segment(data=bLeg,aes(x=x+0.75,xend=xend+0.75,y=y,yend=yend),colour=color[1])+
            geom_segment(data=bLeg,aes(x=x,xend=xend,y=y-0.05,yend=yend-0.05),colour=color[length(color)])+
            geom_segment(data=bLeg,aes(x=x+0.75,xend=xend+0.75,y=y-0.05,yend=yend-0.05),colour=color[length(color)])+
            scale_x_discrete(expand = c(0,0))+
            theme_leg+coord_fixed(ratio = 1)
        if (is.null(limits))
            ggLegMap<-ggLegMap+
                scale_fill_gradientn(colours = color)+
                scale_y_continuous(breaks=bLeg$y,
                                   labels=breaks,
                                   expand = c(0,0),position = "right")
        else
            ggLegMap<-ggLegMap +
                scale_fill_gradientn(colours =  color)+
                geom_rect(aes(xmin=0.5,xmax=1.5,ymin=-1,ymax=-0),fill=color[1],colour="black")+
                geom_rect(aes(xmin=0.5,xmax=1.5,ymin=length(color)+1,ymax=length(color)+2),
                          fill=color[length(color)],colour="black")+
                scale_y_continuous(breaks=c(-0.5,bLeg$y,length(color)+1.5),
                                   labels=c(parse(text = paste("phantom()<",limits[1], sep = "")),
                                            breaks,
                                            parse(text = paste("phantom()>",limits[2], sep = ""))),
                                   expand = c(0,0),position = "right")
    }
    else{
        themDens<-theme_bw()+theme(legend.position="none",
                                   axis.title.y=element_blank(),
                                   axis.text.y=element_blank(),
                                   axis.ticks.y=element_blank(),
                                   axis.title.x=element_blank(),
                                   panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.ticks = element_blank(),
                                   axis.text.x = element_text(size = 7),
                                   plot.margin=unit(c(0,0,0,0),"mm"))
        dt <- density(na.omit(data$fill),n=1000)
        dt <- data_frame(x=dt$x,y=dt$y/max(dt$y))
        dt <- data_frame(x=dt$x,xmin=dt$x-diff(dt$x)[1]/2,xmax=dt$x+diff(dt$x)[1]/2,y=dt$y)
        if (!is.null(limits)){
            if (min(dt$x)> limits[1])
                dt<-rbind(dt,
                          data_frame(x=seq(limits[1],min(dt$x),by=diff(dt$x)[1]))%>%
                                     mutate(xmin=x-diff(dt$x)[1]/2,xmax=x+diff(dt$x)[1]/2,y=0),
                          data_frame(x=limits[1],xmax=limits[1],xmin=limits[1]-min(diff(bLeg$breaks))/2,y=0))
            if (max(dt$x)< limits[2])
                dt<-rbind(dt,
                          data_frame(x=seq(max(dt$x),limits[2],by=diff(dt$x)[1]))%>%
                          mutate(xmin=x-diff(dt$x)[1]/2,xmax=x+diff(dt$x)[1]/2,y=0),
                          data_frame(x=limits[2],xmin=limits[2],xmax=limits[2]+min(diff(bLeg$breaks))/2,y=0))
            }
        ggLegMap<-ggplot(dt, aes(x, y)) +
            themDens+
            scale_y_continuous(expand = c(0, 0),limits=c(-0.05,1.05))
        if (is.null(limits))
            ggLegMap<-ggLegMap+
                geom_rect(aes(xmin=xmin,xmax=xmax,ymin=-0.05,ymax=y,fill = x))+
                scale_fill_gradientn(colours=color)+
                scale_x_continuous(breaks=bLeg$breaks,
                                   labels=bLeg$breaks,
                                   expand = c(0,0))+
                geom_hline(aes(yintercept=0),colour=color[length(color)])+
                geom_hline(aes(yintercept=-0.01),colour=color[1])+
                geom_segment(data=bLeg,aes(x=breaks,xend=breaks,y=-0.05,yend=0),linetype=5,colour=color[1])+
                geom_segment(data=bLeg,aes(x=breaks-diff(range(dt$x))/500,xend=breaks-diff(range(dt$x))/500,y=-0.05,yend=0),
                             linetype=5,colour=color[length(color)])+
                geom_line()
        else
            ggLegMap<-ggLegMap+
                geom_rect(data=dt%>%filter(x>limits[1],x<limits[2]),aes(xmin=xmin,xmax=xmax,ymin=-0.05,ymax=y,fill = x))+
                scale_fill_gradientn(colours=color)+
                geom_rect(data=dt%>%filter(x<=limits[1]),aes(xmin=xmin,xmax=xmax,ymin=-0.05,ymax=y),fill = color[1])+
                geom_rect(data=dt%>%filter(x>=limits[2]),aes(xmin=xmin,xmax=xmax,ymin=-0.05,ymax=y),fill = color[length(color)])+
                scale_x_continuous(breaks=c(limits[1]-(limits[1]-min(dt$x))/2,
                                            bLeg$breaks[-c(1,length(bLeg$breaks))],
                                            limits[2]+(max(dt$x)-limits[2])/2),
                                   labels=c(parse(text = paste("phantom()<",limits[1], sep = "")),
                                            bLeg$breaks[-c(1,length(bLeg$breaks))],
                                            parse(text = paste("phantom()>",limits[2], sep = ""))),
                                   expand = c(0,0))+
                geom_hline(aes(yintercept=0),colour=color[length(color)])+
                geom_hline(aes(yintercept=-0.01),colour=color[1])+
                geom_segment(data=bLeg,aes(x=breaks,xend=breaks,y=-0.05,yend=0),linetype=5,colour=color[1])+
                geom_segment(data=bLeg,aes(x=breaks-diff(range(dt$x))/500,xend=breaks-diff(range(dt$x))/500,y=-0.05,yend=0),linetype=5,
                             colour=color[length(color)])+
                geom_line()
    }
    if (!is.null(legTitle)) ggLegMap<-ggLegMap+ggtitle(legTitle)
    ggLegMap+theme(plot.background = element_blank())
}


##' Plot variables values on a map
##'
##'
##' \code{ggMap} create a map from a data frame with observations on spatial units.
##' It allows to represent values in continuous (with possible lower and upper bounds)
##' or categorial scale, and to hatch non significant areas.
##'
##' @param var Variable name from the data \code{data} to map, or
##' vector sorted in the order of the map spatial units
##' @param data Data with column to plot. Must contain an \code{id} variable
##' that identifies the district (jointure variable with \code{map} data)
##' @param map Map shape ("polygg" type, set to ShapDep if NULL)
##' @param limits Set limits if a continuous map is chosen
##' @param breaks Values to cut \code{var} into classes
##' @param legend List of legend parameters \describe{
##' \item{breaks}{Breaks in the legend}
##' \item{density}{Set to \code{TRUE} to represent the legend as the density distribution
##' of the observed values}
##' \item{title}{Legend title (variable name otherwise)}
##' \item{pretty}{If \code{var} is cut in classe with \code{cut} function, convert classes names to pretty labels}
##' \item{placement}{Lengend placement (x0,x1,y0,y1) coordinates,
##' in proportion of the plot zone} }
##' @param color Name of the brewer palette
##' \code{\link[RColorBrewer]{brewer.pal}}
##' @param rev Reverse color palette order (\code{TRUE,FALSE})
##' @param na.action What to do with missing values
##' @param ns Variable name in the data set that identifies non significant areas (\code{TRUE} when
##' non significant). In that case, the corresponding geographical areas are hatched.
##' @param path Draw geographical boundaries (default to \code{TRUE})
##' @return A ggplot2 graph
##' @author Edouard Chatignoux
##' @export
##' @examples
##' ## Represent crude hospitalisation rate in France by district for LP cancer in men, 2007-2011
##' library(tidyverse)
##' data(lopm.Fr)
##' crude.rate<-lopm.Fr%>%group_by(dist)%>%summarise(rate=100000*sum(H)/sum(py))%>%mutate(id=dist)
##' ggMap(var=rate,data=crude.rate,color= "YlOrRd",legend=list(title="Hosp. rate\nfor 100 000"))+
##' ggtitle("Crude hospitalisation rate for LOP, men, 2007-2011")
ggMap <-
    function (var, data = NULL, map = NULL, limits = NULL, breaks = NULL,
              legend=list(breaks=NULL,density=F,title = NULL, pretty=T,
                          placement=NULL),
              color = "Greys",rev = FALSE, na.action=na.pass,
              ns=NULL,path=TRUE){
        ## Parametres de la legende
        if (!is.null(legend))
            legend<-list(breaks=legend$breaks,
                         density=ifelse(is.null(legend$density),F,legend$density),
                         title = ifelse(is.null(legend$title),deparse(substitute(var)),legend$title),
                         pretty=ifelse(is.null(legend$pretty),T,legend$pretty),
                         placement=if(is.null(legend$placement)){
                           if (!is.null(legend$density) && legend$density==T) c(0.01,0.2,0.30,0.55)
                           else c(0,.15,0.15,.6) } else legend$placement
                         )
        ## Mef des donnees
        if (is.null(map))
            map <- ShapDep
        if (!is.null(data)) {
            if (!"id" %in% names(data))
                data$id<-data$dept
            dt <- merge(map, data, all = T)
            var <- deparse(substitute(var))
        }
        else {
            dt <- merge(data.frame(id = unique(map$id), var = var),
                        map)
            var <- "var"
        }
        dt<-eval(call(deparse(substitute(na.action)),dt))
        if (!is.null(breaks))
            dt[, var] <- cut(dt[, var], breaks)
        dt$fill <- dt[, var]
        ##ns?
        ns<-deparse(substitute(ns))
        if (ns!="NULL"){
            dt$ns<-dt[,ns]
            if (is.null(attributes(map)$hach))
              hach<-hach.poly(attributes(map)$poly)
          else hach<-attributes(map)$hach
            hach<-hach%>%left_join(dt,by=c("id"))
            }
        dt <- dt[order(dt$order), ]%>%tbl_df
        ## Theme general du graphe
        theme_map <- theme(axis.text = element_blank(), axis.line = element_blank(),
                           axis.title = element_blank(), axis.ticks = element_blank(),
                           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(),
                           plot.margin = unit(rep(0, 4), "cm"))
        ## Gestion du cas continu
        if (is.numeric(dt$fill%>%unlist)) {
            if (length(color)==1) color = RColorBrewer:::brewer.pal(9, color)
            if (rev) color = rev(color)
            ## Graphe general
            p <- ggplot(dt)+coord_quickmap()+theme_map
                p<-p+geom_polygon(aes(x = long, y = lat,
                        group = group, fill = fill))
            ## Si limites
            if (!is.null(limits)) {
                require(gridExtra)
                p <- p+scale_fill_gradientn(name = "",
                                            colours = color,
                                            limits = limits,
                                            na.value = "red")
                if (nrow(subset(dt, fill > limits[2])) > 0){
                        p <- p +
                            geom_polygon(aes(x = long, y = lat,group = group),
                                         data = subset(dt, fill >limits[2]), fill = I(color[length(color)]))
                }
                if (nrow(subset(dt, fill < limits[1])) > 0){
                        p <- p +
                            geom_polygon(aes(x = long, y = lat,group = group),data = subset(dt, fill <limits[1]), fill = I(color[1]))
                    }
            }
            ## Si pas de limites
            else
                p<-p + scale_fill_gradientn(colours = color)

            ## Gestion des depts non significatifs
            if (ns!="NULL"){
                hach<-hach%>%filter(ns==TRUE)
                p<-p+
                geom_line(data=hach,aes(x,y,group=l),color=color[1],lwd=0.1)+
                geom_line(data=hach,aes(x+1000,y,group=l),color=color[length(color)],lwd=0.1)
                }
            ## Gestion de la legende
            if (!is.null(legend)){
                g = ggplotGrob(ggLegMap(legend=legend,
                                        p=p,
                                        color=color,limits=limits))
                xR<-ggplot_build(p)$layout$panel_ranges[[1]]$x.range
                yR<-ggplot_build(p)$layout$panel_ranges[[1]]$y.range
                p<-p + theme(legend.position="none")+
                    annotation_custom(grob = g,
                                      xmin = xR[1]+(xR[2]-xR[1])*legend$placement[1],
                                      xmax = xR[1]+(xR[2]-xR[1])*legend$placement[2],
                                      ymin = yR[1]+(yR[2]-yR[1])*legend$placement[3],
                                      ymax = yR[1]+(yR[2]-yR[1])*legend$placement[4])
            }
        }
        ## Gestion du cas discret
        else {
            dt$fill <- factor(dt$fill, levels = rev(levels(dt$fill)))
            labs <- levels(dt$fill)
            if (length(color)==1) color = RColorBrewer:::brewer.pal(name = color, n = nlevels(dt$fill))
            if (rev)
                color = rev(color)
            if (legend$pretty)
                labs <- unlist(sapply(levels(dt$fill), function(x) {
                    l <- deparse(substitute(group(a, b, c),
                                            list(a = gsub("\\(","]",substr(x, 1, 1)),
                                                 b = substr(x, 2, nchar(x) -1),
                                                 c = substr(x, nchar(x), nchar(x)))))
                    if (length(grep("Inf", x)) > 0) {
                        l <- strsplit(x, ",")[[1]][1]
                        l <- paste(ifelse(length(grep("\\[", l)), "phantom()>=",
                                                             "phantom()>"), gsub("\\(", "", gsub("\\[","", l)))
                    }
                    if (length(grep("-Inf", x)) > 0) {
                        l <- strsplit(x, ",")[[1]][2]
                        l <- paste(ifelse(length(grep("]", l)), "phantom()<=",
                                          "phantom()<"), gsub(")", "", gsub("]", "",l)))
                    }
                    l
                }))
            p <- ggplot(dt) +geom_polygon(aes(x = long, y = lat,group = group, fill = fill))

            p<-p+scale_fill_manual(values = rev(color),drop = FALSE, labels = parse(text = labs), guide = guide_legend(label.hjust = 0))+
                coord_quickmap() +
                theme_map+theme(legend.position=legend$placement[c(2,4)],legend.justification=c(1,1))+
                guides(fill = guide_legend(title = legend$title))
            if (is.null(legend)) p<-p+ theme(legend.position="none")
            ## Gestion des depts non significatifs
            if (ns!="NULL"){
                hach<-hach%>%filter(id %in% (dt$id[dt$ns]%>%unique))
                p<-p+
                geom_line(data=hach,aes(x,y,group=l),color=color[1],lwd=0.1)+
                geom_line(data=hach,aes(x+1000,y,group=l),color=color[length(color)],lwd=0.1)
                }
        }
        if (path) p <- p + geom_path(aes(x = long, y = lat, group = group),
                                     size = I(1/20))
        p
    }


##' Convert a spatialpolygon to a data frame for a ggplot2 use
##'
##'
##' @param file : File path to shape file
##' @param poly : names of a SpatialPolygon if \code{NULL} file
##' @param id : id variable for  polygons
##' @return A "polygg" object, to be used with  ggplot2 (ggMap for example)
##' @author Edouard Chatignoux
##' @export
poly2ggplot <-function(file=NULL,poly=NULL,id=NULL)
{
    poly@data<-as.data.frame(poly@data)
    if (!is.null(file))  poly <- readShapeSpatial::readShapeSpatial(file)
    idn<-"ID"
    if (is.character(id) & length(id)==1)
    {
        idn<-id
        if(! idn %in% names(poly@data)) id<-NULL
    }
    if (is.null(id))
    {
        poly@data[,idn]<- sapply(slot(poly, "polygons"), function(x) slot(x, "ID"))
    }
    else
    {
        if(is.numeric(id)) poly@data[,idn]<-id
        if(is.list(id)) {idn<-names(id)
            poly@data[,idn]<-id[[1]]}
    }
    poly <- maptools::unionSpatialPolygons(poly,as.factor(poly@data[,idn]))
    map<-fortify(poly,region=idn)
    map[,idn]<-map$id
    ## On attache le polygone ? map
    attr(map,"poly")<-poly
    class(map)<-c("polygg",class(map))
    ##Carte des fronti?res
    map
}


##' Print a polygg object
##'
##'
##' @param x
##' @return Print a summary of the polygg object
##' @author Edouard Chatignoux
##' @export
##' @keywords internal
##' @examples
##' print(ShapDep)
print.polygg <- function(x) {
    cat("Map with ",length(unique(x$id))," areas\n")
}


