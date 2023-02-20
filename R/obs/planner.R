dd = 'C:/Users/gastauer/Documents/Projekte/HERAS Planning/'

tdef = read.csv(paste0(dd,'HERAS_Strata_2020.csv'))
ts = read.csv(paste0(dd,'HERAS2022_TransectSpacing.csv'))
tdf = ts%>%left_join(tdef)
unique(tdf$Stratum)

library(sf)
library(purrr)
library(dplyr)
library(sp)
library(tidyr)

line_stretchntrim <- function(line, polygon) {
  if (st_crs(line) != st_crs(polygon))
    return("CRS not matching")
  bb <- st_bbox(polygon)
  bbdiagLength <-
    as.numeric(sqrt((bb$xmin - bb$xmax) ^ 2 + (bb$ymin - bb$ymax) ^ 2))
  xy <- st_coordinates(line)[, 1:2]
  npairs <- nrow(xy) / 2
  etline <- NULL
  for (i in 1:npairs) {
    ii <- (i - 1) * 2 + 1
    x <- as.numeric(xy[ii:(ii + 1), 1])
    y <- as.numeric(xy[ii:(ii + 1), 2])
    dxline <- diff(x)
    dyline <- diff(y)
    d <- sqrt(dxline ^ 2 + dyline ^ 2)
    scale <- abs(as.numeric(bbdiagLength)) # * extra if need be
    signx <- sign(dxline)
    signy <- sign(dyline)
    theta <- atan(dxline / dyline)
    #  expand
    if (signy == 1) {
      dx1 <-  -sin(theta) * scale #* d
      dy1 <-  -cos(theta) * scale #* d
      dx2 <-    sin(theta) * scale #* d
      dy2 <-    cos(theta) * scale #* d
    }
    if (signy == -1) {
      dx1 <-    sin(theta) * scale# * d
      dy1 <-    cos(theta) * scale# * d
      dx2 <-  -sin(theta) * scale# * d
      dy2 <-  -cos(theta) * scale# * d
    }


    ## Cases when dxline == 0 or dyline == 0
    # dxline == 0
    if ((dxline == 0) * (signy == -1)) {
      dx1 <-  0
      dy1 <-  cos(theta) * scale# * d
      dx2 <-  0
      dy2 <-  -cos(theta) * scale# * d
    }

    if ((dxline == 0) * (signy ==  1)) {
      dx1 <-  0
      dy1 <-  -cos(theta) * scale# * d
      dx2 <-  0
      dy2 <-    cos(theta) * scale# * d
    }
    if ((signx == 1) * (dyline == 0)) {
      dx1 <-  -sin(theta) * scale# * d
      dy1 <-  0
      dx2 <-    sin(theta) * scale# * d
      dy2 <-  0
    }

    if ((signx == -1) * (dyline == 0)) {
      dx1 <-    sin(theta) * scale# * d
      dy1 <-  0
      dx2 <-  -sin(theta) * scale# * d
      dy2 <-  0
    }


    x1 <- x[1] + dx1
    y1 <- y[1] + dy1
    # second point shift
    x2 <- x[2] + dx2
    y2 <- y[2] + dy2
    # construct spatial line
    sline <- st_linestring(matrix(c(x1, y1, x2, y2),
                                  byrow = TRUE, ncol = 2))
    slineSf <- st_sf(geom = st_sfc(sline), crs = st_crs(polygon))
    # Now trim to polygon
    stline <-  st_intersection(slineSf, polygon)

    etline <- if (i == 1)
      stline
    else
      rbind(etline, stline)
  }
  etline
}



get_surv <- function(poly,spacing, Stratum, orient=0, axrot='x'){
  ##############################################
  #coordinate transformations
  spol = data.frame(x= poly$Lon, y=poly$Lat)
  coordinates(spol) <- c('x','y')
  proj4string(spol) <- suppressWarnings(CRS("+proj=longlat +datum=WGS84"))
  #get utm zone
  UTMzone <- (floor((mean(spol$x) + 180)/6) %% 60) + 1
  #make transformation
  sp_utm <- suppressWarnings(spTransform(spol, CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))))
  polygon <- st_polygon(list(as.matrix(coordinates(sp_utm))),
                        CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84")))
  ###############################################

  #rotation
  rotang = orient
  rot = function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
  tran = function(geo, ang, center) (geo - center) * rot(ang * pi / 180) + center


  inpoly <- polygon%>%sf::st_geometry()

  center <- st_centroid(st_union(inpoly))


  if(orient>0){ext = -1e5}else{ext=1e5}
  #generate grid
  if(axrot=='x'){
    center[[1]][2] = center[[1]][2] - ext
    grd <- st_make_grid(polygon%>%sf::st_geometry(),cellsize = c(1e7,as.numeric(spacing)))
    grd_rot <- tran(grd, rotang, center)

    cgrd = data.frame(st_coordinates(grd_rot))#[,1:2][!duplicated(st_coordinates(grd_rot)[,1:2]),])
    cgrd$Y = cgrd$Y + sample(seq(0,as.numeric(spacing)/2),1)
  }else{

    center[[1]][1] = center[[1]][1] + ext
    grd <- st_make_grid(polygon%>%sf::st_geometry(),cellsize = c(as.numeric(spacing),1e7))
    grd_rot <- tran(grd, rotang, center)
    # plot(polygon)
    # plot(grd, add=T)
    # plot(grd_rot, add=T)
    # points(cgrd)

    cgrd = data.frame(st_coordinates(grd_rot))
    #cgrd = data.frame(st_coordinates(grd_rot)[,1:2][!duplicated(st_coordinates(grd_rot)),])
    cgrd$X = cgrd$X + sample(seq(0,as.numeric(spacing)/2),1)
  }

  #cgrd$Ind = rep(seq(1,nrow(cgrd)/2), each=2)

  ccc = split(cgrd[,1:2],cgrd[,3])
  ccc=lapply(ccc, FUN=function(x)as.matrix(x))
  ccc=st_multilinestring(ccc)
  l2 =st_intersection(ccc, polygon)

  if(class(l2)[2]=="GEOMETRYCOLLECTION"){
    l3 = lapply(l2, FUN=function(x) if(class(x)[2]=="LINESTRING"){st_coordinates(x)})
    l3 = lapply(l3[lapply(l3,length)>0],as.matrix)
    l2=st_multilinestring(l3)
  }

  l2 = line_stretchntrim(l2, polygon)
  plot(polygon)
  plot(l2, add=T, col='red')
  #plot(ccc,add=T)
  if(class(l2[[1]])=="sfc_GEOMETRY"){
    cl=data.frame()
    for(i in 1:length(l2[[1]])){
      tmp = st_coordinates(l2[[1]][i])
      tmp=as.data.frame(tmp[,1:2])
      tmp$Transect=i
      tmp = tmp[!duplicated(tmp),]
      cl=rbind(cl, tmp)
    }
  # }
  #   cl=data.frame(do.call('rbind',lapply(1:length(l2[[1]]),
  #                             FUN= function(x){
  #                               out = st_coordinates(l2[[1]][x])
  #
  #                               if(ncol(out)>2){
  #
  #                                 return(out[, 1:3])
  #                               }else(return())})))
  }else{
    cl = data.frame(st_coordinates(l2))
    cl$Transect=rep(seq(1:(nrow(cl)/2)),each=2)
  }
  cl = cl%>%
    mutate(X=round(X,1),Y=round(Y,1))%>%
    group_by(X,Y)%>%
    summarise(X=mean(X), Y=mean(Y), Transect=min(Transect))

  coordinates(cl) = c('X','Y')
  proj4string(cl) = suppressWarnings(CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84")))
  cc=spTransform(cl, CRS(paste0("+proj=longlat +datum=WGS84")))
  out = data.frame(coordinates(cc))%>%
           mutate(Transect = cc$Transect,
                  Stratum = Stratum)
  tl = levels=seq(1, length(unique(out$Transect)))
  out$Transect = factor(out$Transect)
  levels(out$Transect)= tl
  out$Transect=as.numeric(as.character(out$Transect))
  if(axrot=='x'){
    out$sort = out$X
    out$sort[out$Transect%%2==1] = - out$X[out$Transect%%2==1]
  }else{
    out=out%>%group_by(Transect)%>%mutate(sort = row_number())
    out$sort[out$Transect%%2==0] =  -(out$sort[out$Transect%%2==0])
  }
  out=out%>%arrange(Transect,sort)
  out$sort=1:nrow(out)
  ggplot(data=out, aes(x=X,y=Y))+geom_path()+geom_text(aes(label=sort))

  plot(poly$Lon, poly$Lat, type='l', main=Stratum)
  points(out$X, out$Y)
  points(out$X, out$Y, type='l', col='red')

  return(out)
}

library(ggplot2)
library(mapdata)
xlims=c(-12,13)
ylims=c(52,62.5)
 p=ggplot(data=tdef, aes(x=Lon, y=Lat, group=Stratum, col=factor(Stratum)))+
   geom_polygon()+
  borders(database='worldHires', fill="darkgrey",colour="black",xlim=xlims,ylim=ylims)+
 coord_quickmap(xlim = xlims,ylim=ylims)

 p
#p+geom_path(data=allt, aes(x=X,y=Y, group=Stratum))

allt=data.frame()
for(i in 1:nrow(ts)){
sel=ts$Stratum[i]
if(sel==42){
  next
}
spacing = ts%>%filter(Stratum==sel)%>%select(Spacing) * 1000 %>% as.numeric() * 1.852
psel=sel
tdef$Stratum[tdef$Stratum==42] = 41
poly=tdef%>%filter(Stratum == psel)

orient=0
axrot='x'
if(sel==21){orient=-28.25}
if(sel==31){orient=-28.25}
if(sel==41){orient=-17.5; axrot='y'}
if(sel==121){orient=-32;axrot='x'}
if(sel==61){orient=0; axrot='y'}
if(sel==101){orient=-32; axrot='y'}
#if(sel==41){axrot='y'; orient=20}
tl = get_surv(poly=poly,spacing=spacing, Stratum=sel, orient=orient, axrot=axrot)
#tl$sort = tl$X
#tl$sort[tl$Transect%%2==0] = - tl$X[tl$Transect%%2==0]
#tl=tl%>%arrange(Transect,sort)
allt=rbind(allt, tl)
#p=p+geom_path(data=tl, aes(x=X, y=Y, group=Stratum),col='white')
}



allsp = lapply(unique(allt$Stratum),
               FUN=function(x)
                 SpatialLines(
                   list(
                     Lines(
                       list(
                         Line(
                           allt[allt$Stratum==x,1:2])),
                       x))))
library(leaflet)
tdf=tdf%>%filter(Stratum!=42)
pols=SpatialPolygons(
  lapply(unique(tdf$Stratum),
                     FUN=function(x)
                       Polygons(list(Polygon(tdef%>%
                                      filter(Stratum==x)%>%
                                      select('Lon','Lat'))),ID=x)),
             proj4string=CRS(paste0("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))
plot(pols)

p2 <- pols %>%
  group_by(ID) %>%
  summarise(geometry = st_union(geometry),
            AREA = sum(AREA))

popup=paste0('Stratum: ', unique(tdf$Stratum))
map = leaflet() %>%
  addTiles() %>%
  addPolygons(data=pols,popup=popup)
for(i in seq(1,length(allsp))){
  map=map%>%addPolylines(data=allsp[[i]], col='black', opacity=1, weight=1)
}

map
