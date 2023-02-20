
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



get_surv <- function(poly,spacing, Stratum, orient=0, axrot='x',mode='parallel'){
  ##############################################
  #coordinate transformations
  spol = data.frame(x= poly$Lon, y=poly$Lat)
  #spol = spol[!duplicated(spol),]
  #spol = spol[order(spol$y,spol$x),]

  #plot(spol)
  #polygon(spol)
  #if(spol[1,1] != spol[length(spol),1] & spol[1,2] != spol[length(spol),2]){spol=rbind(spol, spol[,1])}
  coordinates(spol) <- c('x','y')
  proj4string(spol) <- suppressWarnings(CRS("+proj=longlat +datum=WGS84"))

  #get utm zone
  UTMzone <- (floor((mean(spol$x) + 180)/6) %% 60) + 1
  #crs="+proj=tmerc +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  #make transformation
  #sp_utm <- suppressWarnings(spTransform(spol, CRS(crs)))
  sp_utm <- suppressWarnings(spTransform(spol, CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))))

  sp_utm=rbind(sp_utm,sp_utm[1])

  polygon <- st_polygon(list(as.matrix(coordinates(sp_utm))),
                        CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84")))
  #polygon <- st_polygon(list(as.matrix(coordinates(sp_utm))),
  #                      CRS(crs))
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
    cgrd$X = cgrd$X + sample(seq(0,as.numeric(spacing)/2),1)
    # coordinates(cgrd) = c('X','Y')
    # proj4string(cgrd) = CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))
    # llgrad=spTransform(cgrd, CRS(paste0("+proj=longlat +datum=WGS84")))
    # ccl = coordinates(llgrad)
    # br = 90-geosphere::bearing(ccl[1,1:2], ccl[2,1:2])
    #
    # while(abs(br)>0.01){
    #   cgrd = data.frame(st_coordinates(grd_rot))#[,1:2][!duplicated(st_coordinates(grd_rot)[,1:2]),])
    #   coordinates(cgrd) = c('X','Y')
    #   proj4string(cgrd) = CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))
    #   llgrad=spTransform(cgrd, CRS(paste0("+proj=longlat +datum=WGS84")))
    #   ccl = coordinates(llgrad)
    #   br = 90-geosphere::bearing(ccl[1,1:2], ccl[2,1:2])-rotang
    #   print(br)
    #   grd_rot <- tran(grd_rot, br, center)
    # }
    # cgrd = data.frame(st_coordinates(grd_rot))


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
    # coordinates(cgrd) = c('X','Y')
    # proj4string(cgrd) = CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))
    # llgrad=spTransform(cgrd, CRS(paste0("+proj=longlat +datum=WGS84")))
    # ccl = coordinates(llgrad)
    # br = 90-geosphere::bearing(ccl[1,1:2], ccl[2,1:2])
    #
    # while(abs(br)>0.01){
    #   cgrd = data.frame(st_coordinates(grd_rot))
    #   coordinates(cgrd) = c('X','Y')
    #   proj4string(cgrd) = CRS(paste0("+proj=utm +zone=",UTMzone," ellps=WGS84"))
    #   llgrad=spTransform(cgrd, CRS(paste0("+proj=longlat +datum=WGS84")))
    #   ccl = coordinates(llgrad)
    #   br = 90-geosphere::bearing(ccl[1,1:2], ccl[2,1:2])-rotang
    #   plot(llgrad)
    #   grd_rot <- tran(grd_rot, br, center)
    # }
    # cgrd = data.frame(st_coordinates(grd_rot))


  }
  #cgrd$Ind = rep(seq(1,nrow(cgrd)/2), each=2)
  # tmp = cgrd%>%
  #   group_by(L2)%>%
  #   mutate(Xf=first(X[which(X==min(X))]),
  #          Yf=first(Y[which(X==min(X))]),
  #          Xl=last(X[which(X==max(X))]),
  #          Yl=last(Y[which(X==max(X))]))%>%
  #   summarise(Xf=mean(Xf),Yf=mean(Yf),Xl=mean(Xl),Yl=mean(Yl))%>%
  #   mutate(X=Xf,Y=Yf)%>%
  #   mutate(X=replace(X,L2%%2 == 1, NA),
  #          Y=replace(Y,L2%%2 == 1, NA))
  # tmp$X[is.na(tmp$X)] = tmp$Xl[is.na(tmp$X)]
  # tmp$Y[is.na(tmp$Y)] = tmp$Yl[is.na(tmp$Y)]
  # ccc=list(as.matrix(tmp[,c('X','Y')]))

  ccc = split(cgrd[,1:2],cgrd[,3])
  ccc=lapply(ccc, FUN=function(x)as.matrix(x))


  ccc=st_multilinestring(ccc)
  plot(ccc)
  l2 =st_intersection(ccc, polygon)

  if(class(l2)[2]=="GEOMETRYCOLLECTION"){
    l3 = lapply(l2, FUN=function(x) if(class(x)[2]=="LINESTRING"){st_coordinates(x)})
    l3 = lapply(l3[lapply(l3,length)>0],as.matrix)
    l2=st_multilinestring(l3)
  }

  l2 = line_stretchntrim(l2, polygon)
  plot(polygon)
  plot(ccc,add=T)
  plot(l2, add=T, col='red')
  #plot(ccc,add=T)
  if(class(l2[[1]])[1]=="sfc_GEOMETRY"){
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

  if(rotang==0){
    if(axrot=='x'){
      out=out[order(out$Transect),]
      for(i in seq(2, nrow(out), by=2)){
        out$Y[i]=out$Y[i-1]
      }
    }else{
      for(i in seq(2, nrow(out), by=2)){
        out$X[i]=out$X[i-1]
      }
    }
  }
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
  #ggplot(data=out, aes(x=X,y=Y))+geom_path()+geom_text(aes(label=sort))

  #plot(poly$Lon, poly$Lat, type='l', main=Stratum)
  #points(out$X, out$Y)
  #points(out$X[seq(1,nrow(out),2)], out$Y[seq(1,nrow(out),2)], type='l', col='red')
  if(mode=='zigzag'){
    out=out[seq(1, nrow(out),2),]
  }

  return(out)
}

