library(shiny)
library(leaflet)
library(DT)
library(ggplot2)

source('./R/helper.R')
source('./R/plot_trans.R')
allt=NULL
getwd()
ui <- fluidPage(
  titlePanel(title=div(img(src="logo.jpg"), "Thünen Transect Planner")),
  img(src = "logo.png", height = 140, width = 400),
  fluidRow(column(3,fileInput("polfn", "Choose Polygon CSV File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"))),
  column(3,fileInput("spacefn", "Choose Spacing CSV File",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"))),
  column(3,actionButton("cmptr", "Compute Transects")),
  column(3,downloadButton("getTrans", "Download Transects"))),
  hr(),
  h3('Map'),
  br(),
  leafletOutput("survmap", width = "70%", height = 800),
  p(),
  br(),
  h3('Transect Spacing Table'),
  dataTableOutput('tsp')
)

server <- function(input, output, session) {

  values <- shiny::reactiveValues()
  if(exists("values$allt") == F){values$allt=NULL} #transects
  if(exists("values$pols") == F){values$pols=NULL} # polygon definitions
  if(exists("values$space") == F){values$space=NULL} #spacings and orientation

  ##################################################
  # POLYGONS
  ##################################################

  #get polygons / strata from csv
  getpols <- reactive({
    req(input$polfn)
    message('Adding polygons')
    read.csv(input$polfn$datapath)
  })



  #create polygons / strata
  polmap<- reactive({
    values$pols = getpols()

    tdef = values$pols

    strata=SpatialPolygons(
      lapply(unique(tdef$Stratum),
             FUN=function(x)
               Polygons(list(Polygon(tdef%>%
                                       filter(Stratum==x)%>%
                                       select('Lon','Lat'))),ID=x)),
      proj4string=CRS(paste0("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))
    message('Created Polygon projection')
    strata


  })
  #create polygon map
  output$survmap <- renderLeaflet({

    pols <- polmap()
    tdef = values$pols
    if(is.null(tdef)==F){
      popup=paste0('Stratum: ', unique(tdef$Stratum))
    }

    map = leaflet() %>%
      addMeasure(
        position = "bottomleft",
        primaryLengthUnit = "meters",
        primaryAreaUnit = "sqmeters",
        activeColor = "#3D535D",
        completedColor = "#7D4479",
        localization = "en"
        #showUnitControl= T,         # Show a control to change the units of measurements
      )%>%
      addWMSTiles(
        baseUrl = paste0(
          "https://drive.emodnet-geology.eu/geoserver/gtk/wms?"
        ),
        layers = c("gtk:seabed_substrate_multiscale"),
        options = WMSTileOptions(format = "image/png", transparent = FALSE),group="Substrate") %>%

      addWMSTiles(
        "https://ows.emodnet-bathymetry.eu/wms?",
        layers = c("emodnet:mean_atlas_land"),
        options = WMSTileOptions(format = "image/png", transparent = FALSE),
        attribution = "" ,group="EMODNET") %>%

      addProviderTiles(providers$Esri.NatGeoWorldMap, group="Nat. Geo WorldMap")%>%
      addProviderTiles(providers$OpenStreetMap, group = 'OpenStreetMap') %>%
      addProviderTiles(providers$OpenTopoMap, group='OpenTopoMap') %>%
      addProviderTiles(providers$Esri.WorldImagery, group='ESRI Imagery') %>%
      addProviderTiles(providers$Stamen.Toner, group='Stamen BW') %>%
      addLayersControl(baseGroups = c("EMODNET","Nat. Geo WorldMap",
                                      "Substrate",
                                      'OpenStreetMap',
                                      'OpenTopoMap', 'ESRI Imagery',"Stamen BW"))
      #addProviderTiles(providers$Stamen.TonerLite,
      #                 options = providerTileOptions(noWrap = TRUE)
      #)
    if(is.null(pols) == F){
      map =map%>%
        addPolygons(data=pols,popup=popup, layerId = unique(tdef$Stratum),
                    weight=1,
                    opacity = 0.8,  fillOpacity = 0.0, color="darkgray")
    }

    allt = values$allt

    if(is.null(allt)==FALSE){
      allt$Stratum_Trans=paste(allt$Stratum,allt$Transect, sep="_")
      allt$Stratum_Trans[allt$mode=="zigzag"]=allt$Stratum[allt$mode=="zigzag"]
      allt%>%group_by(Stratum_Trans)%>%summarise(n())
      allsp = lapply(unique(allt$Stratum_Trans),
                   FUN=function(x)
                     SpatialLines(
                       list(
                         Lines(
                           list(
                             Line(
                               allt[allt$Stratum_Trans==x,1:2])),
                           x))))
      for(i in seq(1,length(allsp))){
        map=map%>%addPolylines(data=allsp[[i]], col='red', opacity=1, weight=1)
      }
    }
    message('Generated Map')
    map
  })

  observeEvent(input$survmap_shape_click, { # update the location selectInput on map clicks
    p <- input$survmap_shape_click
    print(p)
  })
  observe({
    event <- input$survMap_shape_click
    if (is.null(event))
      return()
    print(event)
  })


  ###############################################################
  # SPACINGS
  ##############################################################
  #get spacings / strata from csv
  getspace <- reactive({
    req(input$spacefn)
    tmp = read.csv(input$spacefn$datapath)
    if('StratumPlan' %in% names(tmp) == F){
      tmp$StratumPlan = tmp$Stratum
    }
    if('Orientation' %in% names(tmp) ==F){
      tmp$Orientation = 0
    }
    if('Mode' %in% names(tmp) ==F){
      tmp$Mode = 'Parallel'
    }
    if(is.null(values$pols) == F){
      tmp = left_join(tmp,data.frame(Stratum=unique(values$pols$Stratum)))
      tmp$Orientation[is.na(tmp$Orientation)] = -1
    }
    values$space = tmp
    output$tsp =renderDataTable(values$space, server=FALSE,
                                rownames= FALSE,
                                editable='cell',
                                extensions = 'Buttons',
                                options = list(dom = 'Bfrtip',
                                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
    message('Adding Spacing Information')
  })

  output$tsp = reactive({
    getspace()
    #print(values$space)
    renderDataTable(values$space,
                    rownames= FALSE,
                    server=FALSE,
                    editable='cell',
                    extensions = 'Buttons',
                    options = list(dom = 'Bfrtip',
                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  })

  observeEvent(input$tsp_cell_edit, {
    message(Sys.time(), "Editing: ", input$tsp_cell_edit$col)
    values$space[input$tsp_cell_edit$row,input$tsp_cell_edit$col+1] <<- input$tsp_cell_edit$value
  })

  ###################################################################
  # COMPUTE TRANSECTS
  ###################################################################


  observeEvent(input$cmptr, {
    #get spacing
    #getspace()
    ts = values$space
    ts=ts%>%filter(ts$Spacing!=-1)
    #get polygon defs
    tdef = values$pols

    #empty df for results
    allt=data.frame()

    #loop through all polygons in spacing
    #ts=ts[1:3,]

    message('Looping through strata...')

    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Computing Transects...", value = 0)
    #ts$Stratum = as.numeric(ts$StratumPlan)


    for(stp in unique(ts$StratumPlan)){

      # if(i>1){
      #   if(as.numeric(ts$StratumPlan)[i] == as.numeric(ts$StratumPlan)[i-1]){
      #     next}
      # }
      i = which(ts$StratumPlan == stp)
      print(ts[i,])

      sel=ts$Stratum[i]
      sel2 = as.numeric(ts$Stratum[ts$StratumPlan == sel])
      if(length(sel2)<1){next}
      #tdef$Stratum[tdef$Stratum %in% sel] = stp


      ts$Spacing = as.numeric(ts$Spacing)

      spacing = ts%>%filter(Stratum==stp)%>%select(Spacing) * 1000 %>% as.numeric() * 1.852
      spacing = as.numeric(spacing)

      psel=as.numeric(sel)
      progress$inc(1/nrow(ts), detail = paste("Processing", sel, '-', i, 'of', nrow(ts)))

      message('Processing Stratum ',stp)

      poly=tdef%>%filter(Stratum %in% sel)

      orient=as.numeric(ts$Orientation[i])[1]
      if(abs(orient)>=90){
        axrot='y'
        mult = ifelse(orient < abs(orient),-1,1)
        orient = mult * (abs(orient) - 90)

      }else{
        axrot='x'
      }

      tl = get_surv(poly=poly,spacing=spacing, Stratum=stp, orient=orient, axrot=axrot, mode=ts$Mode[i[1]])
      tl = tl%>%group_by(Transect) %>% slice(c(1,n()))
      allt=rbind(allt, tl)

    }

    values$allt = allt
  })

  output$getTrans <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(values$allt, file, row.names = FALSE)
    }
  )


}
runApp(list(ui=ui,server=server), launch.browser=TRUE)

