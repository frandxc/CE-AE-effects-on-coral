datacache <- new.env(hash=TRUE, parent=emptyenv())

org.P_pukoensis.eg <- function() showQCData("org.P_pukoensis.eg", datacache)
org.P_pukoensis.eg_dbconn <- function() dbconn(datacache)
org.P_pukoensis.eg_dbfile <- function() dbfile(datacache)
org.P_pukoensis.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.P_pukoensis.eg_dbInfo <- function() dbInfo(datacache)

org.P_pukoensis.egORGANISM <- "Porites _pukoensis"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.P_pukoensis.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.P_pukoensis.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.P_pukoensis.eg_dbconn())
}

