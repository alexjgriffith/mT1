## #' main
## #'
## #' returns the pertanent information in an object
## #' @param x obj containing information to extract
## #' @param ... variables pertaining to the extraction
## #' @examples
## #' ## load(system.file("extdata","objMT1.RData",package="mT1"))
## #' sampleMT1<-mT1_sampleMT1
## #' main(1)
## #' main(sampleMT1)
## #' @rdname main
## #' @export
## main<-function(x,...){
##     UseMethod("main",x)
## }

## #' @return x unchanged
## #' @rdname main
## #' @method main default
## #' @export
## main.default<-function(x,...){
##     x
## }
