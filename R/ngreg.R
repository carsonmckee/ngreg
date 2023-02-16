rcppeigen_hello_world <- function() {
    .Call('_ngreg_rcppeigen_hello_world', PACKAGE = 'ngreg')
}

rcppeigen_outerproduct <- function(x) {
    .Call('_ngreg_rcppeigen_outerproduct', PACKAGE = 'ngreg', x)
}

rcppeigen_innerproduct <- function(x) {
    .Call('_ngreg_rcppeigen_innerproduct', PACKAGE = 'ngreg', x)
}

rcppeigen_bothproducts <- function(x) {
    .Call('_ngreg_rcppeigen_bothproducts', PACKAGE = 'ngreg', x)
}

