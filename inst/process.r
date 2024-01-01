#! /usr/bin/Rscript --vanilla

library(getopt)
library(mltools)

spec <- matrix(c(
##    'verbose', 'v', 2, "integer",
    'help'              , 'h', 0, "logical",    "",
    'in.name'           , 'i', 1, "character",  "",
    'out.name'          , 'o', 1, "character",  "",
    'ncores'            , 'n', 1, "integer",    "default: 96"
), byrow=TRUE, ncol=5);

opt <- getopt(spec, opt = commandArgs(TRUE));

## if help was asked for, print a friendly message and exit with a
## non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

## get rid of spurious "ARGS" list element
if (names(opt)[1] == "ARGS")
    opt <- opt[-1];

ret <- do.call(process, opt);

q();

