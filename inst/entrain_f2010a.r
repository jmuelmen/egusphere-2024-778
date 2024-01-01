#! /usr/bin/env -S "Rscript" "--vanilla"

library(getopt)
library(mltools)

spec <- matrix(c(
##    'verbose', 'v', 2, "integer",
    'help'              , 'h', 0, "logical",    "",
    'path'              , 'i', 1, "character",  "",
    'pattern'           , 'p', 1, "character",  "",
    'time.min'          , 'n', 1, "character",  "",
    'time.max'          , 'x', 1, "character",  "",
    'h.entrain'         , '9', 1, "character",  "",
    'h.tendenc'         , '8', 1, "character",  "",
    'h.aerocom'         , '7', 1, "character",  "",
    'parallel'          , 'm', 0, "logical",    ""
), byrow=TRUE, ncol=5);

opt <- getopt(spec, opt = commandArgs(TRUE));

## if help was asked for, print a friendly message and exit with a
## non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    cat("Sample invocation: salloc --nodes 1 --qos regular --time 08:00:00 --constraint cpu --account=e3sm inst/entrain_f2010a.r -i /global/cfs/cdirs/e3sm/jmuelmen/rEn0aNcDpDo2/run -p '.*.eam.h3.20[01][019].*.nc' -n 2010-01-01 -x 2012-01-01 -9 eam.h3 -8 eam.h2 -7 eam.h1 -m")
    q(status=1);
}

## get rid of spurious "ARGS" list element
if (names(opt)[1] == "ARGS")
    opt <- opt[-1];

if (is.null(opt$parallel)) { opt$parallel = FALSE; }

ret <- do.call(process.F2010A.entrain, opt);

q();

