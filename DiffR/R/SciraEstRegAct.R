#' Estimate tissue-specific transcription-factor regulatory activity
#'
#' `SciraEstRegAct` estimates transcription-factor (TF) activity from a named
#' expression vector or a gene-by-sample expression matrix. The `tissue`
#' argument automatically selects one of the epithelial regulon matrices
#' distributed with DiffR.
#'
#' @param exp.m A named numeric vector for one sample/cell, or a numeric dense
#'   or sparse matrix (or data frame) with Entrez Gene IDs in rows and
#'   samples/cells in columns. Sparse matrices must inherit from
#'   `Matrix::Matrix`.
#' @param tissue Tissue used to select the packaged regulon. Supported values
#'   are `"stomach"`, `"skin"`, `"esophagus"`, `"liver"`, `"lung"`,
#'   `"pancreas"`, `"colon"`, `"breast"`, and `"kidney"`. Matching is
#'   case-insensitive and ignores spaces and punctuation.
#' @param norm Matrix normalization method. `"z"` performs gene-wise
#'   z-scoring and `"c"` performs gene-wise mean centering. A vector input is
#'   not normalized.
#' @param ncores Number of parallel workers used for matrix input.
#'
#' @return For vector input, a named numeric vector containing one slope
#'   t-statistic per TF. For matrix input, a numeric matrix with TFs in rows
#'   and samples/cells in columns.
#'
#' @details
#' Each TF activity is the t-statistic of the regulon-predictor slope from a
#' simple linear regression with an intercept. It is not the expression value
#' of the TF itself. Expression and network genes are matched by Entrez ID.
#'
#' @references
#' Teschendorff AE and Wang N. Improved detection of tumor suppressor events
#' in single-cell RNA-Seq data. *npj Genomic Medicine* 5 (2020): 43.
#' \doi{10.1038/s41525-020-00151-y}
#'
#' @examples
#' \dontrun{
#' test_data_dir <- Sys.getenv("DIFFR_TEST_DATA_DIR")
#' load(file.path(test_data_dir, "dataStomach.rda"))
#' act.m <- SciraEstRegAct(
#'     scStomachSparse.m, tissue = "stomach", ncores = 1
#' )
#' }
#'
#' @export
SciraEstRegAct <- function(exp.m, tissue, norm = c("z", "c"), ncores = 4) {
    norm <- match.arg(norm)
    if (!is.numeric(ncores) || length(ncores) != 1L ||
        !is.finite(ncores) || ncores < 1L ||
        ncores > .Machine$integer.max || ncores %% 1 != 0) {
        stop("ncores must be a single integer greater than or equal to 1.")
    }
    ncores <- as.integer(ncores)
    if (.Platform$OS.type == "windows" && ncores > 1L) {
        warning("Parallel processing is unavailable on Windows; using ncores = 1.")
        ncores <- 1L
    }

    regnet.m <- .SciraLoadTissueNetwork(tissue)

    if (is.data.frame(exp.m)) {
        exp.m <- as.matrix(exp.m)
    }

    if (is.atomic(exp.m) && is.null(dim(exp.m))) {
        if (!is.numeric(exp.m) || is.null(names(exp.m))) {
            stop("Vector input must be numeric and named with Entrez Gene IDs.")
        }
        if (anyDuplicated(names(exp.m))) {
            stop("Vector Entrez Gene IDs must be unique.")
        }
        common.v <- intersect(names(exp.m), rownames(regnet.m))
        if (length(common.v) < 3L) {
            stop("Fewer than three genes overlap between expression and the regulon network.")
        }
        act.v <- .SciraInferTFActivity(
            as.numeric(exp.m[common.v]),
            regnet.m[common.v, , drop = FALSE]
        )
        names(act.v) <- colnames(regnet.m)
        return(act.v)
    }

    sparse.input <- inherits(exp.m, "sparseMatrix")
    dense.numeric <- is.matrix(exp.m) && is.numeric(exp.m)
    sparse.numeric <- sparse.input && inherits(exp.m, "dMatrix")
    if (!(dense.numeric || sparse.numeric)) {
        stop(paste0(
            "exp.m must be a named numeric vector, numeric dense or sparse ",
            "matrix, or numeric data frame."
        ))
    }
    if (is.null(rownames(exp.m))) {
        stop("Matrix row names must contain Entrez Gene IDs.")
    }
    if (anyDuplicated(rownames(exp.m))) {
        stop("Matrix Entrez Gene IDs must be unique.")
    }
    if (ncol(exp.m) < 1L) {
        stop("Matrix input must contain at least one sample/cell.")
    }
    if (is.null(colnames(exp.m))) {
        colnames(exp.m) <- paste0("sample_", seq_len(ncol(exp.m)))
    }

    common.v <- intersect(rownames(exp.m), rownames(regnet.m))
    if (length(common.v) < 3L) {
        stop("Fewer than three genes overlap between expression and the regulon network.")
    }
    exp.m <- exp.m[common.v, , drop = FALSE]
    regnet.m <- regnet.m[common.v, , drop = FALSE]
    if (sparse.input) {
        normalization.l <- .SciraSparseRowStats(exp.m)
        nexp.m <- exp.m
    } else {
        normalization.l <- NULL
        nexp.m <- .SciraNormalizeExpression(exp.m, norm)
    }

    idx.l <- as.list(seq_len(ncol(nexp.m)))
    out.l <- parallel::mclapply(
        idx.l,
        .SciraEstRegActPRL,
        exp.m = nexp.m,
        regnet.m = regnet.m,
        norm = norm,
        normalization.l = normalization.l,
        mc.cores = ncores
    )

    act.m <- do.call(cbind, out.l)
    rownames(act.m) <- colnames(regnet.m)
    colnames(act.m) <- colnames(nexp.m)
    return(act.m)
}


#### Auxiliary functions

.SciraTissueNetworks <- c(
    stomach = "netGAST2.m",
    skin = "netSKIN.m",
    esophagus = "netESOPH.m",
    liver = "netLIV22.m",
    lung = "netLUNGm.m",
    pancreas = "netPANC.m",
    colon = "netCOL.m",
    breast = "netBREAST.m",
    kidney = "netKID.m"
)


.SciraCanonicalTissue <- function(tissue) {
    supported <- paste(names(.SciraTissueNetworks), collapse = ", ")
    if (!is.character(tissue) || length(tissue) != 1L || is.na(tissue)) {
        stop("tissue must be one string. Supported tissues: ", supported, ".")
    }

    key <- gsub("[^[:alnum:]]", "", tolower(tissue))
    lookup <- stats::setNames(
        names(.SciraTissueNetworks),
        gsub("[^[:alnum:]]", "", names(.SciraTissueNetworks))
    )
    if (!key %in% names(lookup)) {
        stop(
            "No packaged SCIRA network is available for tissue '", tissue,
            "'. Supported tissues: ", supported, "."
        )
    }
    unname(lookup[[key]])
}


.SciraLoadTissueNetwork <- function(tissue) {
    tissue <- .SciraCanonicalTissue(tissue)
    object.name <- unname(.SciraTissueNetworks[[tissue]])
    network.env <- new.env(parent = emptyenv())

    suppressWarnings(utils::data(
        list = object.name,
        package = "DiffR",
        envir = network.env
    ))
    if (!exists(object.name, envir = network.env, inherits = FALSE)) {
        stop("Packaged SCIRA network could not be loaded: ", object.name, ".")
    }

    regnet.m <- get(object.name, envir = network.env, inherits = FALSE)
    if (!is.matrix(regnet.m) || !is.numeric(regnet.m)) {
        stop("Packaged SCIRA network is not a numeric matrix: ", object.name, ".")
    }
    if (is.null(rownames(regnet.m)) || is.null(colnames(regnet.m))) {
        stop("Packaged SCIRA network lacks gene or TF identifiers: ", object.name, ".")
    }
    if (anyDuplicated(rownames(regnet.m)) || anyDuplicated(colnames(regnet.m))) {
        stop("Packaged SCIRA network identifiers are not unique: ", object.name, ".")
    }
    return(regnet.m)
}


.SciraNormalizeExpression <- function(exp.m, norm) {
    finite.m <- exp.m
    finite.m[!is.finite(finite.m)] <- NA_real_
    rowmean.v <- rowMeans(finite.m, na.rm = TRUE)
    cexp.m <- sweep(exp.m, 1L, rowmean.v, FUN = "-")
    if (norm == "c") {
        return(cexp.m)
    }

    rowsd.v <- apply(finite.m, 1L, stats::sd, na.rm = TRUE)
    nexp.m <- cexp.m
    variable.idx <- which(is.finite(rowsd.v) & rowsd.v > 0)
    constant.idx <- which(is.finite(rowsd.v) & rowsd.v == 0)
    if (length(variable.idx) > 0L) {
        nexp.m[variable.idx, ] <- sweep(
            cexp.m[variable.idx, , drop = FALSE],
            1L,
            rowsd.v[variable.idx],
            FUN = "/"
        )
    }
    if (length(constant.idx) > 0L) {
        constant.m <- nexp.m[constant.idx, , drop = FALSE]
        constant.m[is.finite(constant.m)] <- 0
        nexp.m[constant.idx, ] <- constant.m
    }
    return(nexp.m)
}


.SciraSparseRowStats <- function(exp.m) {
    n.gene <- nrow(exp.m)
    sum.v <- numeric(n.gene)
    sumsq.v <- numeric(n.gene)
    count.v <- integer(n.gene)

    for (idx in seq_len(ncol(exp.m))) {
        exp.v <- as.numeric(exp.m[, idx])
        finite.idx <- is.finite(exp.v)
        sum.v[finite.idx] <- sum.v[finite.idx] + exp.v[finite.idx]
        sumsq.v[finite.idx] <- sumsq.v[finite.idx] + exp.v[finite.idx]^2
        count.v[finite.idx] <- count.v[finite.idx] + 1L
    }

    mean.v <- sum.v / count.v
    sd.v <- rep(NA_real_, n.gene)
    estimable.idx <- which(count.v > 1L)
    if (length(estimable.idx) > 0L) {
        variance.v <- (
            sumsq.v[estimable.idx] -
                sum.v[estimable.idx]^2 / count.v[estimable.idx]
        ) / (count.v[estimable.idx] - 1L)
        # Protect constant rows from tiny negative values caused by rounding.
        variance.v[variance.v < 0 &
            abs(variance.v) < .Machine$double.eps^0.5] <- 0
        sd.v[estimable.idx] <- sqrt(variance.v)
    }
    mean.v[count.v == 0L] <- NaN
    return(list(mean = mean.v, sd = sd.v))
}


.SciraEstRegActPRL <- function(idx, exp.m, regnet.m, norm,
                              normalization.l = NULL) {
    exp.v <- as.numeric(exp.m[, idx])
    if (!is.null(normalization.l)) {
        exp.v <- exp.v - normalization.l$mean
        if (norm == "z") {
            variable.idx <- which(
                is.finite(normalization.l$sd) & normalization.l$sd > 0
            )
            constant.idx <- which(
                is.finite(normalization.l$sd) & normalization.l$sd == 0
            )
            exp.v[variable.idx] <-
                exp.v[variable.idx] / normalization.l$sd[variable.idx]
            constant.finite.idx <- constant.idx[
                is.finite(exp.v[constant.idx])
            ]
            exp.v[constant.finite.idx] <- 0
        }
    }
    act.v <- .SciraInferTFActivity(exp.v, regnet.m)
    return(act.v)
}


.SciraInferTFActivity <- function(exp.v, regnet.m) {
    act.v <- apply(regnet.m, 2L, function(predictor.v) {
        valid.idx <- which(is.finite(exp.v) & is.finite(predictor.v))
        if (length(valid.idx) < 3L ||
            length(unique(predictor.v[valid.idx])) < 2L) {
            return(NA_real_)
        }
        fit.o <- stats::lm(exp.v[valid.idx] ~ predictor.v[valid.idx])
        coef.m <- summary(fit.o)$coefficients
        if (nrow(coef.m) < 2L) {
            return(NA_real_)
        }
        return(unname(coef.m[2L, 3L]))
    })
    return(act.v)
}
