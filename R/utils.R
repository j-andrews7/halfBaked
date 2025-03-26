# nocov start
.package_check <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        stop(paste(
            "Package", pkg, "is not installed, but you are",
            "attempting to use functionality provided by it. Please install it first."
        ))
    }
}
# nocov end
