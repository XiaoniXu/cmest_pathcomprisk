#' Internal helper: get variable names from a model call
#'
#' @export
getvarnames <- function(formula, data = NULL) {
  if (is.character(formula)) {
    return(list(varnames = formula, xvar = formula, yvar = NULL))
  }
  if (is.null(formula)) {
    return(list(varnames = NULL, xvar = NULL, yvar = NULL))
  }

  # If it is a call, try to extract formula (or just treat as formula if R allows)
  # But safer to just be robust.
  # In the original script, it was:
  formula <- formula(formula)

  lyv <- NULL
  lxv <- lvnm <- all.vars(formula[1:2])
  if (length(formula) == 3) {
    lyv <- lxv
    lxv <- all.vars(formula[-2])
    if ("." %in% lxv) {
      if (length(data) == 0) {
        stop("!getvarnames! '.' in formula and no 'data'")
      }
      lform <- formula(terms(formula, data = data))
      lxv <- all.vars(lform[-2])
    }
    lvnm <- c(lxv, lvnm)
  }
  list(varnames = lvnm, xvar = lxv, yvar = lyv)
}
