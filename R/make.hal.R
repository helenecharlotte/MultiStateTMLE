### make.hal.library.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: May  1 2026 (13:35) 
## Version: 
## Last-Updated: May  4 2026 (11:19) 
##           By: Helene
##     Update #: 151
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

make.hal <- function(hal.jj) {

    if (sum(names(hal.jj) %in% c("cut.one.way", "max.count", "T.vars", "N.vars",
                                 "cut.Tk", "cut.Tk.values", "two.way"))>0) {
        return(hal.jj) 
    }

    if (sum(names(hal.jj) %in% c("components", "interactions")) == 0) {
        if (hal.jj == "default") {
            hal.jj <- list(
                components = list(
                    baseline = list(complexity = "medium"),
                    history  = list(complexity = "medium", time.mode = "auto")
                ),
                interactions = list(
                    time.with = c("history", A0.name),
                    baseline.with = "baseline"
                )
            )
        } else if (hal.jj == "no.interactions") {
            hal.jj <- list(
                components = list(
                    baseline = list(complexity = "medium"),
                    history  = list(complexity = "medium", time.mode = "auto")
                )
            )
        } else if (hal.jj == "baseline") {
            hal.jj <- list(
                components = list(
                    baseline = list(complexity = "medium"),
                    history  = list(complexity = "none", time.mode = "none")
                )
            )
        } else if (hal.jj == "sparse") {
            hal.jj <- list(
                components = list(
                    baseline = list(complexity = "medium"),
                    history  = list(complexity = "medium", time.mode = "none")
                )
            )
        } 
    }

    if (sum(names(hal.jj) %in% c("components", "interactions"))>0) {
        hal.out <- c()
        if (length(baseline <- hal.jj[["components"]][["baseline"]])>0) {
            if (baseline$complexity == "medium") {
                hal.out <- c(hal.out, cut.one.way = 20)
            } else if (baseline$complexity == "high") {
                hal.out <- c(hal.out, cut.one.way = 35)
            } else {
                hal.out <- c(hal.out, cut.one.way = 10)
            }
        }
        if (length(history <- hal.jj[["components"]][["history"]])>0) {
            if (history$complexity == "medium") {
                hal.out <- c(hal.out, max.count = 3)
            } else if (history$complexity == "high") {
                hal.out <- c(hal.out, max.count = 10)
            } else {
                hal.out <- c(hal.out, max.count = 0)
            }
            if (history$time.mode == "auto") {
                hal.out <- c(hal.out, T.vars = grep(paste0(".", 1:hal.out["max.count"], collapse = "|"),
                                                    T.vars, value = TRUE))                             
            } else if (history$time.mode == "first") {
                hal.out <- c(hal.out, T.vars = grep(".1", T.vars, value = TRUE))
            } else {
                hal.out <- c(hal.out, T.vars = "")
            }
        }
        if (FALSE) { ## not supported yet
            if (length(time <- hal.jj[["components"]][["time"]])>0) {
                if (time$complexity == "medium") {
                    hal.out <- c(hal.out, cut.time = 20)
                } else if (time$complexity == "high") {
                    hal.out <- c(hal.out, cut.time = 35)
                } else {
                    hal.out <- c(hal.out, cut.time = 10)
                }
            }
        }
        if (sum(names(hal.jj) %in% c("interactions"))>0) {

            two.way.out <- list()

            if (length(time.with <- hal.jj[["interactions"]][["time.with"]])>0) {
                if ("history" %in% time.with) {
                    for (T.var in hal.out[grep("T.vars", names(hal.out))]) {
                        two.way.out[[length(two.way.out)+1]] <-
                            c("time", T.var, 10)
                    }
                    for (N.var in setdiff(N.vars, grep("T.", N.vars, value = TRUE))) {
                        two.way.out[[length(two.way.out)+1]] <-
                            c("time", N.var, 10)
                    }
                }
                if (A0.name %in% time.with) {
                    two.way.out[[length(two.way.out)+1]] <-
                        c("time", A0.name, 10)
                }
                if ("baseline" %in% time.with) {
                    for (baseline.var in setdiff(baseline.vars, A0.name)) {
                        two.way.out[[length(two.way.out)+1]] <-
                            c("time", baseline.var, 10)
                    }
                }
            }

            if (length(baseline.with <- hal.jj[["interactions"]][["baseline.with"]])>0) {
                if ("baseline" %in% baseline.with) {
                    for (baseline.j1 in 1:(length(baseline.covars <- setdiff(baseline.vars, A0.name))-1)) {
                        for (baseline.j2 in (baseline.j1+1):length(baseline.covars)) {
                            two.way.out[[length(two.way.out)+1]] <-
                                c(baseline.covars[baseline.j1], baseline.covars[baseline.j2], 10)
                        }
                    }
                }
                if (A0.name %in% baseline.with) {
                    for (baseline.j1 in 1:(length(baseline.covars <- setdiff(baseline.vars, A0.name)))) {
                        two.way.out[[length(two.way.out)+1]] <-
                            c(baseline.covars[baseline.j1],  A0.name, 10)
                    }
                }
            }

            hal.out <- c(hal.out, two.way.out)
        }

        return(hal.out)
    } 
}


hal_spec(
    components = list(
        baseline = list(complexity = "medium"),
        history  = list(complexity = "medium", time.mode = "auto"),
        time     = list(complexity = "high")
    ),
    interactions = list(
        time.with = "history"
    )
)

add_auto_Tvars <- function(tmp.long, N.vars, process.names,
                           T.vars = character(),
                           mode = c("auto", "first", "none"),
                           q = 0.975) {
    mode <- match.arg(mode)
    if (mode == "none") return(unique(T.vars))

    proc.vars <- intersect(N.vars, process.names)

    for (N.var in proc.vars) {
        N.var.count <- tmp.long[, .(count = get(N.var)[.N]), by = "id"]
        NK <- ceiling(quantile(N.var.count$count, probs = q, na.rm = TRUE))
        if (is.na(NK) || NK < 1) next
        if (mode == "first") NK <- 1L

        T.vars <- unique(c(T.vars, paste0("T.", N.var, ".", seq_len(NK))))
    }

    unique(T.vars)
}


resolve_hal_spec <- function(spec, tmp.long, varnames, state.names, process.names) {
    components <- spec$components %||% list()
    history <- components$history %||% list()

    N.vars <- state.names
    T.vars <- grep("^T\\.", varnames, value = TRUE)

    T.vars <- add_auto_Tvars(
        tmp.long = tmp.long,
        N.vars = N.vars,
        process.names = process.names,
        T.vars = T.vars,
        mode = history$time_mode %||% "auto",
        q = history$q %||% 0.975
    )

    N.vars <- unique(c(N.vars, T.vars))
    baseline.vars <- setdiff(varnames, N.vars)

    list(
        N.vars = N.vars,
        T.vars = T.vars,
        baseline.vars = baseline.vars
    )
}

complexity_to_hal <- function(level = c("low", "medium", "high")) {
    level <- match.arg(level)
    switch(level,
           low = list(max.count = 1, cut.Tk = 1),
           medium = list(max.count = 5, cut.Tk = 5),
           high = list(max.count = 10, cut.Tk = 10)
           )
}

build_two_way_from_spec <- function(spec, groups, A0.name) {
    interactions <- spec$interactions %||% list()
    out <- list()

    time_with <- interactions$time_with %||% "none"
    if (identical(time_with, "none")) return(NULL)

    targets <- character()

    if (identical(time_with, "history")) {
        targets <- c(targets, groups$T.vars)
    } else if (identical(time_with, "baseline")) {
        targets <- c(targets, groups$baseline.vars)
    } else if (identical(time_with, "all")) {
        targets <- c(targets, groups$baseline.vars, groups$N.vars, groups$T.vars)
    }

    targets <- setdiff(unique(targets), c("", A0.name))

    for (v in targets) {
        out[[length(out) + 1]] <- c("time", v, 10)
    }

    out
}

compile_hal_spec_to_hal_sl <- function(spec, tmp.long, varnames, state.names, process.names, A0.name) {
    groups <- resolve_hal_spec(
        spec = spec,
        tmp.long = tmp.long,
        varnames = varnames,
        state.names = state.names,
        process.names = process.names
    )

    comp <- spec$components %||% list()
    baseline.level <- comp$baseline$complexity %||% "medium"
    history.level  <- comp$history$complexity %||% "medium"
    time.level     <- comp$time$complexity %||% "medium"

    base.par <- complexity_to_hal(baseline.level)
    hist.par <- complexity_to_hal(history.level)
    time.par <- complexity_to_hal(time.level)

    two.way <- build_two_way_from_spec(spec, groups, A0.name)

    # a small library for SL; you can expand this later
    hal.sl <- list(
        list(
            cut.one.way = 35,
            max.count = base.par$max.count,
            cut.Tk = base.par$cut.Tk,
            two.way = two.way,
            T.vars = groups$T.vars
        ),
        list(
            cut.one.way = 35,
            max.count = hist.par$max.count,
            cut.Tk = hist.par$cut.Tk,
            two.way = two.way,
            T.vars = groups$T.vars
        ),
        list(
            cut.one.way = 35,
            max.count = time.par$max.count,
            cut.Tk = time.par$cut.Tk,
            two.way = two.way,
            T.vars = groups$T.vars
        )
    )

    list(
        hal.sl = hal.sl,
        N.vars = groups$N.vars,
        T.vars = groups$T.vars,
        baseline.vars = groups$baseline.vars
    )
}

if (inherits(hal.sl, "hal_spec")) {
    hal.comp <- compile_hal_spec_to_hal_sl(
        spec = hal.sl,
        tmp.long = tmp.long,
        varnames = varnames,
        state.names = state.names,
        process.names = process.names,
        A0.name = A0.name
    )

    hal.sl <- hal.comp$hal.sl
    N.vars <- hal.comp$N.vars
    T.vars <- hal.comp$T.vars
    baseline.vars <- hal.comp$baseline.vars
}




######

hal_spec("default")
hal_spec("time_rich")
hal_spec("sparse")
hal_spec("no_interactions")

## also dependence on how many of the times..? all? should it be tied to the count? 

hal_spec(
    components = list(
        baseline = list(complexity = "medium"),
        history  = list(vars = c("counts", "times"), complexity = "medium"),
        time     = list(complexity = "high")
    ),
  
  interactions = list(
      time_with = c("history"),        # or "all"
      baseline_with = "none"
  )
)

hal_spec(
    baseline = c(complexity = , # if not specified, will just pick "medium"
                 interactions = "none" | "with_baseline" | "with_treatment"),  # if not specified, will just pick "none"
    history = list(c("counts", "times"),
                   complexity = ,
                   interaction = "none"  | "with_baseline" | "with_treatment"), 
    complexity = , # this is overall level, which overrides? if not specified, will just pick
    time = c(complexity = , # if not specified, will just pick
             interactions = "none" | "with_all" | "with_history" | "with_baseline")
)

hal_spec(
    baseline = "medium",
    history = "medium",
    time = list(complexity = c(),
                interactions = "history")
)

# hal_spec(
#  main.effects = c("baseline", "counts", "times"),
#  interactions = c("time_counts", "time_treatment"),
#  complexity = "medium"
# )

# structure = list(
#  main = c("baseline", "counts", "times"),
#  interactions = c("time_counts", "time_treatment")
# )

# hal.sl = make_hal_library(c(
#  "main_only",
#  "main_plus_time",
#  "rich_time_model",
#  "sparse"
# ))

# complexity = "low" | "medium" | "high"

# main.effects = c("baseline", "counts", "times")

# interactions = c("baseline", "time_times", "time_counts", "time_treatment"

# complexity =
#  "high" = list(),
#  "medium" = list(),
#  "low" = list()

# complexity.interactions =
#  "high" = list(),
#  "medium" = list(),
#  "low" = list()


make.hal.library <- function(names, data = NULL) {
    
    lapply(names, function(name) {
        switch(name,
               "main_effects" = list(
                   cut.one.way = 35,
                   max.count = 10,
                   cut.Tk = 10,
                   two.way = NULL
               ),
      
               "sparse_main" = list(
                   cut.one.way = 35,
                   max.count = 1,
                   cut.Tk = 1,
                   cut.Tk.values = 0,
                   two.way = NULL
               ),
      
               "time_interactions" = list(
                   cut.one.way = 35,
                   max.count = 10,
                   cut.Tk = 10,
                   two.way = generate.time.interactions(data)
               ),
      
               stop(paste("Unknown HAL type:", name))
               )
    })
}


generate.time.interactions <- function(data, max_degree = 10) {
    vars <- names(data)
    other_vars <- setdiff(vars, "time")
  
  lapply(other_vars, function(v) {
    c("time", v, max_degree)
  })
}

generate.interactions <- function(data, include = NULL, exclude = NULL, degree = 5) {
  vars <- names(data)
  
  if (!is.null(include)) {
    vars <- vars[vars %in% include]
  }
  if (!is.null(exclude)) {
    vars <- vars[!vars %in% exclude]
  }
  
  combn(vars, 2, simplify = FALSE, FUN = function(pair) {
    c(pair[1], pair[2], degree)
  })
}


######################################################################
### make.hal.library.R ends here
