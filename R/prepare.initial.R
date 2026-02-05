### prepare.initial.R --- 
#----------------------------------------------------------------------
## Author: Helene
## Created: Feb  4 2026 (08:47) 
## Version: 
## Last-Updated: Feb  5 2026 (09:13) 
##           By: Helene
##     Update #: 248
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

prepare.initial <- function(dt,
                            tau = 1.2, 
                            fit.types = list(
                                Z = list(model = "Surv(tstart, tstop, delta == 2)~L0+L+A0",
                                         fit = "cox",
                                         at.risk = function(dt) (dt[["Z"]] == 0)),
                                L = list(model = "Surv(tstart, tstop, delta == 3)~L0+Z+A0",
                                         fit = "cox",
                                         at.risk = function(dt) (dt[["L"]] == 0)),
                                #W = list(model = "Surv(tstart, tstop, delta == 4)~L0+Z",
                                #         fit = "cox",
                                #         at.risk = function(dt) (dt[["W"]] == 0)),
                                outcome1 = list(model = "Surv(tstart, tstop, delta == 1)~L0+Z+L+A0",
                                                fit = "cox"),
                                #death = list(model = "Surv(tstart, tstop, delta == 5)~L0+Z+L",
                                #             fit = "cox"),
                                censoring = list(model = "Surv(tstart, tstop, delta == 0)~L0+Z+L+A0",
                                                 fit = "cox")
                            ),
                            browse = FALSE,
                            max.cap.events = 10,
                            a = NULL,
                            fit.treatment = NULL,
                            verbose = FALSE) {

    dt <- copy(dt)

    n <- length(unique(dt[["id"]]))

    #-- if no "fit" is specified for types, then specify it as cox:
    for (nm in names(fit.types)) {
        if (!is.list(fit.types[[nm]])) fit.types[[nm]] <- list(model = fit.types[[nm]], fit = "cox")
    }

    #-- if no "fit" is specified for treatment model, then specify it as glm: 
    if (!is.list(fit.treatment)) {
        fit.treatment <- list(model = fit.treatment, fit = "glm")
    }

    at.risk.ids <- (1:length(fit.types))[sapply(fit.types, function(fit.type) ("at.risk" %in% names(fit.type)))]

    #-- covariates/predictors extracted from models: 
    varnames <- unique(c(
        unlist(lapply(fit.types, function(fit.type) {
            if (length(fit.type[["model"]])>0) return(strsplit(strsplit(fit.type[["model"]], "~")[[1]][2], "\\+|\\*")[[1]])
        })),
        unlist(as.character(lapply(fit.types[at.risk.ids], function(fit.type) {
            setdiff(as.character(body(fit.type[["at.risk"]])[[2]][[2]]), c("[[", "dt"))
        })))
    ))

    #-- in case of interactions: 
    if (length(interaction.names <- grep(":", varnames, value = TRUE))>0) {
        varnames <- unique(c(varnames[!varnames %in% interaction.names],
                             str_split(interaction.names, ":")[[1]]))
    }

    #-- corresponding delta values: 
    process.deltas <- c(as.numeric(unlist(lapply(fit.types, function(fit.type) {
        tmp.delta <- strsplit(fit.type[["model"]], "~")[[1]][1]
        as.numeric(substr(tmp.delta, nchar(tmp.delta)-1, nchar(tmp.delta)-1))
    }))))

    
    #-- names of types
    process.names <- names(fit.types)
    state.names <- process.names[process.names %in% varnames[varnames %in% process.names]]
    state.deltas <- process.deltas[process.names %in% state.names]

    #-- name of baseline treatment variable:
    if (length(fit.treatment[["model"]])>0) {
        A0.name <- strsplit(fit.treatment[["model"]], "~")[[1]][1]
        varnames <- c(varnames, A0.name)
    } else {
        A0.name <- NULL 
    }

    if (typeof(dt[["time"]]) != "double") {
        warning("NB: the time variable is not numeric - will be converted")
        dt[, time := as.numeric(time)]
    }

    dt[, tstart := c(0, time[-.N]), by = "id"]
    dt[, tstop := time]
   
    for (varname in state.names) {
        dt[, (varname) := c(0, get(varname)[-.N]), by = "id"]
    }

    #--------------------------------
    #-- intervention part; for clever weight estimation:
    # (we start with cox models, if HAL is specified this is fitted later)

    cens.process.id <- (1:length(process.deltas))[process.deltas == 0]

    dt[, idN := 1:.N, by = "id"]
    if (fit.treatment[["fit"]] == "glm" & length(fit.treatment[["model"]])>0) {
        fit.A0 <- glm(as.formula(fit.treatment[["model"]]), data=dt[idN == 1], family=binomial)
        if (verbose) print(summary(fit.A0))
    } else if (length(fit.treatment[["model"]])>0)  {
        print("NB: need to incorporate other estimations methods than glm for treatment")
    }

    #--------------------------------
    #-- outcome / clever covariate part:
    # (we start with cox models, if HAL is specified this is fitted later)
    fit.cox.types <- lapply(1:length(fit.types), function(fit.type.jj) {
        if (fit.type.jj %in% at.risk.ids) {
            tmp.cox <- coxph(as.formula(fit.types[[fit.type.jj]][["model"]]),
                             data = dt[fit.types[[fit.type.jj]][["at.risk"]](dt)], 
                             control = coxph.control(timefix = FALSE))
        } else {
            tmp.cox <- coxph(as.formula(fit.types[[fit.type.jj]][["model"]]),
                             data = dt,
                             control = coxph.control(timefix = FALSE))
        }
        if (verbose) message("------------------------------")
        if (verbose) message(paste0("fit.type = ", names(fit.types)[fit.type.jj]))
        if (verbose) print(tmp.cox)
        tmp.type <- suppressWarnings(setDT(basehaz(tmp.cox, centered=TRUE)))[, (paste0("dhazard.", names(fit.types)[fit.type.jj])) := c(hazard[1],diff(hazard))][, (paste0("hazard.", names(fit.types)[fit.type.jj])) := hazard][, -"hazard", with = FALSE]
        return(list(fit.cox = tmp.cox, tmp.type = tmp.type[tmp.type[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]>0]))
    })

    #-- get all unique times; 
    unique.times <- sort(unique(dt[["time"]]))
    unique.times <- unique.times[unique.times <= tau]
    all.times <- data.table(expand.grid(time = unique.times,
                                        id = unique(dt[["id"]])))

    #-- collect data with all time-points;
    tmp.inner <- merge(dt[, time.obs := time], all.times, by = c("id", "time"), all = TRUE)[order(id, time.obs)]

    for (varname in varnames) {
        tmp.inner[, (varname) := na.locf(get(varname)), by = "id"]
    }

    #--------------------------------
    #-- expanded dataset to work with for TMLE:

    tmp.long <- tmp.inner[order(id, time)][is.na(delta), delta := 0][, -c("tstop"), with = FALSE]

    for (fit.type.jj in 1:length(fit.types)) {
        tmp.long <- merge(tmp.long, fit.cox.types[[fit.type.jj]]["tmp.type"][[1]], by = "time", all = TRUE)
    }

    tmp.long <- tmp.long[order(id,time)][!is.na(id)]

    tmp.long[, time.obs := nafill(time.obs, "nocb"), by = "id"]
    tmp.long[is.na(time.obs), time.obs := -Inf]

    tmp.long[, final.time := max(time.obs), by = "id"]

    for (fit.type.jj in 1:length(fit.types)) {
        tmp.long[is.na(get(paste0("dhazard.", names(fit.types)[fit.type.jj]))), (paste0("dhazard.", names(fit.types)[fit.type.jj])) := 0]
    }

    for (process.jj in 1:length(state.names)) {
        tmp.long[, (state.names[process.jj]) := cumsum(1*(delta == state.deltas[process.jj])), by = "id"]
    }

    for (varname in state.names) {
        tmp.long[, (varname) := c(0, get(varname)[-.N]), by = "id"]
    }

    #--------------------------------
    #-- compute needed quantities for (initial) estimation and targeting
    
    for (fit.type.jj in 1:length(fit.types)) {
        if (fit.type.jj %in% at.risk.ids) {
            tmp.long[, (paste0("at.risk.", names(fit.types)[fit.type.jj])) :=
                           fit.types[[fit.type.jj]][["at.risk"]](tmp.long)]
        }
        tmp.long[, (paste0("exp.", names(fit.types)[fit.type.jj])) := exp(predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.long, type="lp"))]
        tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj])) :=
                       tmp.long[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.long[[paste0("exp.", names(fit.types)[fit.type.jj])]]]
        if (fit.type.jj == cens.process.id) {
            tmp.long[, surv.0 := exp(-cumsum(get(paste0("P.", names(fit.types)[fit.type.jj])))), by = "id"]
            tmp.long[, surv.0.1 := c(1, surv.0[-.N]), by = "id"]
        }
    }
    
    #--------------------------------    
    #-- training part is over
    
    tmp.long <- tmp.long[time <= tau]

    #--------------------------------    
    #-- to handle dependence on jumps in the past:

    if (browse) browser()

    grid.list <- lapply(state.names,
                        function(varname) 0:min(max.cap.events, max(unique(tmp.long[[varname]]))))

    depend.matrix <- as.data.table(do.call(expand.grid, grid.list))

    names(depend.matrix) <- state.names

    depend.matrix[, state := 1:.N]

    if (length(a) == 0) {
        intervene.A0 <- FALSE
        a <- 1
    } else {
        intervene.A0 <- TRUE
        setnames(tmp.long, A0.name, "A0.obs")
        tmp.long[, pi.A0.1 := predict(fit.A0, newdata = tmp.long, type = "response")]
    }

    #--------------------------------    
    #-- initial steps to compute clever weights:
    
    tmp.long[, C := cumsum(1*(time == time.obs & delta == 0 & time %in% dt[["time"]])), by = "id"]
    tmp.long[, C.1 := c(0, C[-.N]), by = "id"]
    tmp.long[, clever.weight := (C.1 == 0)/surv.0.1]
    
    #--------------------------------    
    #-- predict for different values of state and intervention on baseline treatment:
   
    for (aa in a) {

        if (intervene.A0) {
            tmp.long[, (A0.name) := aa]
            tmp.long[, (paste0("clever.weight.a", aa)) := clever.weight*
                           (A0.obs == aa)/((pi.A0.1^(A0.obs == 1)*(1-pi.A0.1)^(A0.obs == 0)))]
            for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {
                tmp.long[, (paste0("exp.a", a, ".", names(fit.types)[fit.type.jj])) := exp(predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.long, type="lp"))]
                tmp.long[, (paste0("P.a", a, ".", names(fit.types)[fit.type.jj])) :=
                               tmp.long[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.long[[paste0("exp.a", aa, ".", names(fit.types)[fit.type.jj])]]]
            }
        }
       
        for (state.jj in depend.matrix[, unique(state)]) {

            tmp.long.jj <- copy(tmp.long)

            which.jj <- TRUE

            for (varname in state.names) {
                tmp.long.jj[[varname]] <- depend.matrix[state == state.jj][[varname]]
                which.jj <- which.jj*(tmp.long[[varname]] == tmp.long.jj[[varname]])
            }

            for (fit.type.jj in (1:length(fit.types))[-cens.process.id]) {
                tmp.long.jj[, (paste0("exp.", names(fit.types)[fit.type.jj])) := exp(predict(fit.cox.types[[fit.type.jj]]["fit.cox"][[1]], newdata=tmp.long.jj, type="lp"))]
                tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)) :=
                               tmp.long.jj[[paste0("dhazard.", names(fit.types)[fit.type.jj])]]*tmp.long.jj[[paste0("exp.", names(fit.types)[fit.type.jj])]]]
                if (fit.type.jj %in% at.risk.ids) { # <- FIX: may just want to remove again
                    tmp.long[, (paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)) :=
                                   tmp.long[[paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj)]]*
                                   fit.types[[fit.type.jj]][["at.risk"]](depend.matrix[state == state.jj])]
                }
                if (intervene.A0) {
                    setnames(tmp.long, paste0("P.", names(fit.types)[fit.type.jj], ".", state.jj),
                             paste0("P.", "a", aa, ".", names(fit.types)[fit.type.jj], ".", state.jj))
                }
            }

            tmp.long[
                which.jj == 1,
                state := state.jj]
        }
    }

    at.risks <- lapply(fit.types[at.risk.ids], function(fit.type) fit.type[["at.risk"]])
    names(at.risks) <- names(fit.types)[at.risk.ids]

    return(list(
        tmp.long = tmp.long,
        depend.matrix = depend.matrix,
        process.names = process.names,
        process.deltas = process.deltas,
        cens.process.id = cens.process.id,
        at.risks = at.risks
    ))
}

######################################################################
### prepare.initial.R ends here
