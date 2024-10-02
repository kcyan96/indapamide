library(data.table)
library(readxl)
library(tableone)
library(survival)
library(pbapply)
library(writexl)
library(ggplot2)
library(ggpubr)

M <- function(d) year(d)*12 + month(d)
until <- function(panel, cond) panel[, shift(cummax(eval(parse(text=cond))), fill=0), id][, V1==0]
fmt <- function(x,n) format(round(x,n), nsmall=n)
merge_ <- function(L, by=NULL, all=F, ...) Reduce(function(x,y) merge(x, y, by, all=all, sort=F, ...), L)
table1 <- function(cohort, wts_col=NULL, strata="INDA") {
  cohort_ <- copy(cohort)
  cohort_[, (grep("^(dx|px|ip|rx|hx|mm|tx)",names(cohort_),value=T)) := lapply(.SD, as.logical), .SDcols=grep("^(dx|px|ip|rx|hx|mm|tx)",names(cohort_),value=T)]
  if(!is.null(wts_col)) {cohort_ <- svydesign(ids = ~ 1, data = cohort_, weights = as.formula(paste0("~ ", wts_col))); tblfun <- svyCreateTableOne}
  else {tblfun <- CreateTableOne}
  tb1_def <- setDT(read_excel("codes.xlsx", sheet="table1"))
  t1 <- tblfun(tb1_def[!is.na(Name), Name], c(strata), data=cohort_, test=F)
  print(t1, smd=T, test=F)
  return(t1)
}
as.data.table.glm <- function(m, robust=F, cluster=NULL, vcov.=NULL) {
  if(!is.null(vcov.)) m.vcov <- vcov.
  else if(!is.null(cluster)) m.vcov <- multiwayvcov::cluster.vcov(m, m$data[[cluster]])
  else if(robust) m.vcov <- sandwich::vcovHC(m, type="HC1")
  else m.vcov <- vcov(m)
  t <- lmtest::coeftest(m, vcov.=m.vcov)
  data.table(fmt(cbind(Estimate=exp(t[,"Estimate"]), exp(lmtest::coefci(m, vcov.=m.vcov)), p=t[,"Pr(>|z|)"]), 3), keep.rownames=T)
}

panel_data <- readRDS("cohort/6.panel_data.RDS")
covar_def <- setDT(read_excel("codes.xlsx", sheet="covariates"))[!is.na(Name)]
cohort <- panel_data[month==format(index.date,"%Y%m")]
cohort[, index.year := as.character(year(index.date))]
cohort[, time.since.htn := as.numeric(index.date-Date.1st.HTN)/365.25]
INCL_START <- as.Date("2005-01-01"); INCL_END <- as.Date("2020-12-31")
cohort <- cohort[index.date>=INCL_START & index.date<=INCL_END]                 
cohort <- cohort[!is.na(INDA)]                                                  
cohort <- cohort[age>=18]
cohort <- cohort[Date.1st.HTN <= index.date]
cohort <- cohort[dx.mi==0 & dx.stroke_isch==0 & dx.hf==0]
cohort <- cohort[startdte <= (index.date-365)]
cohort[, smoking := addNA(smoking)][, townsend_quintile := addNA(townsend_quintile)]
cohort[, complete := as.integer(!is.na(bmi) & !is.na(sbp) & !is.na(dbp) & !is.na(smoking) & !is.na(townsend_quintile))]
cohort <- cohort[complete==1]
tb1 <- table1(cohort)

ipw_data <- setorder(panel_data[id %in% cohort$id], id, month); stopifnot(ipw_data[, month >= format(index.date,"%Y%m")])
ipw_data[, smoking := as.integer(smoking)][, smoking := nafill(smoking, type="locf"), id][, smoking := addNA(factor(smoking, 1:3, labels=c("non-smoker","ex-smoker","current smoker")))][, townsend_quintile := addNA(townsend_quintile)]
OUTCOMES <- c("oc.mi_stroke","oc.mi","oc.stroke_isch","dead")
ipw_data[, oc.mi_stroke := as.integer(oc.mi | oc.stroke_isch)]
cohort[, dx.mi_stroke.date := pmin(dx.mi.date, dx.stroke_isch.date, na.rm=T)]
ipw_data <- merge(setnames(cohort[, c("id", "index.date", "INDA", paste0("discon.date"), paste0("switch.date"), covar_def[Baseline==1, Name]), with=F], covar_def[Baseline & Varying, Name], paste0(covar_def[Baseline & Varying, Name], ".bl")), ipw_data[, c("id", "month", "xfrd", unique(c(covar_def[Varying==1, Name], OUTCOMES))), with=F], by="id")
setorder(ipw_data, id, month); stopifnot(ipw_data[, uniqueN(month)==.N, id]$V1, ipw_data[, all(month==cummax(month)), id]$V1)
ipw_data[, deviate_any := cummax(as.integer(month >= format(pmin(get(paste0("discon.date")), get(paste0("switch.date")), na.rm=T),"%Y%m"))), id]
ipw_data[, t := M(as.Date(paste0(month,"01"), format="%Y%m%d"))-M(index.date)][, t2 := t^2]

gen_ipw <- function(pd, oc, var) {
  cat(oc, var, "...\n"); stopifnot(var %in% c("deviate_any", "xfrd"))
  ipw_spec <- list(n1 = list("INDA=='MOD'", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)]),
                   d1 = list("INDA=='MOD'", c(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], covar_def[Varying==1, Name])),
                   n0 = list("INDA=='STD'", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)]),
                   d0 = list("INDA=='STD'", c(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], covar_def[Varying==1, Name])))
  ipw_end <- list(deviate_any = paste0(oc, "|dead|xfrd|deviate_any"), xfrd = paste0(oc, "|dead|xfrd"))[[var]]
  pd[, in_ipw_model := as.integer(until(pd, ipw_end) & t!=0)]
  ipw_model <- pblapply(ipw_spec, function(s) glm(as.formula(paste0(var, " ~ t + t2 + ", paste0(s[[2]], collapse=" + "))), data = pd[in_ipw_model==1][eval(parse(text=s[[1]]))], family = binomial(link="logit"), model=F, y=F))
  pd[in_ipw_model==1 & INDA=="MOD", `:=`(n = ipw_model$n1$fitted.values, d = ipw_model$d1$fitted.values)][in_ipw_model==1 & INDA=="STD", `:=`(n = ipw_model$n0$fitted.values, d = ipw_model$d0$fitted.values)]
  pd[, w := fifelse(in_ipw_model==1, (1-n)/(1-d), 1)]  
  stopifnot(pd[, all(t==cummax(t)), id]$V1)
  pd[, ws := pmin(pmax(w, quantile(w, 0.01)), quantile(w, 0.99))][, wt := cumprod(ws), id][, wts := pmin(pmax(wt, quantile(wt, 0.01)), quantile(wt, 0.99))]
  return(pd$wts)
}
for(oc in OUTCOMES) {for(typ in c("deviate_any","xfrd")) {ipw_data[, c(paste0("wts_", typ, ".", oc)) := gen_ipw(copy(ipw_data), oc, typ)]; gc()}; rm(typ)}; rm(oc)
for(oc in OUTCOMES) ipw_data[, c(paste0("wts.", oc)) := get(paste0("wts_deviate_any.", oc)) * get(paste0("wts_xfrd.", oc))]; rm(oc)

ipw_data <- merge(ipw_data[, INDAMOD := as.integer(INDA=="MOD")], cohort[, c("id", paste0(sub("^oc.","dx.",setdiff(OUTCOMES, "dead")),".date")), with=F], by="id")
est_fmla <- paste0(" ~ INDAMOD + t + t2 + ", paste0(covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], collapse=" + "))
popu_data <- ipw_data[t==0, .SD, .SDcols = c("id", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)])]
popu_data <- setorder(rbindlist(lapply(0:119, function(n) copy(popu_data)[, t := n]))[, t2 := t^2], id, t)
res_itt <- pbsapply(OUTCOMES, function(oc) {cat("\n", oc, "...\n")
  est_data <- ipw_data[until(ipw_data, paste0(oc,"|dead|xfrd")), .SD, .SDcols = c("id", "INDAMOD", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts_xfrd.", oc))]
  evt <- est_data[, .(N=uniqueN(id), Events=sum(get(oc)), FU_median=as.numeric(median(.SD[,.N,id]$N)), FU_total=.N), keyby=-INDAMOD][, lapply(.SD, function(c) paste0(c, collapse="/"))][, INDAMOD := NULL]
  m <- glm(as.formula(paste0(oc, est_fmla)), data = est_data, weights = est_data[, get(paste0("wts_xfrd.", oc))], family = quasibinomial(link="logit"))
  res <- setDT(cbind(evt, as.data.table(m, robust=T)[rn=="INDAMOD"]))
print(res); rm(m, est_data); gc(); return(res)}, simplify=F)
res_itt <- rbindlist(res_itt, id="Outcome")[, rn := NULL]
abs_risk_itt <- pbsapply(OUTCOMES, function(oc) {cat("\n", oc, "...\n")
  est_data <- ipw_data[until(ipw_data, paste0(oc,"|dead|xfrd")), .SD, .SDcols = c("id", "INDAMOD", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts_xfrd.", oc))]
  m <- glm(as.formula(paste0(oc, paste0(est_fmla, " + I(INDAMOD * t) + I(INDAMOD * t2)"))), data = est_data, weights = est_data[, get(paste0("wts_xfrd.", oc))], family = quasibinomial(link="logit"))
  arc <- merge_(lapply(c(1,0), function(treat) setnames(cbind(popu_data[, .(id, t)], prob = predict(m, copy(popu_data)[, INDAMOD := treat], type = "response"))[, surv := cumprod(1-prob), id][, mean(1-surv), t], "V1", if(treat==1) "treated" else "control")), by="t")
rm(m, est_data); gc(); return(arc)}, simplify=F)
abs_risk_itt <- melt(rbindlist(abs_risk_itt, id="Outcome")[, `:=`(diff = treated-control, ratio = treated/control)], id.vars=c("Outcome","t"), variable.name="Type", value.name="Estimate")
res_pp <- pbsapply(OUTCOMES, function(oc) {cat("\n", oc, "...\n")
  est_data <- ipw_data[until(ipw_data, paste0(oc,"|dead|xfrd|deviate_any")), .SD, .SDcols = c("id", "INDAMOD", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts.",oc), if(startsWith(oc,"oc")) c("deviate_any", paste0(sub("^oc.","dx.",oc),".date"), paste0("discon.date"), paste0("switch.date")) else c())]
  evt <- est_data[, .(N=uniqueN(id), Events=sum(get(oc)), FU_median=as.numeric(median(.SD[,.N,id]$N)), FU_total=.N), keyby=-INDAMOD][, lapply(.SD, function(c) paste0(c, collapse="/"))][, INDAMOD := NULL]
  m <- glm(as.formula(paste0(oc, est_fmla)), data = est_data, weights = est_data[, get(paste0("wts.", oc))], family = quasibinomial(link="logit"))
  res <- setDT(cbind(evt, as.data.table(m, robust=T)[rn=="INDAMOD"]))
print(res); rm(m, est_data); gc(); return(res)}, simplify=F)
res_pp <- rbindlist(res_pp, id="Outcome")[, rn := NULL]
abs_risk_pp <- pbsapply(OUTCOMES, function(oc) {cat("\n", oc, "...\n")
  est_data <- ipw_data[until(ipw_data, paste0(oc,"|dead|xfrd|deviate_any")), .SD, .SDcols = c("id", "INDAMOD", "t", "t2", covar_def[Baseline==1, ifelse(Varying==1, paste0(Name,".bl"), Name)], oc, paste0("wts.",oc), if(startsWith(oc,"oc")) c("deviate_any", paste0(sub("^oc.","dx.",oc),".date"),  paste0("discon.date"),  paste0("switch.date")) else c())]
  m <- glm(as.formula(paste0(oc, paste0(est_fmla, " + I(INDAMOD * t) + I(INDAMOD * t2)"))), data = est_data, weights = est_data[, get(paste0("wts.", oc))], family = quasibinomial(link="logit"))
  arc <- merge_(lapply(c(1,0), function(treat) setnames(cbind(popu_data[, .(id, t)], prob = predict(m, copy(popu_data)[, INDAMOD := treat], type = "response"))[, surv := cumprod(1-prob), id][, mean(1-surv), t], "V1", if(treat==1) "treated" else "control")), by="t")
rm(m, est_data); gc(); return(arc)}, simplify=F)
abs_risk_pp <- melt(rbindlist(abs_risk_pp, id="Outcome")[, `:=`(diff = treated-control, ratio = treated/control)], id.vars=c("Outcome","t"), variable.name="Type", value.name="Estimate")

