
* ==============================================================================  
* Title: Impact of national health insurance on health services utilization in Lao People's Democratic Republic (Lao PDR): a quasi-experimental evaluation using longitudinal administrative data
* Author: Sekeun Daniel Yu (yus109@mcmaster.ca), Michel Grignon, Godefroy Emmanuel Guindon, Jean-Éric Tarride
* Date: July/2025
* ==============================================================================  

clear all
global dir "/enter-directory-here/"
cd "$dir"

capture mkdir "$dir/stata"
capture mkdir "$dir/stata/output"
capture mkdir "$dir/stata/graph"
capture mkdir "$dir/stata/graph/temp"
capture mkdir "$dir/stata/temp"

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Count total number of districts, excluding Provincial Hospitals (PH)
preserve
keep if faclevel != "PH"
egen district_tag = tag(district)
count if district_tag	// total district is 148
restore

* Count total number of districts with NHI (treat == 1), excluding PH
preserve
keep if faclevel != "PH" & treat == 1
egen district_tag = tag(district)
count if district_tag	// NHI district is 138
restore

* Number of districts per province, excluding PH
preserve
keep if faclevel != "PH"
duplicates drop province district, force
bysort province: gen count = _N
bysort province (district): keep if _n == 1
list province count, clean
restore

* Number of observations per province and district
preserve
contract province district
list province district _freq, sepby(province)
restore

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"
//drop if inlist(quarter, 2015.1, 2015.2, 2015.3)

* Districts in final analysis (treat == 1 and excluding PH)
preserve
keep if faclevel != "PH" & treat == 1
egen district_tag = tag(district)
count if district_tag	// Final analysis districts is 129
restore

* Number of health facilities at event time == 0 (by faclevel)
preserve
keep if faclevel != "PH" & treat == 1 & event == 0
egen tag = tag(district)
egen district_count = total(tag), by(faclevel)
egen total_facnum = total(facnum), by(faclevel)
bysort faclevel: keep if _n == 1
list faclevel district_count total_facnum
restore

/*
# 129 districts
# (116 district with district hospitals + 13 districts without district hospital)
# (945-116+91 = 920 health centers)
*/


/*
* ==============================================================================  
* Set-up
* ==============================================================================  

# Provinces and the month of inception
# -----------------------------------
# VTC No intervention (never-treated control group);
# PSL 2017/02/01; LNT 2016/08/01; ODX 2017/02/01; BKO 2017/10/01; LPG 2017/10/01;
# HPN 2017/01/01; XYB 2017/10/01; XKG 2017/01/01; VTP 2017/10/01; BKX 2017/01/01;
# KMN 2017/10/01; SVK 2017/09/01; SRV 2016/12/01; SEK 2017/06/01; CPS 2017/09/01;
# ATP 2016/07/01; XSB 2016/08/01

# Intervention sequence:
# -----------------------------------
# 1st treated (Q3/2016): 17 ATP, 03 LNT, 18 XSB
# 2nd treated (Q4/2016): 14 SRV
# 3rd treated (Q1/2017): 11 BKX, 07 HPN, 09 XKG, 04 ODX, 02 PSL
# 4th treated (Q3/2017): 13 SVK, 16 CPS
# 5th treated (Q4/2017): 05 BKO, 12 KMN, 06 LPG, 08 XYB, 10 VTP
*/

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* ==============================================================================  
* Figure 3. Difference-in-differences estimates by utilization outcome
* ==============================================================================  

/*
Group-Time ATT (gt-ATT) estimator (Sant'Anna and Callaway 2021)
Dynamic treatment effect with simultaneous confidence intervals (sci)
Intervention group: Q3/2016, Q4/2016, Q1/2017, Q3/2017
Control group: not-yet-treated (Q4/2016, Q1/2017, Q3/2017) + last-treated (Q4/2017)
*/

ssc install drdid, all replace
ssc install csdid, all replace
// help csdid

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
    csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
        method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) ///
        rseed(555) saverif(_csdidrif_) replace // sci applied
	estimates store tbl_fig3_`m'
	estat pretrend			// sci not applied

	csdid_plot, title("`t'", size(small)) legend(off) ///
    xlabel(-7(1)4, labsize(vsmall)) ///
    ylabel(, labsize(vsmall)) ///
    xtitle("Time since NHI introduction", size(vsmall)) ///
    ytitle("ATT", size(vsmall))
    graph save "stata/graph/temp/fig_03_did_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_fig3_opo5r tbl_fig3_opu5r tbl_fig3_ipo5r ///
	   tbl_fig3_ipu5r tbl_fig3_iplos tbl_fig3_delr ///
using "stata/output/fig_03_did.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/fig_03_did_opo5r.gph" ///
"stata/graph/temp/fig_03_did_opu5r.gph" ///
"stata/graph/temp/fig_03_did_ipo5r.gph" ///
"stata/graph/temp/fig_03_did_ipu5r.gph" ///
"stata/graph/temp/fig_03_did_iplos.gph" ///
"stata/graph/temp/fig_03_did_delr.gph" ///
, rows(3) cols(2) iscale(1) graphregion(margin(0 0 0 0)) xsize(8) ysize(10)
graph export "stata/graph/fig_03_did_main.pdf", replace

/*
Simultaneous confidence intervals differ due to the use of multiplier bootstrap with 1,000 iterations. Exact replication is not possible because Stata and R use different random number generation processes.

preserve
use _csdidrif_, clear
csdid_stats event, wboot(reps(1000)) rseed(555)	// sci applied
csdid_plot, title("DiD dynamic treatment effect")
restore
*/


* ==============================================================================  
* Appendix Figure A2: DiD analysis by facility type
* ==============================================================================

* ----------------------------------------------------------
* District hospitals
* ----------------------------------------------------------

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
keep if faclevel == "DH"
drop if district == "1210 Khounkham"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig5_dh_`m'
	estat pretrend			// sci not applied

	csdid_plot, title("`t' (DH)", size(small)) legend(off) ///
    xlabel(-7(1)4, labsize(vsmall)) ///
    ylabel(, labsize(vsmall)) ///
    xtitle("Time since NHI introduction", size(vsmall)) ///
    ytitle("ATT", size(vsmall))
    graph save "stata/graph/temp/sfig_05_did_facility_dh_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig5_dh_opo5r tbl_sfig5_dh_opu5r tbl_sfig5_dh_ipo5r ///
	   tbl_sfig5_dh_ipu5r tbl_sfig5_dh_iplos tbl_sfig5_dh_delr ///
using "stata/output/sfig_05_did_facility_dh.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* Health centers
* ----------------------------------------------------------

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
keep if faclevel == "HC"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig5_hc_`m'

	csdid_plot, title("`t' (HC)", size(small)) legend(off) ///
    xlabel(-7(1)4, labsize(vsmall)) ///
    ylabel(, labsize(vsmall)) ///
    xtitle("Time since NHI introduction", size(vsmall)) ///
    ytitle("ATT", size(vsmall))
    graph save "stata/graph/temp/sfig_05_did_facility_hc_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig5_hc_opo5r tbl_sfig5_hc_opu5r tbl_sfig5_hc_ipo5r ///
	   tbl_sfig5_hc_ipu5r tbl_sfig5_hc_iplos tbl_sfig5_hc_delr ///
using "stata/output/sfig_05_did_facility_hc.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_05_did_facility_dh_opo5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_opo5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_dh_opu5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_opu5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_dh_ipo5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_ipo5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_dh_ipu5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_ipu5r.gph" ///
"stata/graph/temp/sfig_05_did_facility_dh_iplos.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_iplos.gph" ///
"stata/graph/temp/sfig_05_did_facility_dh_delr.gph" ///
"stata/graph/temp/sfig_05_did_facility_hc_delr.gph" ///
, rows(6) cols(2) iscale(1) graphregion(margin(0 0 0 0)) xsize(8) ysize(22)
graph export "stata/graph/sfig_05_did_facility.pdf", replace


* ==============================================================================  
* Appendix Figure A3: DiD analysis by district-level poverty status
* ==============================================================================  

ssc install event_plot, replace

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0)

/*
'xthdidregress' command produces results, but according to the output in R, the estimates are only partially identifiable for certain cohorts due to a violation of the overlap condition in the IPW component ("warning: fitted probabilities numerically 0 or 1 occurred"), indicating that the estimates are not doubly robust and may be biased in some cases. This warning appears to be ignored or surpressed in Stata.

Based on the results in R, possible alternatives include estimating an outcome model with/without ovariates, or applying different DiD estimators that do not rely on inverse probability weighting.

In the manuscript, I employed the interaction-weighted estimator of Sun and Abraham (2021) with covariates and the last-treated group as the control group.
*/

* ----------------------------------------------------------
* Interaction-weighted DiD estimator (Sun and Abraham 2021)
* ----------------------------------------------------------

/*
Install user-written package:

net install github, from("https://haghish.github.io/github/")
github install lsun20/eventstudyinteract
github update eventstudyinteract

Resource: https://www.sciencedirect.com/science/article/abs/pii/S030440762030378X
		  https://github.com/lsun20/EventStudyInteract
*/

gen group_na = group
recode group_na (9=.)

gen rt = time - group_na
gen never_treat = (group == 9)

tab rt

* Create dummy
forvalues k = 7(-1)1 {
	gen lead_`k' = rt == -`k'
	}
forvalues k = 0/4 {
	gen lag_`k' = rt == `k'
	}

* ----------------------------------------------------------
* High-poverty area *
* ----------------------------------------------------------

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 1

	eventstudyinteract `m' lead_7-lead_2 lag_0-lag_4, cohort(group_na) control_cohort(never_treat) absorb(i.dnum i.time) vce(cluster pnum) covariates(facnum semp urba litr imps)

	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Time since NHI introduction", size(vsmall)) ytitle("ATT", size(vsmall)) xlabel(-7(1)4, labsize(vsmall)) title("`t' (high-poverty)", size(small))) stub_lead(lead_#) stub_lag(lag_#) trimlag(4) trimlead(7) together
	graph save "stata/graph/temp/sfig_06_did_poverty_high_`m'.gph", replace
	
	matrix b = e(b_iw)
	matrix V = e(V_iw)
	ereturn post b V
	esttab
	estimates store tbl_sfig6_high_`m'
	
	restore
}

* Export NHI estimates
esttab tbl_sfig6_high_opo5r tbl_sfig6_high_opu5r tbl_sfig6_high_ipo5r ///
	   tbl_sfig6_high_ipu5r tbl_sfig6_high_iplos tbl_sfig6_high_delr ///
using "stata/output/sfig_06_did_poverty_high.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* Low-poverty area *
* ----------------------------------------------------------

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 0

	eventstudyinteract `m' lead_7-lead_2 lag_0-lag_4, cohort(group_na) control_cohort(never_treat) absorb(i.dnum i.time) vce(cluster pnum) covariates(facnum semp urba litr imps)

	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Time since NHI introduction", size(vsmall)) ytitle("ATT", size(vsmall)) xlabel(-7(1)4, labsize(vsmall)) title("`t' (low-poverty)", size(small))) stub_lead(lead_#) stub_lag(lag_#) trimlag(4) trimlead(7) together
	graph save "stata/graph/temp/sfig_06_did_poverty_low_`m'.gph", replace
	
	matrix b = e(b_iw)
	matrix V = e(V_iw)
	ereturn post b V
	esttab
	estimates store tbl_sfig6_low_`m'
	
	restore
}

* Export NHI estimates
esttab tbl_sfig6_low_opo5r tbl_sfig6_low_opu5r tbl_sfig6_low_ipo5r ///
	   tbl_sfig6_low_ipu5r tbl_sfig6_low_iplos tbl_sfig6_low_delr ///
using "stata/output/sfig_06_did_poverty_low.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_06_did_poverty_high_opo5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_opo5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_high_opu5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_opu5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_high_ipo5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_ipo5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_high_ipu5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_ipu5r.gph" ///
"stata/graph/temp/sfig_06_did_poverty_high_iplos.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_iplos.gph" ///
"stata/graph/temp/sfig_06_did_poverty_high_delr.gph" ///
"stata/graph/temp/sfig_06_did_poverty_low_delr.gph" ///
, rows(6) cols(2) iscale(1) graphregion(margin(0 0 0 0)) xsize(8) ysize(22)
graph export "stata/graph/sfig_06_did_poverty.pdf", replace


* ==============================================================================  
* Appendix Figure A4: DiD analysis by intervention group
* ==============================================================================  

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig7_`m'
	
	csdid_plot, group(4) title("`t' (group 4)", size(vsmall)) legend(off) ///
        xlabel(-2(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_07_did_group_4_`m'.gph", replace
	
	csdid_plot, group(5) title("`t' (group 5)", size(vsmall)) legend(off) ///
        xlabel(-3(1)3, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny)) 
    graph save "stata/graph/temp/sfig_07_did_group_5_`m'.gph", replace
	
	csdid_plot, group(6) title("`t' (group 6)", size(vsmall)) legend(off) ///
        xlabel(-4(1)2, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny)) 
    graph save "stata/graph/temp/sfig_07_did_group_6_`m'.gph", replace

	csdid_plot, group(8) title("`t' (group 8)", size(vsmall)) legend(off) ///
        xlabel(-6(1)0, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny)) 
    graph save "stata/graph/temp/sfig_07_did_group_8_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_sfig7_opo5r tbl_sfig7_opu5r tbl_sfig7_ipo5r ///
	   tbl_sfig7_ipu5r tbl_sfig7_iplos tbl_sfig7_delr ///
using "stata/output/sfig_07_did_group.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_07_did_group_4_opo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_5_opo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_6_opo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_8_opo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_4_opu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_5_opu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_6_opu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_8_opu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_4_ipo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_5_ipo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_6_ipo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_8_ipo5r.gph" ///
"stata/graph/temp/sfig_07_did_group_4_ipu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_5_ipu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_6_ipu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_8_ipu5r.gph" ///
"stata/graph/temp/sfig_07_did_group_4_iplos.gph" ///
"stata/graph/temp/sfig_07_did_group_5_iplos.gph" ///
"stata/graph/temp/sfig_07_did_group_6_iplos.gph" ///
"stata/graph/temp/sfig_07_did_group_8_iplos.gph" ///
"stata/graph/temp/sfig_07_did_group_4_delr.gph" ///
"stata/graph/temp/sfig_07_did_group_5_delr.gph" ///
"stata/graph/temp/sfig_07_did_group_6_delr.gph" ///
"stata/graph/temp/sfig_07_did_group_8_delr.gph" ///
, rows(6) cols(4) iscale(1) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(37) ysize(60)
graph export "stata/graph/sfig_07_did_group.pdf", replace


* ==============================================================================  
* Appendix Figure A5: DiD analysis with adjusted model specifications
* ==============================================================================  

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

* ----------------------------------------------------------
* Base model
* ----------------------------------------------------------

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
        method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) ///
        rseed(555) saverif(_csdidrif_) replace // sci applied
	estimates store tbl_sfig8_base_`m'

	csdid_plot, title("`t' (base)", size(vsmall)) legend(off) ///
		xlabel(-7(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_08_did_adjust_base_`m'.gph", replace

* ----------------------------------------------------------
* Balanced NHI exposure (e=0-2)
* ----------------------------------------------------------

* Equation 3.11 from Callaway and Sant'Anna (2021) is implementable in R but not available through 'xthdidregress' or 'csdid'.

* ----------------------------------------------------------
* Last-treated units as control
* ----------------------------------------------------------

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig8_last_`m'

	csdid_plot, title("`t' (last)", size(vsmall)) legend(off) ///
		xlabel(-7(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_08_did_adjust_last_`m'.gph", replace

* ----------------------------------------------------------
* No covariate (unconditional Parallel Trend Assumption)
* ----------------------------------------------------------

csdid `m', ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig8_nocov_`m'

	csdid_plot, title("`t' (no covarite)", size(vsmall)) legend(off) ///
		xlabel(-7(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_08_did_adjust_nocov_`m'.gph", replace

}

* Export NHI estimates
esttab tbl_sfig8_base_opo5r tbl_sfig8_base_opu5r tbl_sfig8_base_ipo5r ///
	   tbl_sfig8_base_ipu5r tbl_sfig8_base_iplos tbl_sfig8_base_delr ///
	   using "stata/output/sfig_08_did_adjust_base.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

esttab tbl_sfig8_last_opo5r tbl_sfig8_last_opu5r tbl_sfig8_last_ipo5r ///
	   tbl_sfig8_last_ipu5r tbl_sfig8_last_iplos tbl_sfig8_last_delr ///
	   using "stata/output/sfig_08_did_adjust_last.xlsx", ///
	   mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

esttab tbl_sfig8_nocov_opo5r tbl_sfig8_nocov_opu5r tbl_sfig8_nocov_ipo5r ///
	   tbl_sfig8_nocov_ipu5r tbl_sfig8_nocov_iplos tbl_sfig8_nocov_delr ///
	   using "stata/output/sfig_08_did_adjust_nocov.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_08_did_adjust_base_opo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_opo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_opo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_base_opu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_opu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_opu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_base_ipo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_ipo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_ipo5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_base_ipu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_ipu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_ipu5r.gph" ///
"stata/graph/temp/sfig_08_did_adjust_base_iplos.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_iplos.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_iplos.gph" ///
"stata/graph/temp/sfig_08_did_adjust_base_delr.gph" ///
"stata/graph/temp/sfig_08_did_adjust_last_delr.gph" ///
"stata/graph/temp/sfig_08_did_adjust_nocov_delr.gph" ///
, rows(6) cols(3) iscale(1) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) ///
 xsize(35) ysize(65)
graph export "stata/graph/sfig_08_did_adjust.pdf", replace


* ==============================================================================  
* Appendix Figure A6: DiD analysis using alternative heterogeneity-robust estimators
* ==============================================================================  

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* ---------------------------------------------------------
* Extended TWFE DID (Wooldrigde 2021)
* ----------------------------------------------------------

/*
Install user-written package:
ssc install jwdid, all replace

Resource: https://doi.org/10.1093/ectj/utad016
		  https://friosavila.github.io/app_metrics/app_metrics11.html
*/

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

	jwdid `m' facnum_fixed semp urba litr imps, ivar(dnum) tvar(time) gvar(group_0) cluster(pnum)
	estat event, estore(tbl_sfig9_jw_`m')
	
	estat plot, title("`t' (Extended TWFE)", size(vsmall)) legend(off) ///
		xlabel(-7(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
	graph save "stata/graph/temp/sfig_09_did_hetero_jwdid_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_sfig9_jw_opo5r tbl_sfig9_jw_opu5r tbl_sfig9_jw_ipo5r ///
	   tbl_sfig9_jw_ipu5r tbl_sfig9_jw_iplos tbl_sfig9_jw_delr ///
using "stata/output/sfig_09_did_hetero_jwdid.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

/*
The estimates differ between R and Stata, and also between the 'jwdid' and 'xthdidregress' commands within Stata. These discrepancies are likely due to differences in how the commands handle covariates that are collinear with fixed effects. When covariates are excluded, the estimates from R and Stata are identical.
*/

* ----------------------------------------------------------
* Stacked DiD (Wing 2024)
* ----------------------------------------------------------

/*
Code is from:
https://rawcdn.githack.com/hollina/stacked-did-weights/18a5e1155506cbd754b78f9cef549ac96aef888b/stacked-example-r-and-stata.html

Stacked DiD uses a stack of trimmed, cohort-specific panels to estimate treatment effects. Stacked DiD estimator uses trimmed aggregate ATT which is a weighted average of group-time ATTs from a trimmed set of data. This is analogous to the balanced NHI exposure configuration in Equation 3.11 of Callaway and Sant'Anna (2021).

The estimates and confidence intervals are the same in R and Stata.
*/

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if province != "01 Vientiane Capital"

gen group_na = group
recode group_na (9=.)
rename treat treat_org


/* Create sub-experiment data for stack */

* clear programs
capture program drop _all

* start new program
program create_sub_exp
	syntax, ///
		timeID(string) ///
        groupID(string) ///
        adoptionTime(string) ///
        focalAdoptionTime(int) ///
        kappa_pre(numlist) ///
        kappa_post(numlist)
    
    * Suppress output
    qui {
        * Save dataset in memory, so we can call this function multiple times. 
        preserve

        * Determine earliest and latest time in the data. 
        * Used for feasibility check later
        sum `timeID'
        local minTime = r(min)
        local maxTime = r(max)

        * variable to label sub-experiment if treated in focalAdoptionTime
        gen sub_exp = `focalAdoptionTime' if `adoptionTime' == `focalAdoptionTime'

        * Now fill in this variable for states with adoptionTime > focalAdoptionTime + kappa_post
        * note, this will include never treated, because adopt_year is ., which stata counts as infinity
        replace sub_exp = `focalAdoptionTime' if `adoptionTime' > `focalAdoptionTime' + `kappa_post'

        * Keep only treated and clean controls
        keep if sub_exp != .

        * gen treat variable in subexperiment
        gen treat = `adoptionTime' == `focalAdoptionTime'

        * gen event_time 
        gen event_time = time - sub_exp

        * gen post variable
        gen post = event_time >= 0

        * trim based on kappa's: -kappa_pre < event_time < kappa_post
        keep if inrange(event_time, -`kappa_pre', `kappa_post')

        * keep if event_time >= -`kappa_pre' & event_time <= `kappa_post'
        gen feasible = 0 
        replace feasible = 1 if !missing(`adoptionTime')
        replace feasible = 0 if `adoptionTime' < `minTime' + `kappa_pre' 
        replace feasible = 0 if `adoptionTime' > `maxTime' - `kappa_post' 
        drop if `adoptionTime' < `minTime' + `kappa_pre'

        * Save dataset
        compress
        save stata/temp/subexp`focalAdoptionTime', replace
        restore
    }
end

// Build the stack of sub-experiments //

* create the sub-experimental data sets
levelsof group_na, local(alist)
di "`alist'"

* Loop over the events and make a data set for each one
foreach j of local alist { 
    * Preserve dataset
    preserve

    * run function
    create_sub_exp, ///
        timeID(time) ///
        groupID(dnum) ///
        adoptionTime(group_na) ///
        focalAdoptionTime(`j') ///
        kappa_pre(3) ///
        kappa_post(2)

    * restore dataset
    restore
}

* Append the stacks together, but only from feasible stacks
* Determine earliest and latest time in the data. 
* Used for feasibility check later
sum time
local minTime = r(min)
local maxTime = r(max)
local kappa_pre = 3
local kappa_post = 2

gen feasible_time = group_na
replace feasible_time = . if group_na < `minTime' + `kappa_pre'
replace feasible_time = . if group_na > `maxTime' - `kappa_post'
sum feasible_time

local minadopt = r(min)
levelsof feasible_time, local(alist)
clear

foreach j of local alist {
    display `j'
    if `j' == `minadopt' use stata/temp/subexp`j', clear
    else append using stata/temp/subexp`j'
}

* Clean up
* Group 9 is not-yet-treated, not never-treated; thus, sub-experiment 8 must be dropped as it doesn't have a clean control.
erase stata/temp/subexp8.dta
drop if sub_exp == 8

* Summarize
sum dnum time group_na opo5r treat post event_time feasible sub_exp

* Treated, control, and total count by stack
preserve
keep if event_time == 0
gen N_treated = treat
gen N_control = 1 - treat
gen N_total = 1
collapse (sum) N_treated N_control N_total, by(sub_exp)
list sub_exp N_treated N_control N_total in 1/3
restore

// The `compute_weights()` function //

capture program drop _all

program compute_weights
    syntax, ///
        treatedVar(string) ///
        eventTimeVar(string) ///
        groupID(string) ///
        subexpVar(string)

    * Create weights
    bysort `subexpVar' `groupID': gen counter_treat = _n if `treatedVar' == 1
    egen n_treat_tot = total(counter_treat)
    by `subexpVar': egen n_treat_sub = total(counter_treat)

    bysort `subexpVar' `groupID': gen counter_control = _n if `treatedVar' == 0
    egen n_control_tot = total(counter_control)
    by `subexpVar': egen n_control_sub = total(counter_control)

    gen stack_weight = 1 if `treatedVar' == 1
    replace stack_weight = (n_treat_sub/n_treat_tot)/(n_control_sub/n_control_tot) if `treatedVar' == 0
end

* Compute the stacked weights
compute_weights, ///
    treatedVar(treat) ///
    eventTimeVar(event_time) ///
    groupID(dnum) ///
    subexpVar(sub_exp)

* Summarize 
sumup stack_weight if treat == 0 & event_time == 0, by(sub_exp) s(mean)

// Estimate the stacked regression //

* Create dummy variables for event-time
char event_time[omit] -1
xi i.event_time

* Rename
rename _Ievent_tim_1 ievent_lead_3
rename _Ievent_tim_2 ievent_lead_2
rename _Ievent_tim_4 ievent_lag_0
rename _Ievent_tim_5 ievent_lag_1
rename _Ievent_tim_6 ievent_lag_2

foreach i in 2 3 {
    gen lead_`i' = treat * (ievent_lead_`i' == 1)
}
foreach i in 0 1 2 {
    gen lag_`i' = treat * (ievent_lag_`i' == 1)
}


* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

	reghdfe `m' facnum_fixed semp urba litr imps lead_* lag_* [aw = stack_weight], cluster(pnum) absorb(treat event_time)
	estimates store tbl_sfig9_stk_`m'

	* Show results
	esttab tbl_sfig9_stk_`m', keep(lead* lag*) se

	* Compute the average post-treatment effect
	lincom (lag_0 + lag_1 + lag_2)/3

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' (Stacked DiD)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/sfig_09_did_hetero_stkdid_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_sfig9_stk_opo5r tbl_sfig9_stk_opu5r tbl_sfig9_stk_ipo5r ///
	   tbl_sfig9_stk_ipu5r tbl_sfig9_stk_iplos tbl_sfig9_stk_delr ///
using "stata/output/sfig_09_did_hetero_stkdid.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* Dynamic TWFE DiD
* ----------------------------------------------------------

/*
This estimator is possibly biased in the context of treatment effect heterogeneity. Code is from:

https://lost-stats.github.io/Model_Estimation/Research_Design/event_study.html#:~:text=Difference-in-Differences%20Event%20Study,periods%20in%20your%20respective%20study.

The estimates and confidence intervals are the same in R and Stata.
*/

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"

xtset dnum time
xtdidregress (opo5r facnum_fixed semp urba litr imps) (treat), group(dnum) time(time) vce(cluster pnum)

gen group_na = group
recode group_na (9=.)

gen rt = time - group_na
gen never_treat = (group == 9)

tab rt

* Create dummy
forvalues k = 7(-1)1 {
	gen lead_`k' = rt == -`k'
	}
forvalues k = 0/4 {
	gen lag_`k' = rt == `k'
	}
	

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

	reghdfe `m' lead_7-lead_2 lag_0-lag_4 facnum_fixed semp urba litr imps, absorb(dnum time) vce(cluster pnum)
	estimates store tbl_sfig9_twfe_`m'

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' (Dynamic TWFE)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/sfig_09_did_hetero_twfedid_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_sfig9_twfe_opo5r tbl_sfig9_twfe_opu5r tbl_sfig9_twfe_ipo5r ///
	   tbl_sfig9_twfe_ipu5r tbl_sfig9_twfe_iplos tbl_sfig9_twfe_delr ///
using "stata/output/sfig_09_did_hetero_twfedid.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace
	
* ----------------------------------------------------------
* Two-stage DiD (Gardner 2021)
* ----------------------------------------------------------

/*
Install user-written package:
net install did2s, replace from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/")

Resource: https://github.com/kylebutts/did2s_stata

The 'did2s' is available in both Stata and R, but they produce different results, even when using the replication materials provided by the author. The R version produces estimates that are closer to the true treatment effect.
*/


* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

	did2s `m', first_stage(i.dnum i.time facnum_fixed semp urba litr imps) second_stage(lead_* lag_*) treatment(treat) cluster(pnum)
	estimates store tbl_sfig9_2stg_`m'

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' Two-stage DiD)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/sfig_09_did_hetero_2stgdid_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_sfig9_2stg_opo5r tbl_sfig9_2stg_opu5r tbl_sfig9_2stg_ipo5r ///
	   tbl_sfig9_2stg_ipu5r tbl_sfig9_2stg_iplos tbl_sfig9_2stg_delr ///
using "stata/output/sfig_09_did_hetero_2stgdid.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_opo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_opo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_opo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_opo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_opu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_opu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_opu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_opu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_ipo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_ipo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_ipo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_ipo5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_ipu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_ipu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_ipu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_ipu5r.gph" ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_iplos.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_iplos.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_iplos.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_iplos.gph" ///
"stata/graph/temp/sfig_09_did_hetero_jwdid_delr.gph" ///
"stata/graph/temp/sfig_09_did_hetero_stkdid_delr.gph" ///
"stata/graph/temp/sfig_09_did_hetero_twfedid_delr.gph" ///
"stata/graph/temp/sfig_09_did_hetero_2stgdid_delr.gph" ///
, rows(6) cols(4) iscale(1) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(37) ysize(60)
graph export "stata/graph/sfig_09_did_hetero.pdf", replace


* ==============================================================================  
* Appendix Figure A7: DiD analysis using count data
* ==============================================================================  

* ----------------------------------------------------------
* All (PH+DH+HC)
* ----------------------------------------------------------

* Import Lao poverty data
import excel "$dir/data import/lao_poverty_2015.xlsx", sheet("Lao Poverty 2015 by district") firstrow clear

* Keep only selected variables
keep dnum Poverty_He Poverty_Ga Poverty_Se Improved_S Lit_rate_1564old Self_employment_rate Urban_popu Area

* Scale selected variables by 0.01
foreach var in Poverty_He Poverty_Ga Poverty_Se Improved_S Lit_rate_1564old Self_employment_rate Urban_popu {
    replace `var' = `var' / 100
}

* Rename variables
rename Poverty_He pov_he
rename Poverty_Ga pov_ga
rename Poverty_Se pov_se
rename Improved_S imps
rename Lit_rate_1564old litr
rename Self_employment_rate semp
rename Urban_popu urba
rename Area area

* Save poverty dataset
save stata/temp/pvty_d, replace

* Import main dataset
import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali", "1000 PH Vientiane")

* Generate dnum based on pnum and faclevel
replace dnum = 205 if pnum == 2  & faclevel == "PH"
replace dnum = 301 if pnum == 3  & faclevel == "PH"
replace dnum = 401 if pnum == 4  & faclevel == "PH"
replace dnum = 501 if pnum == 5  & faclevel == "PH"
replace dnum = 601 if pnum == 6  & faclevel == "PH"
replace dnum = 701 if pnum == 7  & faclevel == "PH"
replace dnum = 801 if pnum == 8  & faclevel == "PH"
replace dnum = 901 if pnum == 9  & faclevel == "PH"
// pnum == 10 already dropped
replace dnum = 1101 if pnum == 11 & faclevel == "PH"
replace dnum = 1201 if pnum == 12 & faclevel == "PH"
replace dnum = 1301 if pnum == 13 & faclevel == "PH"
replace dnum = 1401 if pnum == 14 & faclevel == "PH"
replace dnum = 1501 if pnum == 15 & faclevel == "PH"
replace dnum = 1601 if pnum == 16 & faclevel == "PH"
replace dnum = 1702 if pnum == 17 & faclevel == "PH"
replace dnum = 1801 if pnum == 18 & faclevel == "PH"

* Collapse (aggregate) the data
collapse (sum) opo5 opu5 ipo5 ipu5 ipday smsrgo5 smsrgu5 anc pnc del facnum, ///
    by(province year quarter pnum dnum time group treat event)

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

sort province dnum quarter

* Merge poverty data by dnum
merge m:1 dnum using stata/temp/pvty_d

* Keep only matched or all if needed
drop if _merge == 2		// PH level data
drop _merge

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig10_all_`m'

	csdid_plot, title("`t' (All)", size(vsmall)) legend(off) ///
    xlabel(-7(1)4, labsize(tiny)) ///
    ylabel(, labsize(tiny)) ///
    xtitle("Time since NHI introduction", size(tiny)) ///
    ytitle("ATT", size(tiny))
    graph save "stata/graph/temp/sfig_10_did_count_all_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig10_all_opo5 tbl_sfig10_all_opu5 tbl_sfig10_all_ipo5 ///
	   tbl_sfig10_all_ipu5 tbl_sfig10_all_ipday tbl_sfig10_all_del ///
using "stata/output/sfig_10_did_count_all.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* District-level (DH+HC)
* ----------------------------------------------------------

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig10_dhhc_`m'

	csdid_plot, title("`t' (DH+HC)", size(vsmall)) legend(off) ///
    xlabel(-7(1)4, labsize(tiny)) ///
    ylabel(, labsize(tiny)) ///
    xtitle("Time since NHI introduction", size(tiny)) ///
    ytitle("ATT", size(tiny))
    graph save "stata/graph/temp/sfig_10_did_count_dhhc_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig10_dhhc_opo5 tbl_sfig10_dhhc_opu5 tbl_sfig10_dhhc_ipo5 ///
	   tbl_sfig10_dhhc_ipu5 tbl_sfig10_dhhc_ipday tbl_sfig10_dhhc_del ///
using "stata/output/sfig_10_did_count_dhhc.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* District hospitals (DH only)
* ----------------------------------------------------------

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
keep if faclevel == "DH"	// facnum is 1 and needs to be dropped
drop if district == "1210 Khounkham"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig10_dh_`m'

	csdid_plot, title("`t' (DH)", size(vsmall)) legend(off) ///
    xlabel(-7(1)4, labsize(tiny)) ///
    ylabel(, labsize(tiny)) ///
    xtitle("Time since NHI introduction", size(tiny)) ///
    ytitle("ATT", size(tiny))
    graph save "stata/graph/temp/sfig_10_did_count_dh_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig10_dh_opo5 tbl_sfig10_dh_opu5 tbl_sfig10_dh_ipo5 ///
	   tbl_sfig10_dh_ipu5 tbl_sfig10_dh_ipday tbl_sfig10_dh_del ///
using "stata/output/sfig_10_did_count_dh.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* ----------------------------------------------------------
* Health centers (HC only)
* ----------------------------------------------------------

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if quarter <= 2017.3
keep if province != "01 Vientiane Capital"
keep if faclevel == "HC"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig10_hc_`m'

	csdid_plot, title("`t' (HC)", size(vsmall)) legend(off) ///
    xlabel(-7(1)4, labsize(tiny)) ///
    ylabel(, labsize(tiny)) ///
    xtitle("Time since NHI introduction", size(tiny)) ///
    ytitle("ATT", size(tiny))
    graph save "stata/graph/temp/sfig_10_did_count_hc_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_sfig10_hc_opo5 tbl_sfig10_hc_opu5 tbl_sfig10_hc_ipo5 ///
	   tbl_sfig10_hc_ipu5 tbl_sfig10_hc_ipday tbl_sfig10_hc_del ///
using "stata/output/sfig_10_did_count_hc.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_10_did_count_all_opo5.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_opo5.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_opo5.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_opo5.gph" ///
"stata/graph/temp/sfig_10_did_count_all_opu5.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_opu5.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_opu5.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_opu5.gph" ///
"stata/graph/temp/sfig_10_did_count_all_ipo5.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_ipo5.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_ipo5.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_ipo5.gph" ///
"stata/graph/temp/sfig_10_did_count_all_ipu5.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_ipu5.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_ipu5.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_ipu5.gph" ///
"stata/graph/temp/sfig_10_did_count_all_ipday.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_ipday.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_ipday.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_ipday.gph" ///
"stata/graph/temp/sfig_10_did_count_all_del.gph" ///
"stata/graph/temp/sfig_10_did_count_dhhc_del.gph" ///
"stata/graph/temp/sfig_10_did_count_dh_del.gph" ///
"stata/graph/temp/sfig_10_did_count_hc_del.gph" ///
, rows(6) cols(4) iscale(1) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(37) ysize(60)
graph export "stata/graph/sfig_10_did_count.pdf", replace


* ==============================================================================  
* Appendix Figure A8: DiD analysis using monthly data
* ==============================================================================  

import excel "$dir/data export/02_clean_8_cln_m_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
describe month
keep if month <= clock("01sep2017 00:00:00", "DMYhms")
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

* ----------------------------------------------------------
* Base model
* ----------------------------------------------------------

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
        method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) ///
        rseed(555) saverif(_csdidrif_) replace // sci applied
	estimates store tbl_sfig11_base_`m'

	csdid_plot, title("`t' (base)", size(vsmall)) legend(off) ///
		xlabel(-23(1)14, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_11_did_month_base_`m'.gph", replace

* ----------------------------------------------------------
* Balanced NHI exposure (e=0-2)
* ----------------------------------------------------------

* Equation 3.11 from Callaway and Sant'Anna (2021) is implementable in R but not available through 'xthdidregress' or 'csdid'.

* ----------------------------------------------------------
* Last-treated units as control
* ----------------------------------------------------------

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig11_last_`m'

	csdid_plot, title("`t' (last)", size(vsmall)) legend(off) ///
		xlabel(-23(1)14, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
    graph save "stata/graph/temp/sfig_11_did_month_last_`m'.gph", replace

* ----------------------------------------------------------
* No covariate (unconditional Parallel Trend Assumption)
* ----------------------------------------------------------

csdid `m', ivar(dnum) time(time) gvar(group_0) method(dripw) agg(event) long2 notyet wboot(reps(1000)) cluster(pnum) rseed(555) replace
	estimates store tbl_sfig11_nocov_`m'

csdid_plot, title("`t' (no covarite)", size(vsmall)) legend(off) ///
	xlabel(-23(1)14, labsize(tiny)) ///
	ylabel(, labsize(tiny)) ///
	xtitle("Time since NHI introduction", size(tiny)) ///
	ytitle("ATT", size(tiny))
	graph save "stata/graph/temp/sfig_11_did_month_nocov_`m'.gph", replace

}

* Export NHI estimates
esttab tbl_sfig11_base_opo5r tbl_sfig11_base_opu5r tbl_sfig11_base_ipo5r ///
	   tbl_sfig11_base_ipu5r tbl_sfig11_base_iplos tbl_sfig11_base_delr ///
using "stata/output/sfig_11_did_month_base.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

esttab tbl_sfig11_last_opo5r tbl_sfig11_last_opu5r tbl_sfig11_last_ipo5r ///
	   tbl_sfig11_last_ipu5r tbl_sfig11_last_iplos tbl_sfig11_last_delr ///
using "stata/output/sfig_11_did_month_last.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

esttab tbl_sfig8_nocov_opo5r tbl_sfig8_nocov_opu5r tbl_sfig8_nocov_ipo5r ///
	   tbl_sfig8_nocov_ipu5r tbl_sfig8_nocov_iplos tbl_sfig8_nocov_delr ///
using "stata/output/sfig_11_did_month_nocov.xlsx", ///
mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") replace

* Generate graph
graph combine ///
"stata/graph/temp/sfig_11_did_month_base_opo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_last_opo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_opo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_base_opu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_last_opu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_opu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_base_ipo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_last_ipo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_ipo5r.gph" ///
"stata/graph/temp/sfig_11_did_month_base_ipu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_last_ipu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_ipu5r.gph" ///
"stata/graph/temp/sfig_11_did_month_base_iplos.gph" ///
"stata/graph/temp/sfig_11_did_month_last_iplos.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_iplos.gph" ///
"stata/graph/temp/sfig_11_did_month_base_delr.gph" ///
"stata/graph/temp/sfig_11_did_month_last_delr.gph" ///
"stata/graph/temp/sfig_11_did_month_nocov_delr.gph" ///
, rows(6) cols(3) iscale(1) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) ///
 xsize(35) ysize(65)
graph export "stata/graph/sfig_11_did_month.pdf", replace

