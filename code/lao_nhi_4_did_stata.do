/*=============================================================================*
 Title: The impact of national health insurance on health services utilization in Lao People's Democratic Republic (Lao PDR)
 Author: Sekeun Daniel Yu (yus109@mcmaster.ca), Michel Grignon, Godefroy Emmanuel Guindon, Jean-Éric Tarride
 Date: September/2026
*=============================================================================*/

clear all
global dir "/enter-directory-here/"
cd "$dir"
version 19

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
*==============================================================================*
* Set-up
*==============================================================================*

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
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 


*==============================================================================*
* Figure 3. Difference-in-differences estimates by utilization outcome
*==============================================================================*

/*
Group-Time ATT (gt-ATT) estimator (Sant'Anna and Callaway 2021)
Dynamic treatment effect with simultaneous confidence intervals (sci)
Intervention group: Q3/2016, Q4/2016, Q1/2017, Q3/2017
Control group: not-yet-treated (Q4/2016, Q1/2017, Q3/2017) + last-treated (Q4/2017)
*/

ssc install drdid, all replace
//ssc install csdid, all replace
net install csdid, from("https://raw.githubusercontent.com/pedrohcgs/csdid-stata/main") replace
// help csdid
// which csdid: *! csdid 2.0.0 08sep2026

set graphics off

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	* No bootstrap
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) ///
		cluster(pnum) analytical pointwise // analytical clustered SE
	estat event, post	// Event-study and post_avg are pointwise CI
	estimates store tbl_fig3_`m'_noboot
	estat tidy, saving("stata/temp/fig_03_did_`m'_noboot_tidy.dta") replace
	
	* Bootstrap
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
		cluster(pnum)	// multiplicative wild bootstrap
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_fig3_`m'
	estat tidy, saving("stata/temp/fig_03_did_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_03_did_`m'.gph", replace
	
	*---------------------------------------------------------*
	* HonestDiD
	*---------------------------------------------------------*
	
	if inlist("`m'", "opo5r", "opu5r", "ipo5r") {

        * Impact at e=0: relative magnitude
        local plotopts xtitle(Mbar) ytitle(95% Robust CI)

        honestdid, ///
            pre(1/6) ///
            post(8/12) ///
            mvec(0.5(0.5)2) ///
            delta(rm) ///
	        matasave(hon_e0_rm) ///
            coefplot `plotopts'

        graph save "stata/graph/temp/honest_`m'_e0_rm.gph", replace
        graph export "stata/graph/temp/honest_`m'_e0_rm.pdf", replace

        * Impact at e=0: smoothness
        local plotopts xtitle(M) ytitle(95% Robust CI)

        honestdid, ///
            pre(1/6) ///
            post(8/12) ///
            mvec(0(5)20) ///
            delta(sd) ///
			matasave(hon_e0_sd) ///
            coefplot `plotopts'

        graph save "stata/graph/temp/honest_`m'_e0_sd.gph", replace
        graph export "stata/graph/temp/honest_`m'_e0_sd.pdf", replace

        * Five-quarter average: relative magnitude
        matrix l_vec = 0.2 \ 0.2 \ 0.2 \ 0.2 \ 0.2
        local plotopts xtitle(Mbar) ytitle(95% Robust CI)

        honestdid, ///
            pre(1/6) ///
            post(8/12) ///
            l_vec(l_vec) ///
            mvec(0.5(0.5)2) ///
            delta(rm) ///
			matasave(hon_avg_rm) ///
            coefplot `plotopts'

        graph save "stata/graph/temp/honest_`m'_avg_rm.gph", replace
        graph export "stata/graph/temp/honest_`m'_avg_rm.pdf", replace

        * Five-quarter average: smoothness
        local plotopts xtitle(M) ytitle(95% Robust CI)

        honestdid, ///
            pre(1/6) ///
            post(8/12) ///
            l_vec(l_vec) ///
            mvec(0(5)20) ///
            delta(sd) ///
	        matasave(hon_avg_sd) ///
            coefplot `plotopts'

        graph save "stata/graph/temp/honest_`m'_avg_sd.gph", replace
        graph export "stata/graph/temp/honest_`m'_avg_sd.pdf", replace
		
		foreach target in e0 avg {
		foreach restriction in rm sd {

        mata: st_matrix("hon_ci", hon_`target'_`restriction'.CI)
        matrix colnames hon_ci = M lb ub

        esttab matrix(hon_ci, fmt(3)) using ///
            "stata/output/honest_`m'_`target'_`restriction'.csv", ///
            csv nomtitles replace
		}
		}
    }
}

* Export NHI estimates (bootstrap simultaneous CI, no bootstrap CI)
preserve
local writeopt replace
foreach m of local varlist {

    * Bootstrap
    use "stata/temp/fig_03_did_`m'_tidy.dta", clear
    export excel using "stata/output/fig_03_did_bootstrap.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'

    * No bootstrap
    use "stata/temp/fig_03_did_`m'_noboot_tidy.dta", clear
    export excel using "stata/output/fig_03_did_noboot.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'

    local writeopt sheetreplace
}
restore

* Export NHI estimates (no bootstrap pointwise CI, post_avg)
esttab tbl_fig3_opo5r_noboot tbl_fig3_opu5r_noboot tbl_fig3_ipo5r_noboot ///
	   tbl_fig3_ipu5r_noboot tbl_fig3_iplos_noboot tbl_fig3_delr_noboot ///
using "stata/output/fig_03_did_pointwise_noboot.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Export NHI estimates (bootstrap pointwise CI, post_avg)
esttab tbl_fig3_opo5r tbl_fig3_opu5r tbl_fig3_ipo5r ///
	   tbl_fig3_ipu5r tbl_fig3_iplos tbl_fig3_delr ///
using "stata/output/fig_03_did_pointwise_bootstrap.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
foreach m of local varlist {
    graph use "stata/graph/temp/fig_03_did_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_03_did_`m'.gph", replace
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_03_did_opo5r.gph" ///
"stata/graph/temp/fig_03_did_opu5r.gph" ///
"stata/graph/temp/fig_03_did_ipo5r.gph" ///
"stata/graph/temp/fig_03_did_ipu5r.gph" ///
"stata/graph/temp/fig_03_did_iplos.gph" ///
"stata/graph/temp/fig_03_did_delr.gph" ///
, rows(3) cols(2) iscale(0.5) graphregion(margin(0 0 0 0)) ///
xsize(7) ysize(8)
graph export "stata/graph/fig_03_did_main.pdf", replace

* Generate graph
graph combine ///
"stata/graph/temp/honest_opo5r_e0_rm.gph" ///
"stata/graph/temp/honest_opo5r_e0_sd.gph" ///
"stata/graph/temp/honest_ipo5r_e0_rm.gph" ///
"stata/graph/temp/honest_ipo5r_e0_sd.gph" ///
"stata/graph/temp/honest_opu5r_e0_rm.gph" ///
"stata/graph/temp/honest_opu5r_e0_sd.gph" ///
, rows(3) cols(2) ///
iscale(0.5) graphregion(margin(0 0 0 0)) ///
xsize(6) ysize(8)
graph save "stata/graph/fig_b_05_honest_e0.gph", replace
graph export "stata/graph/fig_b_05_honest_e0.pdf", replace

* Generate graph
graph combine ///
"stata/graph/temp/honest_opo5r_avg_rm.gph" ///
"stata/graph/temp/honest_opo5r_avg_sd.gph" ///
"stata/graph/temp/honest_ipo5r_avg_rm.gph" ///
"stata/graph/temp/honest_ipo5r_avg_sd.gph" ///
"stata/graph/temp/honest_opu5r_avg_rm.gph" ///
"stata/graph/temp/honest_opu5r_avg_sd.gph" ///
, rows(3) cols(2) ///
iscale(0.5) graphregion(margin(0 0 0 0)) ///
xsize(6) ysize(8)
graph save "stata/graph/fig_b_05_honest_avg.gph", replace
graph export "stata/graph/fig_b_05_honest_avg.pdf", replace

/*
The point estimates are identical in R and Stata, with consistent conclusions regarding statistical significance. The simultaneous confidence intervals likely differ because R and Stata use different bootstrap implementations and random number generation processes.
*/


*==============================================================================*
* Appendix Figure A2: DiD analysis by facility type
*==============================================================================*

*---------------------------------------------------------*
* District hospitals
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
		cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa2_dh_`m'
	estat tidy, saving("stata/temp/fig_a_02_did_facility_dh_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_02_did_facility_dh_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_02_did_facility_dh_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_02_did_facility_dh_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa2_dh_opo5r tbl_figa2_dh_opu5r tbl_figa2_dh_ipo5r ///
	   tbl_figa2_dh_ipu5r tbl_figa2_dh_iplos tbl_figa2_dh_delr ///
using "stata/output/fig_a_02_did_facility_dh_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* Health centers
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
		cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa2_hc_`m'
	estat tidy, saving("stata/temp/fig_a_02_did_facility_hc_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_02_did_facility_hc_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_02_did_facility_hc_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_02_did_facility_hc_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa2_hc_opo5r tbl_figa2_hc_opu5r tbl_figa2_hc_ipo5r ///
	   tbl_figa2_hc_ipu5r tbl_figa2_hc_iplos tbl_figa2_hc_delr ///
using "stata/output/fig_a_02_did_facility_hc_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
foreach m of local varlist {
    graph use "stata/graph/temp/fig_a_02_did_facility_dh_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_02_did_facility_dh_`m'.gph", replace
    
	graph use "stata/graph/temp/fig_a_02_did_facility_hc_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_02_did_facility_hc_`m'.gph", replace
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_02_did_facility_dh_opo5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_opo5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_dh_opu5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_opu5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_dh_ipo5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_ipo5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_dh_ipu5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_ipu5r.gph" ///
"stata/graph/temp/fig_a_02_did_facility_dh_iplos.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_iplos.gph" ///
"stata/graph/temp/fig_a_02_did_facility_dh_delr.gph" ///
"stata/graph/temp/fig_a_02_did_facility_hc_delr.gph" ///
, rows(6) cols(2) iscale(0.5) graphregion(margin(0 0 0 0)) xsize(7) ysize(16)
graph export "stata/graph/fig_a_02_did_facility.pdf", replace


*==============================================================================*
* Appendix Figure A3: DiD analysis by district-level poverty status
*==============================================================================*

ssc install event_plot, replace

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0)

*---------------------------------------------------------*
/** Group-Time DiD (Callaway and Sant'Anna (2021)) **/
*---------------------------------------------------------*

/****** High-poverty area ******/

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 1
	
	csdid `m', ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000)) ///
		cluster(pnum) rseed(555)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa3_cs_high_`m'
	estat tidy, saving("stata/temp/fig_a_03_did_poverty_cs_high_`m'_tidy.dta") replace

	csdid_plot
	graph save "stata/graph/temp/fig_a_03_did_poverty_cs_high_`m'.gph", replace
	
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_03_did_poverty_cs_high_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_03_did_poverty_cs_high_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa3_cs_high_opo5r tbl_figa3_cs_high_opu5r tbl_figa3_cs_high_ipo5r ///
	   tbl_figa3_cs_high_ipu5r tbl_figa3_cs_high_iplos tbl_figa3_cs_high_delr ///
using "stata/output/fig_a_03_did_poverty_cs_high_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

/****** Low-poverty area ******/

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 0
	
	csdid `m', ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000)) ///
		cluster(pnum) rseed(555)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa3_cs_low_`m'
	estat tidy, saving("stata/temp/fig_a_03_did_poverty_cs_low_`m'_tidy.dta") replace

	csdid_plot
	graph save "stata/graph/temp/fig_a_03_did_poverty_cs_low_`m'.gph", replace
	
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_03_did_poverty_cs_low_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_03_did_poverty_cs_low_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa3_cs_low_opo5r tbl_figa3_cs_low_opu5r tbl_figa3_cs_low_ipo5r ///
	   tbl_figa3_cs_low_ipu5r tbl_figa3_cs_low_iplos tbl_figa3_cs_low_delr ///
using "stata/output/fig_a_03_did_poverty_cs_low_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
foreach m of local varlist {
    graph use "stata/graph/temp/fig_a_03_did_poverty_cs_high_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_03_did_poverty_cs_high_`m'.gph", replace
    
	graph use "stata/graph/temp/fig_a_03_did_poverty_cs_low_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_03_did_poverty_cs_low_`m'.gph", replace
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_opo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_opo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_opu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_opu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_ipo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_ipo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_ipu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_ipu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_iplos.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_iplos.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_high_delr.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_cs_low_delr.gph" ///
, rows(6) cols(2) iscale(0.5) graphregion(margin(0 0 0 0)) xsize(7) ysize(16)
graph export "stata/graph/fig_a_03_did_poverty_cs.pdf", replace

*---------------------------------------------------------*
* Interaction-weighted DiD estimator (Sun and Abraham 2021)
*---------------------------------------------------------*
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

/****** High-poverty area ******/

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 1

	eventstudyinteract `m' lead_7-lead_2 lag_0-lag_4, cohort(group_na) ///
		control_cohort(never_treat) absorb(i.dnum i.time) vce(cluster pnum) ///
		covariates(i.time#c.facnum i.time#c.semp i.time#c.urba ///
				   i.time#c.litr i.time#c.imps)

	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Time since NHI introduction", size(vsmall)) ytitle("ATT", size(vsmall)) xlabel(-7(1)4, labsize(vsmall)) title("`t' (high-poverty)", size(small))) stub_lead(lead_#) stub_lag(lag_#) trimlag(4) trimlead(7) together
	graph save "stata/graph/temp/fig_a_03_did_poverty_as_high_`m'.gph", replace
	
	matrix b = e(b_iw)
	matrix V = e(V_iw)
	ereturn post b V
	esttab
	estimates store tbl_figa3_high_`m'
	
	restore
}

* Export NHI estimates
esttab tbl_figa3_high_opo5r tbl_figa3_high_opu5r tbl_figa3_high_ipo5r ///
	   tbl_figa3_high_ipu5r tbl_figa3_high_iplos tbl_figa3_high_delr ///
using "stata/output/fig_a_03_did_poverty_as_high_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

/****** Low-poverty area ******/

local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	preserve
	keep if pov_med == 0

	eventstudyinteract `m' lead_7-lead_2 lag_0-lag_4, cohort(group_na) ///
	control_cohort(never_treat) absorb(i.dnum i.time) vce(cluster pnum) ///
	covariates(i.time#c.facnum i.time#c.semp i.time#c.urba ///
			   i.time#c.litr i.time#c.imps)

	event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Time since NHI introduction", size(vsmall)) ytitle("ATT", size(vsmall)) xlabel(-7(1)4, labsize(vsmall)) title("`t' (low-poverty)", size(small))) stub_lead(lead_#) stub_lag(lag_#) trimlag(4) trimlead(7) together
	graph save "stata/graph/temp/fig_a_03_did_poverty_as_low_`m'.gph", replace
	
	matrix b = e(b_iw)
	matrix V = e(V_iw)
	ereturn post b V
	esttab
	estimates store tbl_figa3_low_`m'
	
	restore
}

* Export NHI estimates
esttab tbl_figa3_low_opo5r tbl_figa3_low_opu5r tbl_figa3_low_ipo5r ///
	   tbl_figa3_low_ipu5r tbl_figa3_low_iplos tbl_figa3_low_delr ///
using "stata/output/fig_a_03_did_poverty_as_low_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_opo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_opo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_opu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_opu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_ipo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_ipo5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_ipu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_ipu5r.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_iplos.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_iplos.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_high_delr.gph" ///
"stata/graph/temp/fig_a_03_did_poverty_as_low_delr.gph" ///
, rows(6) cols(2) iscale(0.5) graphregion(margin(0 0 0 0)) xsize(7) ysize(16)
graph export "stata/graph/fig_a_03_did_poverty_as.pdf", replace


*==============================================================================*
* Appendix Figure A4: DiD analysis by intervention group
*==============================================================================*

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
		method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
		cluster(pnum)
	estimates store tbl_figa4_`m'
	estat tidy, saving("stata/temp/fig_a_04_did_group_`m'_tidy.dta") replace

	preserve
	use "stata/temp/fig_a_04_did_group_`m'_tidy.dta", clear

	foreach g in 4 5 6 8 {
		twoway ///
			(rcap conf_high conf_low time if group == `g' & time < `g', ///
				lcolor(navy)) ///
			(scatter estimate time if group == `g' & time < `g', ///
				mcolor(navy)) ///
			(rcap conf_high conf_low time if group == `g' & time >= `g', ///
				lcolor(maroon)) ///
			(scatter estimate time if group == `g' & time >= `g', ///
				mcolor(maroon)), ///
			yline(0, lpattern(dash) lcolor(gs8)) ///
			xtitle("Time") ytitle("ATT(g,t)") ///
			legend(off)
			
		graph save ///
			"stata/graph/temp/fig_a_04_did_group_`g'_`m'.gph", replace
	}
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_04_did_group_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_04_did_group_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
}
restore
	
* Export NHI estimates (pointwise CI)
esttab tbl_figa4_opo5r tbl_figa4_opu5r tbl_figa4_ipo5r ///
	   tbl_figa4_ipu5r tbl_figa4_iplos tbl_figa4_delr ///
using "stata/output/fig_a_04_did_group_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

graph combine ///
"stata/graph/temp/fig_a_04_did_group_4_opo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_opo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_opo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_opo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_4_opu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_opu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_opu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_opu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_4_ipo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_ipo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_ipo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_ipo5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_4_ipu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_ipu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_ipu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_ipu5r.gph" ///
"stata/graph/temp/fig_a_04_did_group_4_iplos.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_iplos.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_iplos.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_iplos.gph" ///
"stata/graph/temp/fig_a_04_did_group_4_delr.gph" ///
"stata/graph/temp/fig_a_04_did_group_5_delr.gph" ///
"stata/graph/temp/fig_a_04_did_group_6_delr.gph" ///
"stata/graph/temp/fig_a_04_did_group_8_delr.gph" ///
, rows(6) cols(4) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) ///
 xsize(14) ysize(16)
graph export "stata/graph/fig_a_04_did_group.pdf", replace


*==============================================================================*
* Appendix Figure A5: DiD analysis with adjusted model specifications
*==============================================================================*

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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

*---------------------------------------------------------*
* Base model
*---------------------------------------------------------*

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
	method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
	cluster(pnum) // sci applied
estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa5_base_`m'
estat tidy, saving("stata/temp/fig_a_05_did_adjust_base_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_05_did_adjust_base_`m'.gph", replace

*---------------------------------------------------------*
* Balanced NHI exposure (e=0-2)
*---------------------------------------------------------*

* Equation 3.11 from Callaway and Sant'Anna (2021) is implementable in R but not available through 'xthdidregress' or 'csdid'.

*---------------------------------------------------------*
* Last-treated units as control
*---------------------------------------------------------*

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
	method(dr) nevertreated base_period(universal) wboot(reps(1000) rseed(555)) ///
	cluster(pnum)
estat event, dropmissing post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa5_last_`m'
estat tidy, saving("stata/temp/fig_a_05_did_adjust_last_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_05_did_adjust_last_`m'.gph", replace

*---------------------------------------------------------*
* No covariate (unconditional Parallel Trend Assumption)
*---------------------------------------------------------*

csdid `m', ivar(dnum) time(time) gvar(group_0) method(dr) ///
	notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
	cluster(pnum)
estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa5_nocov_`m'
estat tidy, saving("stata/temp/fig_a_05_did_adjust_nocov_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_05_did_adjust_nocov_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_05_did_adjust_base_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_05_did_adjust_base_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'

	use "stata/temp/fig_a_05_did_adjust_last_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_05_did_adjust_last_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'

	use "stata/temp/fig_a_05_did_adjust_nocov_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_05_did_adjust_nocov_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa5_base_opo5r tbl_figa5_base_opu5r tbl_figa5_base_ipo5r ///
	   tbl_figa5_base_ipu5r tbl_figa5_base_iplos tbl_figa5_base_delr ///
	   using "stata/output/fig_a_05_did_adjust_base_pointwise.csv", ///
	       cells("b(fmt(3)) ci(fmt(3))") ///
		   mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

esttab tbl_figa5_last_opo5r tbl_figa5_last_opu5r tbl_figa5_last_ipo5r ///
	   tbl_figa5_last_ipu5r tbl_figa5_last_iplos tbl_figa5_last_delr ///
	   using "stata/output/fig_a_05_did_adjust_last_pointwise.csv", ///
	       cells("b(fmt(3)) ci(fmt(3))") ///
		   mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

esttab tbl_figa5_nocov_opo5r tbl_figa5_nocov_opu5r tbl_figa5_nocov_ipo5r ///
	   tbl_figa5_nocov_ipu5r tbl_figa5_nocov_iplos tbl_figa5_nocov_delr ///
	   using "stata/output/fig_a_05_did_adjust_nocov_pointwise.csv", ///
	       cells("b(fmt(3)) ci(fmt(3))") ///
		   mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
local specs base last nocov
foreach s of local specs {
foreach m of local varlist {
    graph use "stata/graph/temp/fig_a_05_did_adjust_`s'_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_05_did_adjust_`s'_`m'.gph", replace
	}
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_05_did_adjust_base_opo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_opo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_opo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_base_opu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_opu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_opu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_base_ipo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_ipo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_ipo5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_base_ipu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_ipu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_ipu5r.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_base_iplos.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_iplos.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_iplos.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_base_delr.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_last_delr.gph" ///
"stata/graph/temp/fig_a_05_did_adjust_nocov_delr.gph" ///
, rows(6) cols(3) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) ///
 xsize(6) ysize(8)
graph export "stata/graph/fig_a_05_did_adjust.pdf", replace


*==============================================================================*
* Appendix Figure A6: DiD analysis using alternative heterogeneity-robust estimators
*==============================================================================*

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

*---------------------------------------------------------*
* Extended TWFE DID (Wooldrigde 2021)
*---------------------------------------------------------*

/*
Install user-written package:
ssc install jwdid, all replace

Resource: https://doi.org/10.1093/ectj/utad016
		  https://friosavila.github.io/app_metrics/app_metrics11.html
		  
The point estimates are the same in R and Stata, and the confidence intervals are nearly identical.
*/

* Loop to estimate NHI effects
local varlist opo5r opu5r ipo5r ipu5r iplos delr
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient LOS" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	

	jwdid `m' facnum_fixed semp urba litr imps, ivar(dnum) tvar(time) gvar(group_0) cluster(pnum)
	estat event, estore(tbl_figa6_jw_`m')
	
	estat plot, title("`t' (Extended TWFE)", size(vsmall)) legend(off) ///
		xlabel(-7(1)4, labsize(tiny)) ///
        ylabel(, labsize(tiny)) ///
        xtitle("Time since NHI introduction", size(tiny)) ///
        ytitle("ATT", size(tiny))      
	graph save "stata/graph/temp/fig_a_06_did_hetero_jwdid_`m'.gph", replace
	
}

* Export NHI estimates
esttab tbl_figa6_jw_opo5r tbl_figa6_jw_opu5r tbl_figa6_jw_ipo5r ///
	   tbl_figa6_jw_ipu5r tbl_figa6_jw_iplos tbl_figa6_jw_delr ///
using "stata/output/fig_a_06_did_hetero_jwdid.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* Stacked DiD (Wing 2024)
*---------------------------------------------------------*

/*
Code is from:
https://rawcdn.githack.com/hollina/stacked-did-weights/18a5e1155506cbd754b78f9cef549ac96aef888b/stacked-example-r-and-stata.html

Stacked DiD estimator uses trimmed aggregate ATT which is a weighted average of group-time ATTs from a trimmed set of data. This is analogous to the balanced NHI exposure configuration in Equation 3.11 of Callaway and Sant'Anna (2021).

The point estimates and confidence intervals are the same in R and Stata.
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
	estimates store tbl_figa6_stk_`m'

	* Show results
	esttab tbl_figa6_stk_`m', keep(lead* lag*) se

	* Compute the average post-treatment effect
	lincom (lag_0 + lag_1 + lag_2)/3

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' (Stacked DiD)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/fig_a_06_did_hetero_stkdid_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_figa6_stk_opo5r tbl_figa6_stk_opu5r tbl_figa6_stk_ipo5r ///
	   tbl_figa6_stk_ipu5r tbl_figa6_stk_iplos tbl_figa6_stk_delr ///
using "stata/output/fig_a_06_did_hetero_stkdid.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* Dynamic TWFE DiD
*---------------------------------------------------------*

/*
This estimator is possibly biased in the context of treatment effect heterogeneity. Code is from:
https://lost-stats.github.io/Model_Estimation/Research_Design/event_study.html#:~:text=Difference-in-Differences%20Event%20Study,periods%20in%20your%20respective%20study.

The point estimates and confidence intervals are the same in R and Stata.
*/

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	estimates store tbl_figa6_twfe_`m'

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' (Dynamic TWFE)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/fig_a_06_did_hetero_twfedid_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_figa6_twfe_opo5r tbl_figa6_twfe_opu5r tbl_figa6_twfe_ipo5r ///
	   tbl_figa6_twfe_ipu5r tbl_figa6_twfe_iplos tbl_figa6_twfe_delr ///
using "stata/output/fig_a_06_did_hetero_twfedid.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace
	
*---------------------------------------------------------*
* Two-stage DiD (Gardner 2021)
*---------------------------------------------------------*

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
	estimates store tbl_figa6_2stg_`m'

	* Generate graph
	event_plot, default_look graph_opt(xtitle("Time since NHI introduction", size(tiny)) ytitle("ATT", size(tiny)) ylabel(, labsize(tiny)) xlabel(-7(1)4, labsize(tiny)) title("`t' (Two-stage DiD)", size(vsmall))) stub_lead(lead_#)  stub_lag(lag_#) trimlag(4) trimlead(7) together
 
	graph save "stata/graph/temp/fig_a_06_did_hetero_2stgdid_`m'.gph", replace
}

* Export NHI estimates
esttab tbl_figa6_2stg_opo5r tbl_figa6_2stg_opu5r tbl_figa6_2stg_ipo5r ///
	   tbl_figa6_2stg_ipu5r tbl_figa6_2stg_iplos tbl_figa6_2stg_delr ///
using "stata/output/fig_a_06_did_hetero_2stgdid.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_opo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_opo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_opo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_opo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_opu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_opu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_opu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_opu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_ipo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_ipo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_ipo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_ipo5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_ipu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_ipu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_ipu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_ipu5r.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_iplos.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_iplos.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_iplos.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_iplos.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_jwdid_delr.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_stkdid_delr.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_twfedid_delr.gph" ///
"stata/graph/temp/fig_a_06_did_hetero_2stgdid_delr.gph" ///
, rows(6) cols(4) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(6) ysize(8)
graph export "stata/graph/fig_a_06_did_hetero.pdf", replace


*==============================================================================*
* Appendix Figure A7: DiD analysis using count data
*==============================================================================*

*---------------------------------------------------------*
* All (PH+DH+HC)
*---------------------------------------------------------*

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
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_all_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_all_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_07_did_count_all_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_all_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_all_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_all_opo5 tbl_figa7_all_opu5 tbl_figa7_all_ipo5 ///
	   tbl_figa7_all_ipu5 tbl_figa7_all_ipday tbl_figa7_all_del ///
using "stata/output/fig_a_07_did_count_all_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* District-level (DH+HC)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_dhhc_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_dhhc_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_07_did_count_dhhc_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_dhhc_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_dhhc_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_dhhc_opo5 tbl_figa7_dhhc_opu5 tbl_figa7_dhhc_ipo5 ///
	   tbl_figa7_dhhc_ipu5 tbl_figa7_dhhc_ipday tbl_figa7_dhhc_del ///
using "stata/output/fig_a_07_did_count_dhhc_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* District hospitals (DH only)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
keep if faclevel == "DH"	// facnum is 1 and needs to be dropped
drop if district == "1210 Khounkham"
replace ipo5 = . if district == "0202 Mai"
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
	
	csdid `m' semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_dh_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_dh_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_07_did_count_dh_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_dh_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_dh_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_dh_opo5 tbl_figa7_dh_opu5 tbl_figa7_dh_ipo5 ///
	   tbl_figa7_dh_ipu5 tbl_figa7_dh_ipday tbl_figa7_dh_del ///
using "stata/output/fig_a_07_did_count_dh_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

*---------------------------------------------------------*
* Health centers (HC only)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
	
	csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_hc_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_hc_`m'_tidy.dta") replace

	csdid_plot
    graph save "stata/graph/temp/fig_a_07_did_count_hc_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_hc_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_hc_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_hc_opo5 tbl_figa7_hc_opu5 tbl_figa7_hc_ipo5 ///
	   tbl_figa7_hc_ipu5 tbl_figa7_hc_ipday tbl_figa7_hc_del ///
using "stata/output/fig_a_07_did_count_hc_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
local specs all dhhc dh hc
local outcomes opo5 opu5 ipo5 ipu5 ipday del
foreach s of local specs {
    foreach m of local outcomes {
    graph use "stata/graph/temp/fig_a_07_did_count_`s'_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_07_did_count_`s'_`m'.gph", replace
	}
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_07_did_count_all_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_all_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_all_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_all_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_all_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_all_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_dhhc_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_dh_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_hc_del.gph" ///
, rows(6) cols(4) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(6) ysize(8)
graph export "stata/graph/fig_a_07_did_count.pdf", replace


*==============================================================================*
* Appendix Figure A7: DiD analysis using count data - log
*==============================================================================*

*---------------------------------------------------------*
* All (PH+DH+HC)
*---------------------------------------------------------*

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

* Generate logs after aggregation
foreach v in opo5 opu5 ipo5 ipu5 ipday del {
    gen double l`v' = ln(`v') if `v' > 0
}

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
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
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient days" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	* Raw-count mean at e=-1 among treated districts.
	quietly summarize `m' if group_0 > 0 & event == -1, meanonly
	local baseline = r(mean)

	csdid l`m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_all_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_log_all_`m'_tidy.dta") replace
    
	* Transform to original scale
    preserve
    use "stata/temp/fig_a_07_did_count_log_all_`m'_tidy.dta", clear
    generate double avgunit = `baseline'
    generate str12 outcome = "l`m'"
    foreach v in estimate conf_low conf_high point_conf_low point_conf_high {
        rename `v' `v'_log
        generate double `v' = avgunit * (exp(`v'_log) - 1)
		}
    rename std_error std_error_log
    rename statistic statistic_log
    rename p_value p_value_log
    order outcome avgunit type term event_time estimate conf_low conf_high
    save "stata/temp/fig_a_07_did_count_log_all_`m'_tidy.dta", replace

    twoway (rcap conf_low conf_high event_time) ///
        (scatter estimate event_time) if !missing(event_time), ///
        yline(0) xline(0) xlabel(-7(1)4) legend(off) ///
        title("`t'") xtitle("Time since NHI introduction, quarter") ///
        ytitle("ATT (count units)")
    graph save "stata/graph/temp/fig_a_07_did_count_log_all_`m'.gph", replace
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_log_all_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_log_all_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_all_opo5 tbl_figa7_all_opu5 tbl_figa7_all_ipo5 ///
	   tbl_figa7_all_ipu5 tbl_figa7_all_ipday tbl_figa7_all_del ///
using "stata/output/fig_a_07_did_count_log_all_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient days" "Institutional birth") plain replace

*---------------------------------------------------------*
* District-level (DH+HC)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_dl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient days" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	* Raw-count mean at e=-1 among treated districts.
	quietly summarize `m' if group_0 > 0 & event == -1, meanonly
	local baseline = r(mean)

	csdid l`m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_dhhc_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_log_dhhc_`m'_tidy.dta") replace

    * Transform to original scale
    preserve
    use "stata/temp/fig_a_07_did_count_log_dhhc_`m'_tidy.dta", clear
    generate double avgunit = `baseline'
    generate str12 outcome = "l`m'"
    foreach v in estimate conf_low conf_high point_conf_low point_conf_high {
        rename `v' `v'_log
        generate double `v' = avgunit * (exp(`v'_log) - 1)
		}
    rename std_error std_error_log
    rename statistic statistic_log
    rename p_value p_value_log
    order outcome avgunit type term event_time estimate conf_low conf_high
    save "stata/temp/fig_a_07_did_count_log_dhhc_`m'_tidy.dta", replace

    twoway (rcap conf_low conf_high event_time) ///
        (scatter estimate event_time) if !missing(event_time), ///
        yline(0) xline(0) xlabel(-7(1)4) legend(off) ///
        title("`t'") xtitle("Time since NHI introduction, quarter") ///
        ytitle("ATT (count units)")
    graph save "stata/graph/temp/fig_a_07_did_count_log_dhhc_`m'.gph", replace
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_log_dhhc_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_log_dhhc_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_dhhc_opo5 tbl_figa7_dhhc_opu5 tbl_figa7_dhhc_ipo5 ///
	   tbl_figa7_dhhc_ipu5 tbl_figa7_dhhc_ipday tbl_figa7_dhhc_del ///
using "stata/output/fig_a_07_did_count_log_dhhc_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient days" "Institutional birth") plain replace

*---------------------------------------------------------*
* District hospitals (DH only)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
keep if faclevel == "DH"	// facnum is 1 and needs to be dropped
drop if district == "1210 Khounkham"
replace ipo5 = . if district == "0202 Mai"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient days" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	* Raw-count mean at e=-1 among treated districts.
	quietly summarize `m' if group_0 > 0 & event == -1, meanonly
	local baseline = r(mean)

	csdid l`m' semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_dh_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_log_dh_`m'_tidy.dta") replace

	* Transform to original scale
    preserve
    use "stata/temp/fig_a_07_did_count_log_dh_`m'_tidy.dta", clear
    generate double avgunit = `baseline'
    generate str12 outcome = "l`m'"
    foreach v in estimate conf_low conf_high point_conf_low point_conf_high {
        rename `v' `v'_log
        generate double `v' = avgunit * (exp(`v'_log) - 1)
		}
    rename std_error std_error_log
    rename statistic statistic_log
    rename p_value p_value_log
    order outcome avgunit type term event_time estimate conf_low conf_high
    save "stata/temp/fig_a_07_did_count_log_dh_`m'_tidy.dta", replace

    twoway (rcap conf_low conf_high event_time) ///
        (scatter estimate event_time) if !missing(event_time), ///
        yline(0) xline(0) xlabel(-7(1)4) legend(off) ///
        title("`t'") xtitle("Time since NHI introduction, quarter") ///
        ytitle("ATT (count units)")
    graph save "stata/graph/temp/fig_a_07_did_count_log_dh_`m'.gph", replace
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_log_dh_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_log_dh_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_dh_opo5 tbl_figa7_dh_opu5 tbl_figa7_dh_ipo5 ///
	   tbl_figa7_dh_ipu5 tbl_figa7_dh_ipday tbl_figa7_dh_del ///
using "stata/output/fig_a_07_did_count_log_dh_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient days" "Institutional birth") plain replace

*---------------------------------------------------------*
* Health centers (HC only)
*---------------------------------------------------------*

import excel "$dir/data export/02_clean_8_cln_q_fl.xlsx", firstrow clear

* Drop pilot districts and PH
drop if province == "15 Xekong"
drop if inlist(district, "1008 Met", "1009 Viangkham", "1011 Mun", "1306 Nong", "1408 Samouay", "0201 Phongsali")
drop if faclevel == "PH"

* Filter on quarter and drop "01 Vientiane Capital"
keep if inrange(quarter, 2015.4, 2017.3)
keep if province != "01 Vientiane Capital"
keep if faclevel == "HC"
gen group_org = group
gen group_0 = group
recode group_0 (9=0) 

* Loop to estimate NHI effects
local varlist opo5 opu5 ipo5 ipu5 ipday del
local labels `" "Outpatient ≥5" "Outpatient <5" "Inpatient ≥5" "Inpatient <5" "Inpatient days" "Institutional birth"  "'

local vl = wordcount("`varlist'")

forvalues i = 1/`vl' {
	
    local m: word `i' of `varlist'
	local t: word `i' of `labels'	
	
	* Raw-count mean at e=-1 among treated districts.
	quietly summarize `m' if group_0 > 0 & event == -1, meanonly
	local baseline = r(mean)

	csdid l`m' facnum semp urba litr imps, ivar(dnum) time(time) ///
		gvar(group_0) method(dr) notyet base_period(universal) ///
		wboot(reps(1000) rseed(555)) cluster(pnum)
	estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
	estimates store tbl_figa7_hc_`m'
	estat tidy, saving("stata/temp/fig_a_07_did_count_log_hc_`m'_tidy.dta") replace

	* Transform to original scale
    preserve
    use "stata/temp/fig_a_07_did_count_log_hc_`m'_tidy.dta", clear
    generate double avgunit = `baseline'
    generate str12 outcome = "l`m'"
    foreach v in estimate conf_low conf_high point_conf_low point_conf_high {
        rename `v' `v'_log
        generate double `v' = avgunit * (exp(`v'_log) - 1)
		}
    rename std_error std_error_log
    rename statistic statistic_log
    rename p_value p_value_log
    order outcome avgunit type term event_time estimate conf_low conf_high
    save "stata/temp/fig_a_07_did_count_log_hc_`m'_tidy.dta", replace
	
    twoway (rcap conf_low conf_high event_time) ///
        (scatter estimate event_time) if !missing(event_time), ///
        yline(0) xline(0) xlabel(-7(1)4) legend(off) ///
        title("`t'") xtitle("Time since NHI introduction, quarter") ///
        ytitle("ATT (count units)")
    graph save "stata/graph/temp/fig_a_07_did_count_log_hc_`m'.gph", replace
	restore
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_07_did_count_log_hc_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_07_did_count_log_hc_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa7_hc_opo5 tbl_figa7_hc_opu5 tbl_figa7_hc_ipo5 ///
	   tbl_figa7_hc_ipu5 tbl_figa7_hc_ipday tbl_figa7_hc_del ///
using "stata/output/fig_a_07_did_count_log_hc_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient days" "Institutional birth") plain replace

* Remove legend
local specs all dhhc dh hc
local outcomes opo5 opu5 ipo5 ipu5 ipday del
foreach s of local specs {
    foreach m of local outcomes {
    graph use "stata/graph/temp/fig_a_07_did_count_log_`s'_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_07_did_count_log_`s'_`m'.gph", replace
	}
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_07_did_count_log_all_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_opo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_all_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_opu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_all_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_ipo5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_all_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_ipu5.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_all_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_ipday.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_all_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dhhc_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_dh_del.gph" ///
"stata/graph/temp/fig_a_07_did_count_log_hc_del.gph" ///
, rows(6) cols(4) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) xsize(6) ysize(8)
graph export "stata/graph/fig_a_07_did_count_log.pdf", replace


*==============================================================================*
* Appendix Figure A8: DiD analysis using monthly data
*==============================================================================*

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

*---------------------------------------------------------*
* Base model
*---------------------------------------------------------*

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
	method(dr) notyet base_period(universal) wboot(reps(1000) rseed(555)) ///
	cluster(pnum)	// sci applied
estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa8_base_`m'
estat tidy, saving("stata/temp/fig_a_08_did_month_base_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_08_did_month_base_`m'.gph", replace

*---------------------------------------------------------*
* Balanced NHI exposure (e=0-2)
*---------------------------------------------------------*

* Equation 3.11 from Callaway and Sant'Anna (2021) is implementable in R but not available through 'xthdidregress' or 'csdid'.

*---------------------------------------------------------*
* Last-treated units as control
*---------------------------------------------------------*

csdid `m' facnum semp urba litr imps, ivar(dnum) time(time) gvar(group_0) ///
	method(dr) nevertreated base_period(universal) wboot(reps(1000) rseed(555)) ///
	cluster(pnum)
estat event, dropmissing post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa8_last_`m'
estat tidy, saving("stata/temp/fig_a_08_did_month_last_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_08_did_month_last_`m'.gph", replace

*---------------------------------------------------------*
* No covariate (unconditional Parallel Trend Assumption)
*---------------------------------------------------------*

csdid `m', ivar(dnum) time(time) gvar(group_0) method(dr) notyet ///
	base_period(universal) wboot(reps(1000) rseed(555)) cluster(pnum)
estat event, post	// Event-study is simultaneous CI, post_avg is pointwise CI
estimates store tbl_figa8_nocov_`m'
estat tidy, saving("stata/temp/fig_a_08_did_month_nocov_`m'_tidy.dta") replace

csdid_plot
graph save "stata/graph/temp/fig_a_08_did_month_nocov_`m'.gph", replace
}

* Export NHI estimates (bootstrap simultaneous CI)
preserve
local writeopt replace
foreach m of local varlist {
    use "stata/temp/fig_a_08_did_month_base_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_08_did_month_base_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
		
    use "stata/temp/fig_a_08_did_month_last_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_08_did_month_last_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
		
    use "stata/temp/fig_a_08_did_month_nocov_`m'_tidy.dta", clear
    export excel using "stata/output/fig_a_08_did_month_nocov_simultaneous.xlsx", ///
        sheet("`m'") firstrow(variables) `writeopt'
    local writeopt sheetreplace
	}
restore

* Export NHI estimates (pointwise CI)
esttab tbl_figa8_base_opo5r tbl_figa8_base_opu5r tbl_figa8_base_ipo5r ///
	   tbl_figa8_base_ipu5r tbl_figa8_base_iplos tbl_figa8_base_delr ///
using "stata/output/fig_a_08_did_month_base_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

esttab tbl_figa8_last_opo5r tbl_figa8_last_opu5r tbl_figa8_last_ipo5r ///
	   tbl_figa8_last_ipu5r tbl_figa8_last_iplos tbl_figa8_last_delr ///
using "stata/output/fig_a_08_did_month_last_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

esttab tbl_figa8_nocov_opo5r tbl_figa8_nocov_opu5r tbl_figa8_nocov_ipo5r ///
	   tbl_figa8_nocov_ipu5r tbl_figa8_nocov_iplos tbl_figa8_nocov_delr ///
using "stata/output/fig_a_08_did_month_nocov_pointwise.csv", ///
    cells("b(fmt(3)) ci(fmt(3))") ///
	mtitles("Outpatient >=5" "Outpatient <5" "Inpatient >=5" "Inpatient <5" "Inpatient LOS" "Institutional birth") plain replace

* Remove legend
local specs base last nocov
local outcomes opo5r opu5r ipo5r ipu5r iplos delr
foreach s of local specs {
    foreach m of local outcomes {
    graph use "stata/graph/temp/fig_a_08_did_month_`s'_`m'.gph", name(g_edit, replace)
    gr_edit .legend.draw_view.setstyle, style(no)	// hide legend
    graph save "stata/graph/temp/fig_a_08_did_month_`s'_`m'.gph", replace
	}
}

* Generate graph
graph combine ///
"stata/graph/temp/fig_a_08_did_month_base_opo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_opo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_opo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_base_opu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_opu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_opu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_base_ipo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_ipo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_ipo5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_base_ipu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_ipu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_ipu5r.gph" ///
"stata/graph/temp/fig_a_08_did_month_base_iplos.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_iplos.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_iplos.gph" ///
"stata/graph/temp/fig_a_08_did_month_base_delr.gph" ///
"stata/graph/temp/fig_a_08_did_month_last_delr.gph" ///
"stata/graph/temp/fig_a_08_did_month_nocov_delr.gph" ///
, rows(6) cols(3) iscale(0.5) graphregion(margin(0 0 0 0)) imargin(0 0 1 0) ///
 xsize(6) ysize(8)
graph export "stata/graph/fig_a_08_did_month.pdf", replace

set graphics on
