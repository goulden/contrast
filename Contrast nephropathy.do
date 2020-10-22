/*Statistical code for Goulden et al., 'Effect of intravenous contrast on renal function: a quasi-experimental study using regression discontinuity'*/

//-----------------------PREPARE DATA-----------------------//

//D-DIMER
import delimited "lab_ddimer.csv", clear
rename labvalue_ddimer ddimer
rename test_uofm dd_uofm
rename ddimer_time dd_time
//keep earliest ddimer
egen firstDD = min(dd_time), by(studyid)
keep if dd_time==firstDD
drop firstDD
/*Standardize the D-dimer units into FEU (https://unitslab.com/node/83).
All but the 6% measured in ng/ml DDU are measured in increments equivalent to
10 ng/ml FEU, so standardize to this, with each whole number = 10 ng/ml*/
replace dd_uofm="mg/L FEU2" if dd_uofm=="mg/L"
replace dd_uofm="ng/mL DDU" if dd_uofm=="ng/mL"
replace ddimer=ddimer*2 if (strpos(strupper(dd_uofm), "DDU")) //double DDU values to FEU
replace ddimer = ddimer*100 if dd_uofm != "ng/mL DDU"
replace ddimer = ddimer/10 if dd_uofm == "ng/mL DDU"
replace ddimer = round(ddimer,1) //fix floats
drop if dd_uofm=="UG FEU/ml" | dd_uofm=="ug{FEU}/mL" //drop 2 observations as they are the only ones with that unit of measurement
//drop duplicates
sort studyid
by studyid: gen dup = cond(_N==1,0,_n)
drop if dup>1
drop dup
save ddimer, replace

//eGFR
import delimited "lab_gfr.csv", clear
rename labvalue_gfr gfr
drop test_uofm
//get long-term values i.e. latest of 7 days to 6 months
sort studyid gfr_time
collapse (first) gfr gfr_time (last) lastGFR=gfr lastGFRTime=gfr_time, by(studyid)
replace lastGFR=. if lastGFRTime <= (60*24*7)
replace lastGFRTime=. if lastGFRTime <= (60*24*7)
save gfr, replace

//CREATININE
import delimited "lab_creat.csv", clear
rename labvalue_cr creat
//drop non-standard creat measures
keep if test_uofm == "umol/L"
drop test_uofm
//get AKI as highest value 12h to 7days
gen wk1Cr = creat if cr_time >= (60*12) & cr_time <=(60*24*7)
egen wk1MaxCr = max(wk1Cr), by(studyid)
egen wk1MaxCrTime = max(cr_time) if creat==wk1MaxCr & !missing(wk1Cr), by(studyid)
//get long-term values i.e. latest of 7 days to 6 months
sort studyid cr_time
collapse (first) creat cr_time (last) lastCr=creat lastCrTime=cr_time (firstnm) wk1MaxCr wk1MaxCrTime, by(studyid)
replace lastCr=. if lastCrTime <= (60*24*7)
replace lastCrTime=. if lastCrTime <= (60*24*7)
save creat, replace

//MERGE ALL INTO COHORT
import delimited "cohort.csv", clear
save CIN, replace
merge 1:1 studyid using ddimer.dta, keep(match) nogen
merge 1:1 studyid using gfr.dta, keep(match) nogen
merge 1:1 studyid using creat.dta, keep(match) nogen 
save CIN, replace
use CIN, clear

//RECODE VARS
gen male = (sex=="M")
gen charlsonfull = charlson
recode charlson (5/max=5) //5 or more = 5, as small numbers
recode ctas (9=.) //missing
rename dialysis_post_ed dialysis
replace lastGFRTime = lastGFRTime/60/24 //convert to days
gen lastGFRTimeMths = lastGFRTime/30.5 //convert to months
replace wk1MaxCrTime=wk1MaxCrTime/60 //convert to hours

//CENTRE D-DIMER ON CUTOFFS
gen cutoff = 46 if dd_uofm=="mg/L(DDU)" | dd_uofm=="ng/mL DDU" //original cutoff was 230 ng/ml DDU, but already converted to FEU
replace cutoff = 47 if dd_uofm=="mg/L FEU" //Calgary Zone until Feb 2018
replace cutoff = 50 if dd_uofm=="mgFEU/L" | dd_uofm=="ugFEU/mL" | (dd_uofm=="mg/L FEU" & year==2018) | dd_uofm=="mg/L FEU2" | dd_uofm=="mg/L(FEU)"
gen ddimer_c = ddimer-cutoff
drop if ddimer_c==0 //exclude those exactly on the cutoff
gen ddpos = (ddimer_c>0)

//CHECK FOR THOSE NOT HAVING CREAT/GFR CHECKED ON ADMISSION
gen ddgfrtime = gfr_time-dd_time
gen ddcrtime = cr_time-dd_time
//drop if not checked within 2 hrs
drop if ddgfrtime>120
drop if ddcrtime>120

//CREAT OUTCOME VARS
gen aki = ((wk1MaxCr/creat)>=1.5 | wk1MaxCr-creat>=26) if !missing(wk1MaxCr)
gen rrt = dialysis
replace rrt=1 if kt_post_ed==1
gen admit = (disp_nm=="Admit to Other Area" | disp_nm=="Admit to Special Care Unit or OR" | disp_nm=="Transfer to Another Acute Care Facility")

//CREATE SPLINES OF CONTINUOUS EXPOSURE VARS
mkspline creatcat = creat, nknots(4) cubic displayknots
mkspline agecat = age, nknots(4) cubic displayknots
mkspline gfrcat = gfr, nknots(4) cubic displayknots

//MISSING DATA VARS
misstable sum, gen(miss_)

//LABEL VARS
la var age "Age"
la var gfr "Baseline eGFR"
la var male "Male"
la var dm "Diabetes"
la var htn "Hypertension"
la var cad "Coronary artery disease"
la var cancer "Cancer"
la var charlsonfull "Charlson comorbidity index"
la var ddpos "D-dimer in relation to cutoff"
la def ddposlab 0 "Below" 1 "Above"
la val ddpos ddposlab
la var ddimer "D-dimer (ng/ml)"
la def ctpalab 0 "No" 1 "Yes"
la val ctpa ctpalab

//save
save CIN, replace
use CIN, clear
/


//-----------------------TABLES-----------------------//

//Table 1
recode ctas (1/2=1) (4/5=4)
recode charlson (3/max=3) 
table1_mc, by(ctpa) ///
vars( ///
age contn %4.0f \ ///
male bin %4.0f \ ///
gfr contn %4.0f \ ///
dm bin %4.0f \ ///
htn bin %4.0f \ ///
cad bin %4.0f \ ///
cancer bin %4.0f \ ///
ctas cat %4.0f \ ///
charlson cat %4.0f \ ///
) ///
nospace percent onecol total(after) clear
drop pvalue N_* m_*
table1_mc_dta2docx using "Tables/Table1.docx", replace font("Calibri Light",11) datafont("Calibri Light",11) datahalign(left)
use CIN, clear

//Table 2
eststo clear
foreach v of varlist lastGFR rrt aki death{
	eststo: rdrobust `v' ddimer_c, c(0) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) bwselect(msetwo) all
	eststo: rdrobust `v' ddimer_c, c(0) fuzzy(ctpa) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) bwselect(msetwo) all
}
esttab using "Tables/Table2.csv", ci nostar replace nogaps nolines wide plain

//Table 3
matrix cutoffs = (45,60,1,1) //cutoffs for subgroups
matrix orders = (1,0,0,0) //must be included as first variable is flipped in terms of the healthier group
scalar i=0
scalar row=-2
putexcel set Tables/Table3.xlsx, sheet(Sheet1) replace
foreach v of varlist gfr age dm htn{
	scalar row = row+3
	scalar i = i+1
	scalar j = orders[1,i]
	local cutoff = cutoffs[1,i]
	
	//subgroup 1
	local k = row+j
	putexcel A`k' = `"`v'<`cutoff'"'
	rdrobust lastGFR ddimer_c if `v' < `cutoff', c(0) fuzzy(ctpa) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) bwselect(msetwo)
	local effect = round(`e(tau_cl)',0.1)
	local l95 = round(`e(ci_l_rb)',0.1)
	local u95 = round(`e(ci_r_rb)',0.1)
	putexcel B`k' = `"`effect' (`l95' to `u95')"'
	scalar diff1 = `e(tau_cl)'
	scalar se1 = `e(se_tau_rb)'
	
	//subgroup 2
	local l = row+1-j
	putexcel A`l' = `"`v'>=`cutoff'"'
	rdrobust lastGFR ddimer_c if `v' >= `cutoff', c(0) fuzzy(ctpa) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) bwselect(msetwo)
	local effect = round(`e(tau_cl)',0.1)
	local l95 = round(`e(ci_l_rb)',0.1)
	local u95 = round(`e(ci_r_rb)',0.1)
	putexcel B`l' = `"`effect' (`l95' to `u95')"'
	scalar diff2 = `e(tau_cl)'
	scalar se2 = `e(se_tau_rb)'
	
	//calculating p-value for interaction using Altman-Bland method
	scalar diff = diff1-diff2
	scalar SEdiff = sqrt(se1^2+se2^2)
	scalar z = diff/SEdiff
	local p = round(2*(1-normal(abs(z))),0.01)
	local k = row
	local l = row+1
	putexcel C`k':C`l' = `p', merge
}


//-----------------------FIGURES-----------------------//

set scheme plotplain
graph set window fontface "Calibri light"

//Figure 1
matrix fig1 = (0.5,0.1 \ 80,20 \ 100,20 \ 1,0.2 \ 0.5,0.1 \ 0.5,0.1 \ 0.5,0.1 \ 0.5,0.1 \ 5,1)
local AI = "A B C D E F G H I"
local ytitles = "Proportion Years mL/min/1.73m{superscript:2} Proportion Proportion Proportion Proportion Proportion Score"
local i=0
foreach v of varlist ctpa age gfr male dm htn cad cancer ctas{
	local i=`i'+1
	local vtext : variable label `v'
	local ymax = fig1[`i',1]
	local yinc = fig1[`i',2]
	local letter = word("`AI'", `i')
	local ytitle = word("`ytitles'", `i')
	local title = `"`letter'. "`vtext'""'
	rdplot `v' ddimer_c if ddimer_c<150, c(0) nbins(10 20) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) yscale(r(0 `ymax')) ylabel(0(`yinc')`ymax') title(`title',size(med) pos(11)) ytitle(`ytitle') scale(1) play(Fig2style) saving(Graphs/`v',replace)) kernel(triangular)
}
cd Graphs
graph combine ctpa.gph age.gph gfr.gph male.gph dm.gph htn.gph cad.gph cancer.gph ctas.gph, iscale(0.7) imargin(small) graphregion(margin(r=30)) b1("           D-dimer (ng/ml) relative to cutoff          D-dimer (ng/ml) relative to cutoff          D-dimer (ng/ml) relative to cutoff",size(2.5) justification(left) margin(small))
cd ..

//Figure 2
//local
rdplot lastGFR ddimer_c if ddimer_c>=-50 & ddimer_c<=150, c(0) nbins(10 20) binselect(es) p(1) graph_options(graphregion(color(white)) legend(off) yscale(r(0 100)) ylab(0(20)100) title("A",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Long-term eGFR (mL/min/1.73m{superscript:2})") xline(-8) xline(119) play(Fig2style) saving("Graphs/PrimaryLocal",replace)) kernel(triangular) h(8.111 117.682)
//global
rdplot lastGFR ddimer_c if ddimer_c>=-50 & ddimer_c<=150, c(0) nbins(10 20) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) title("B",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Long-term eGFR (mL/min/1.73m{superscript:2})") yscale(r(0 100)) ylab(0(20)100) play(Fig2style) saving("Graphs/PrimaryGlobal",replace)) kernel(triangular)
//combine
cd Graphs
graph combine PrimaryLocal.gph PrimaryGlobal.gph, iscale(0.9) graphregion(margin(t=30))
cd ..

//-------------------IN-TEXT SUMMARY STATS AND TESTS-------------------------------//

//exact size of discontinuity in Figure 1A
rdrobust ctpa ddimer_c, c(0) h(50 100) p(3) all

//median (IQR) time to GFR mesaurement
gen lastGFRmonths = lastGFRTime/30.5
codebook lastGFRmonths

//Bertanha and Imbens F-test of external validity of complier average causal effect
rdexo lastGFR ddimer_c ctpa, cut(0) b(1000) cov(gfrcat* agecat* male dm htn cad cancer ctas charlson)


//-----------------------SUPPLEMENTARY TABLES AND FIGURES-----------------------//

//Table S1
table1_mc, ///
vars( ///
miss_lastGFR bin %4.0f \ ///
miss_aki bin %4.0f \ ///
miss_ctas bin %4.0f \ ///
) ///
nospace percent onecol total(after) clear
drop N_* m_*
table1_mc_dta2docx using "Tables/TableS1.docx", replace font("Calibri Light",11) datafont("Calibri Light",11) datahalign(left)
use CIN, clear

//Table S2
eststo clear
foreach v of varlist miss_lastGFR miss_aki lastGFRTime wk1MaxCrTime{
	eststo: rdrobust `v' ddimer_c, c(0) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) bwselect(msetwo) all
}
esttab using "Tables/TableS2.csv", ci nostar replace nogaps nolines wide plain

//Table S3
eststo clear
foreach x in 5 10 25 50{
	eststo: rdrobust lastGFR ddimer_c, c(0) fuzzy(ctpa) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) h(`x') all
}
foreach x in 2 3 4 5{
	eststo: rdrobust lastGFR ddimer_c, c(0) fuzzy(ctpa) covs(gfrcat* agecat* male dm htn cad cancer ctas charlson) h(50 20000) p(`x') all
}
esttab using "Tables/TableS3.csv", ci nostar replace nogaps nolines wide plain

//Figure S2
rdplot miss_lastGFR ddimer_c if ddimer_c<150, c(0) nbins(10 20) binselect(es) p(4) graph_options(graphregion(color(white)) legend(off) yscale(r(0 1)) ylabel(0(0.2)1) title("A",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Proportion") scale(1) play(Fig2style) saving("Graphs/GFRmiss",replace)) kernel(triangular)
rdplot lastGFRTime ddimer_c if ddimer_c<150, c(0) nbins(10 20) binselect(es) p(4) graph_options(graphregion(color(white)) legend(off) yscale(r(0 200)) ylabel(0(40)200) title("B",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Days since ED visit") scale(1) play(Fig2style) saving("Graphs/GFRtime",replace)) kernel(triangular)
cd Graphs
graph combine GFRmiss.gph GFRtime.gph, iscale(1) graphregion(margin(t=30))
cd ..

//Figure S3
rdplot rrt ddimer_c if ddimer_c<150, c(0) nbins(3 6) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) yscale(r(0 0.01)) ylabel(0(0.002)0.01) title("A",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("RRT (proportion)") scale(1) play(Fig2style) saving("Graphs/rrt",replace)) kernel(triangular)
rdplot aki ddimer_c if ddimer_c<150, c(0) nbins(5 10) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) yscale(r(0 0.5)) ylabel(0(0.1)0.5) title("B",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("AKI (proportion)") scale(1) play(Fig2style) saving("Graphs/aki",replace)) kernel(triangular)
rdplot death ddimer_c if ddimer_c<150, c(0) nbins(5 10) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) yscale(r(0 0.5)) ylabel(0(0.1)0.5) title("C",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Death (proportion)") scale(1) play(Fig2style) saving("Graphs/death",replace)) kernel(triangular)
cd Graphs
graph combine rrt.gph aki.gph death.gph, row(1) imargin(small) iscale(0.8) graphregion(margin(t=50))
cd ..

//Figure S4
rdplot miss_aki ddimer_c if ddimer_c<150, c(0) nbins(10 20) binselect(es) p(4) graph_options(graphregion(color(white)) legend(off) yscale(r(0 1)) ylabel(0(0.2)1) title("A",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Proportion") scale(1) play(Fig2style) saving("Graphs/akimiss",replace)) kernel(triangular)
rdplot wk1MaxCrTime ddimer_c if ddimer_c<150, c(0) nbins(10 20) binselect(es) p(4) graph_options(graphregion(color(white)) legend(off) yscale(r(0 120)) ylabel(0(20)120) title("B",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Hours since ED arrival") scale(1) play(Fig2style) saving("Graphs/akitime",replace)) kernel(triangular)
cd Graphs
graph combine akimiss.gph akitime.gph, iscale(1) graphregion(margin(t=30))
cd ..

//Figure S5
rdplot lastGFR ddimer_c if ddimer_c<=150 & ctpa==0, c(0) nbins(10 20) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) title("A",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Long-term eGFR (mL/min/1.73m{superscript:2})") yscale(r(0 100)) ylab(0(20)100) play(Fig1style) saving(Graphs/untreated,replace)) kernel(triangular)
rdplot lastGFR ddimer_c if ddimer_c<=150 & ctpa==1, c(0) nbins(10 20) binselect(es) p(3) graph_options(graphregion(color(white)) legend(off) title("B",size(med) pos(11)) xtitle("D-dimer (ng/ml) relative to cutoff") ytitle("Long-term eGFR (mL/min/1.73m{superscript:2})") yscale(r(0 100)) ylab(0(20)100) play(Fig2style) saving(Graphs/treated,replace)) kernel(triangular)
cd Graphs
graph combine untreated.gph treated.gph, iscale(0.9) graphregion(margin(t=30))
cd ..
