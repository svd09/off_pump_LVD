* long term survival all cause mortality 
* parametric survival model 
* flexible RP model with 5 DF

macro drop _all // remove macros from previous work, if any
capture log close // Close any open logs. Capture will ignore a command that gives 
//                   an error. So if there isn't an open log, instead of giving you 
//                   an error and stopping here, it'll just move onto the next line.
clear all // clean the belfries
drop _all // get rid of everything!

version 15 // Every version of Stata is slightly different, but all are backwards 
//            compatible with previous ones. If you open up this do file with a way 
//            newer version, it'll run it in version 14 compatibility mode. Change 
//            this to the current version of Stata that you are using. This will 
//            also keep your code from running on older versions of stata that will 
//            break with new code that it isn't designed to handle. 

set more off, permanently // so you don't have to keep clicking through stata to 
//                           keep it running

set linesize 255 // this keeps longer lines from getting clipped. Helpful for making 
//                  tables.

log using "F:\opcab\opcab_paper\long_term_results\l2_followup\mortality_psm", replace text

import delim "D:\opcab\opcab_paper\long_term_results\l2_followup\l2_lt_latest.csv"

gen fupyears = (surv_d_n + 1)/365.25

stset fupyears, fail(died_n) id(pat_id)

// choosing the parametric model that fits best 

* plot the cumulative hazard for our data.

. grstyle init

. grstyle color background white

. grstyle color major_grid gs8

. grstyle linewidth major_grid thin

. grstyle linepattern major_grid dot

. grstyle yesno draw_major_hgrid yes

. grstyle yesno grid_draw_min yes

. grstyle yesno grid_draw_max yes

. grstyle anglestyle vertical_tick horizontal

. grstyle gsize axis_title_gap tiny

sts graph, cumhaz ///
xtitle("Follow-up:Years") ///
ytitle("Cumulative Hazard") ///
legend(off)
saving("G:\opcab\opcab_paper\long_term_results\l2_followup\supp_figures\figure1.pdf"), replace 

* plot null model 

stpm2, df(5) scale(hazard)

predict f, failure

sts gen surv2 = s

gen fail = 1 - surv2

set scheme plotplain

graph twoway (line f _t if _t <= 11, sort) ///
(line fail _t if _t <= 11, sort) 


// the model with hazard scale and df 5 fits the data perfectly.


* spline for age_at_surgery

mkspline age_sp = age_at_surgery, cubic display nknots(3)
mat age_knots = r(knots)

* spline for lv function

mkspline lv_sp = preop_lvef, cubic display nknots(3)
mat lv_knots = r(knots)

tab coronary_disease_detail

gen tvd = .
replace tvd = 1 if coronary_disease_detail > 2
replace tvd = 0 if tvd == .

tab tvd

replace preop_mi = 0 if preop_mi == -9

gen bmi_opcab = bmi*opcab

// 10-25-2020 -- had to redo the model with adding gender , dont know how I forgot it earlier 

stpm2 age_sp* lv_sp* i.opcab i.incomplete_revascularisation ///
i.prior_cardiac_surgery  i.preop_mi i.diabetes bmi i.copd   i.dialysis_dependence ///
i.peripheral_vascular_disease   i.tvd  i.cohort i.gender bmi_opcab, df(5) scale(hazard) eform


predict hr, hrnum(opcab 1) hrdenom(opcab 0) timevar(fupyears) ci level(68)


graph twoway rarea hr_uci hr_lci  _t if _t < 10, sort color(dknavy%10) || ///
line hr hr_uci hr_lci _t if _t < 10, sort lcolor(black black black) lp(solid dot dot) ///
ylab(0(0.25)1.5) ///
xlab(0(2)10) ///
ytitle("Hazard Ratio for OPCAB") ///
xtitle("Follow-up Years") ///
yline(1, lcolor(red) lp(dash)) ///
legend(off)
  
  
drop hr hr_lci hr_uci


**********************
* with hazard ratios *
**********************

xi:stpm2 age_sp* lv_sp* opcab incomplete_revascularisation ///
prior_cardiac_surgery  preop_mi diabetes bmi copd   dialysis_dependence ///
peripheral_vascular_disease   tvd  i.cohort gender, df(5) scale(hazard) eform 


* plot age and lv spline 

* age 

. grstyle init

. grstyle color background white

. grstyle color major_grid gs8

. grstyle linewidth major_grid thin

. grstyle linepattern major_grid dot

. grstyle yesno draw_major_hgrid yes

. grstyle yesno grid_draw_min yes

. grstyle yesno grid_draw_max yes

. grstyle anglestyle vertical_tick horizontal

. grstyle gsize axis_title_gap tiny

xbrcspline age_sp , values(50(0.5)85) ref(60) ///
eform matknots(age_knots) gen(ctn hr lb ub) level(68)

graph twoway rarea lb ub ctn, sort color(gs14%60) || ///
line (hr lb ub ctn), lcolor(red dknavy dknavy) lp(- - -) ///
legend(off) ///
yline(1, lp(dash) lcolor(dknavy)) ///
xtitle(Age at surgery) ///
ytitle(Relative Hazard Ratio) ///
ylab(0(1)6) ///
xlab(50(5)85)

drop ctn hr lb ub

* lv function 

xbrcspline lv_sp , values(20(0.5)35) ref(30) ///
eform matknots(lv_knots) gen(ctn hr lb ub) level(68)

graph twoway rarea lb ub ctn, sort color(gs14%60) || ///
line (hr lb ub ctn), lcolor(red dknavy dknavy) lp(- - -) ///
legend(off) ///
yline(1, lp(dash) lcolor(dknavy)) ///
xtitle(LVEF %) ///
ytitle(Relative Hazard Ratio) ///
ylab(0.5(0.25)1.75) ///
xlab(20(5)35)

drop ctn hr lb ub


// At the very end of your .do file: 
log close
// Fin.

* some code; not in log file 

sts graph, cumhaz tmax(15)

sts graph, hazard kernel(gaussian) tmax(15) /* going to present this smoothed hazard  curve to specify why we chose the PSM model with flexible splines */

// using the stmp2 to graph basline hazard 

predict haz_b, hazard ci

twoway line haz_b haz_b_uci haz_b_lci _t if _t < 10, sort

predict cumhaz_b, cumhaz

twoway line cumhaz_b _t if _t < 10, sort

* using stpm2_standsurv

range timevar 0 10 1000

stpm2_standsurv, failure at1(opcab  0) at2(opcab  1) ///
ci timevar(timevar)  level(68) atvars(opcab_p oncab_p)

/* fitted nomogram for failure opcab/oncab */

graph twoway (line oncab_p opcab_p timevar, sort lcolor(blue red)) ///
(line oncab_p_uci oncab_p_lci timevar, sort color(blue blue) lp(dash dash) lw(thin thin)) ///
(line opcab_p_lci opcab_p_uci timevar, sort lcolor(red red) lp(dash dash) lwidth(thin thin)), ///
legend(off) ///
ytitle("Cumulative Mortality Rate") ///
xtitle("Follow-up Time:Years") ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) 

stpm2_standsurv if opcab == 0, failure at1(incomplete_revascularisation 0) at2(incomplete_revascularisation 1) atvars(op_inc) level(68) ci

* rename ig 

rename incomplete_revascularisation inc_g

xi:stpm2 age_sp* lv_sp* opcab inc_g ///
prior_cardiac_surgery  preop_mi diabetes bmi copd   dialysis_dependence ///
peripheral_vascular_disease   tvd  cohort gender, df(5) scale(hazard) eform 

stpm2_standsurv,  failure at1(opcab 1 inc_g 0) at2(opcab 1 inc_g 1) ///
timevar(timevar) ci level(68) atvars(opcab_no 	opcab_yes) 

stpm2_standsurv, failure at1(opcab 0 inc_g 0) at2(opcab 0 inc_g 1) ///
timevar(timevar) ci level(68) atvars(oncab_no	oncab_yes)

twoway (rarea _at1_lci _at1_uci timevar, color(blue%20)) ///
(rarea _at2_lci _at2_uci timevar, color(gray%20)) ///
(line _at1 timevar, sort lcolor(blue)) ///
(line _at2 timevar, sort lcolor(blue) lp(dash)) ///
(rarea opcab_no_lci opcab_no_uci timevar, color(red%20)) ///
(rarea opcab_yes_lci opcab_yes_uci timevar, color(gray%20)) ///
(line opcab_no timevar, sort lcolor(red)) ///
(line opcab_yes timevar, sort lcolor(red) lp(dash))



/* opcab no inc_g vs oncab inc_g yes  - fitted model*
figure for the paper */ 

twoway (rarea opcab_no_lci opcab_no_uci timevar, color(blue%20)) ///
(line opcab_no timevar, sort lcolor(blue))  ///
(rarea oncab_yes_uci oncab_yes_lci timevar, color(red%20)) ///
(line oncab_yes timevar, sort lcolor(red)), ///
legend(off) ///
ytitle("Cumulative Mortality Rate") ///
xtitle("Follow-up Time:Years") ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) ///
text(0.6 6 "ONCAB/Inc. Revasc.", color(red)) ///
text(0.2 6 "OPCAB/Comp. Revasc.", color(blue))

/* 4 lines figure */

twoway (line opcab_no timevar, sort lcolor(blue))  ///
 (line opcab_yes timevar, sort lcolor(blue) lp(dash)) ///
(line oncab_no timevar, sort lcolor(red)) ///
(line oncab_yes timevar, sort lcolor(red) lp(dash)), ///
legend(off) ///
ytitle("Cumulative Mortality Rate") ///
xtitle("Follow-up Time:Years") ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) ///
text(0.1 9 "--- Incomp Revasc") ///
text(0.4 8 "OPCAB", color(blue)) ///
text(0.7 8 "ONCAB", color(red)) 

* am going to save this new dataaset as stpm2_predict

save "F:\opcab\opcab_paper\long_term_results\l2_followup\l2_stmp2_predict.dta", replace 

// using predict to obtain the values 
// also see if contrast can be used.

margins opcab#inc_g

// Cox model and PH test

stcox age_sp* lv_sp* i.opcab i.incomplete_revascularisation ///
i.prior_cardiac_surgery  i.preop_mi i.diabetes bmi i.copd   i.dialysis_dependence ///
i.peripheral_vascular_disease   i.tvd  i.cohort i.gender

estat phtest, detail

// cumulative hazard

sts graph if _t < 10, cumhaz 

// creating model without splines


stpm2 age_at_surgery preop_lvef i.opcab i.incomplete_revascularisation ///
i.prior_cardiac_surgery  i.preop_mi i.diabetes bmi i.copd   i.dialysis_dependence ///
i.peripheral_vascular_disease   i.tvd  i.cohort i.gender bmi_opcab, df(5) scale(hazard) eform


/*

-----------------------------------------------------------------------------------------------
                               |     exp(b)   Std. Err.      z    P>|z|     [95% Conf. Interval]
-------------------------------+----------------------------------------------------------------
xb                             |
                age_at_surgery |   1.069653   .0064307    11.20   0.000     1.057123    1.082331
                    preop_lvef |     .97501   .0076285    -3.23   0.001     .9601723    .9900769
                       1.opcab |   .4868847   .3335869    -1.05   0.294     .1271244    1.864761
1.incomplete_revascularisation |   1.182175   .1310908     1.51   0.131     .9512454    1.469167
       1.prior_cardiac_surgery |    .998439   .2079219    -0.01   0.994     .6638369    1.501695
                    1.preop_mi |   1.074879   .1091768     0.71   0.477     .8808499    1.311647
                    1.diabetes |   1.325208   .1283605     2.91   0.004     1.096065    1.602256
                           bmi |   .9701618   .0147028    -2.00   0.046     .9417686    .9994111
                        1.copd |   1.459478   .2113454     2.61   0.009     1.098844     1.93847
         1.dialysis_dependence |   1.655682   .5677547     1.47   0.141     .8454522    3.242388
 1.peripheral_vascular_disease |   1.241139   .1278737     2.10   0.036     1.014196    1.518865
                         1.tvd |   .9808495   .1454396    -0.13   0.896     .7334776     1.31165
                               |
                        cohort |
                            2  |   1.130918   .1442782     0.96   0.335     .8807207    1.452193
                            3  |   1.185983   .1553091     1.30   0.193     .9175088    1.533016
                               |
                      2.gender |   .8054193   .1113479    -1.57   0.118     .6142495    1.056086
                     bmi_opcab |     1.0191   .0255651     0.75   0.451     .9702054    1.070459
                         _rcs1 |   2.933016   .1332714    23.68   0.000     2.683102    3.206208
                         _rcs2 |    .964317   .0320556    -1.09   0.274     .9034922    1.029237
                         _rcs3 |   .8055437   .0176252    -9.88   0.000     .7717292    .8408399
                         _rcs4 |   .9533476   .0157105    -2.90   0.004     .9230475    .9846423
                         _rcs5 |   1.014574   .0145099     1.01   0.312     .9865302    1.043415
                         _cons |   .0082787   .0053362    -7.44   0.000     .0023405    .0292831
------------------------------------------------------------------------------------------------
Note: Estimates are transformed only in the first equation.

*/

// create variable tar


gen tar = .

replace tar = 1 if bima == 1 & veins == 0
replace tar = 0 if tar == .

tab tar
