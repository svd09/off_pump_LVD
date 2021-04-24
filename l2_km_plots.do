*************************************************************
* Excellent KM cumulative failure graphs *************
**************************************************************

* do file for KM plots
* publication ready
* according to Pocock, better to graph failure function
* am going to develop template for failure function graphs

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

log using "F:\opcab\opcab_paper\long_term_results\l2_followup\kmplots", replace text

* stset the data , choose the scale that is needed
clear 

import delimited "D:\opcab\opcab_paper\long_term_results\l2_followup\l2_lt_latest.csv"

codebook, compact

* to ensure that we have 1161 patients at the t_0, we need to 
* add a small value to followup time

univar survival_days fupdays surv_d_n

gen surv_dm = surv_d_n + 1 // using the most recent followup time data 

tab1 died died_n // died = 319 ; new died_n = 447 - this is the current one that is being used for the long term analysis.

stset surv_dm, fail(died_n== 1) id(pat_id) scale(365.25) // stset the data.

stdescribe

sts graph, cumhaz tmax(10) // cumulative hazard curve for whole group.

sts list,  failure at(5 10) // results at 5, 10 years for whole group.

* make graph
* keep tmax 10 years
* change y label to %
* make risk table
* add CI

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

sts graph, failure risktable tmax(10) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci ciopts(recast(rline) lp(-.-)) legend(off)  ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" ") 

* now to plot according to OPCAB
* cumulative failure plot
* scheme sj
* colors set for each group
* change size to text small
* remove title, remove legend, change axes titles
* change y axis to 90 deg

/* results according to opcab and logrank test */

sts list, failure at(5 10) by(opcab) 

sts test opcab, logrank


sts graph, failure by(opcab) tmax(10) /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(ONCAB)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(OPCAB)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off) ciopts(recast(rline) lp(--)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" ") 



* graph according to gender and then separately for opcab and oncab

* stset for only 1 group at a time
* graph and save that graph with a name

* only for males

stset surv_dm if gender == 1, f(died_n == 1) id(pat_id) scale(365.24)



* look at stset to make sure that the # are correct

sts graph, failure by(opcab) tmax(10)  /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(ONCAB)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(OPCAB)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off) ciopts(recast(rline) lp(dash)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" Male patients: OPCAB vs ONCAB") 



* only for females


stset surv_dm if gender == 2, fail(died_n == 1) id(pat_id) scale(365.24)

sts test opcab, logrank

* look at stset to make sure that the # are correct

sts graph, failure by(opcab) tmax(10)  /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(ONCAB)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(OPCAB)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off) ciopts(recast(rline) lp(dash)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" Female patients: OPCAB vs ONCAB") 




* graph according to incomplete revascularisation
* according to opcab/oncab
* again make sure that stset is correct and for the whole group


stset surv_dm, fail(died_n== 1) id(pat_id) scale(365.25) // stset the data.


* look at stset to make sure that the # are correct

sts graph, failure by(incomplete_revascularisation) tmax(10)  /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(No)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(Yes)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off) ciopts(recast(rline) lp(dash)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" Incomplete Revascularisation") 



* now look at the impact of incomplete_revascularisation
* incomplete grafting for ONCAB 

stset surv_dm if opcab == 0, fail(died_n ==1) id(pat_id) scale(365.24)

sts list, failure at(5 10) by(incomplete_revascularisation)

sts test incomplete_revascularisation, logrank

sts graph, failure by(incomplete_revascularisation) tmax(10)  /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(No)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(Yes)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off)  ciopts(recast(rline) lp(dash)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" Incomplete Revascularisation: ONCAB") 



* now incomplete_revascularisation for OPCAB

stset surv_dm if opcab == 1, f(died_n ==1) id(pat_id) scale(365.24)

sts list, failure at(5 10) 

sts test incomplete_revascularisation, logrank 

sts graph, failure by(incomplete_revascularisation) tmax(10)  /// 
plot1(lcolor(dknavy)) ///
plot2(lcolor(maroon)) ///
risktable(, color(dknavy) size(small) group(#1) rowtitle(No)) ///
risktable(, color(maroon) size(small) group (#2) rowtitle(Yes)) ///
ylab( 0 "0%" .20 "20%" .4 "40%" .60 "60%"  .80  "80%" 1 "100%", angle(90)) /// 
ci  legend(off) ciopts(recast(rline) lp(dash)) ///
xtitle(Follow-up:Years, size(small)) ///
ytitle(Cumulative Mortality, size(small)) ///
title(" Incomplete Revascularisation: OPCAB") 




* now need to obtain the log rank test for all the group comparisons

* again need to stset the data 

stset surv_dm, fail(died_n == 1) id(pat_id) scale(365.25)

*opcab vs oncab

sts test opcab, logrank

* male vs female 

sts test gender, logrank

* incomp grafting yes vs no 

sts test incomplete_revascularisation, logrank

* according to gender
* opcab vs oncab

* for males opcab vs oncab

sts test opcab if gender == 1, logrank

* for females opcab vs oncab

sts test opcab if gender == 2, logrank

* incomp grafting for opcab and then oncab

* ONCAB

sts test incomplete_revascularisation if opcab == 0, logrank

* OPCAB

sts test incomplete_revascularisation if opcab == 1, logrank

// obtain estimates.

sts list,fail at(5 10)

sts list, fail by(opcab) at(5 10)

gen fupyears = surv_dm/365.25

univar fupyears

histogram  fupyears // skewed distribution 