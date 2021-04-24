/* some code to complete the tables and other stuff for the OPCAB paper */

import delim "G:\opcab\opcab_paper\long_term_results\l2_followup\l2_lt_latest.csv"

univar age_at_surgery

tab gender

tab1 diabetes dialysis_dependence systemic_hypertension smoker hyperlipidemia copd peripheral_vascular_disease, freq

tab1 preop_mi preop_cva prior_cardiac_surgery, freq

univar preop_lvef, by (opcab)

median preop_lvef, by(opcab)

tab opcab systemic_hypertension , chi2 cell row

univar logistic_euroscore

univar logistic_euroscore, by(opcab)

median logistic_euroscore,  by(opcab)

univar length_of_surgery, by(opcab)

median length_of_surgery, by(opcab)

tab opcab coronary_disease_detail 

univar no_of_revascularizedvessels, by(opcab)

ttest no_of_revascularizedvessels, by(opcab)

tab conversion_to_onpump if opcab == 1

tab conversion_to_onpump

tab opcab lima, chi2 row cell

tab opcab bima, chi2 row cell

tab opcab radialis, chi2 row cell


tab1 incomplete_revascularisation

tab1 lima bima radialis 

tab opcab bima, chi2 row cell

* postoperative outcomes 

tab1 died_30days postop_iabp postop_ecmo postop_resuscitation re_intubation  postop_new_dialysis

tab opcab postop_low_cardiac_output, chi2 row cell

univar  hospitallenthof_stay

foreach var of varlist postop_iabp postop_ecmo postop_resuscitation re_intubation {
	tab opcab `var', chi2 row cell
}

tab opcab postop_new_dialysis, chi2 row cell

univar ekrbcs hospitallenthof_stay, by(opcab)

ttest ekrbcs, by(opcab)

ttest hospitallenthof_stay, by(opcab)

univar no_of_revascularizedvessels

// number of arterial grafts 


