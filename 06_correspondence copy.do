//-----------------------------------------------------------------------------//
// Correspondence Test Function 
//-----------------------------------------------------------------------------//

//---------------------------------------//
// Example
//---------------------------------------//


/*** Syntax *** /

//------------------------------//
correspondence tau se0 e_th ref
//------------------------------//

tau <- Treatment effect vector (stata matrix)
se0 <- Standard error vector (stata matrix)
alpha <- significance level 
e_th <- equivalence threshold
ref <- reference study 



matrix tau = (1.652768 , 1.3433318 , 1.3652537 , .57561028 , 1.6656503)
matrix se0 = (.22021569 , .23572236 , .22588602 , .36349321 , .20146118)
correspondence tau se0 .10 .5 2
*/

//---------------------------------------//
//Helper Functions
//---------------------------------------//
capture program drop row_to_col
program define row_to_col, rclass
	args tau 
	mata: tau = st_matrix("`tau'")
	mata:row_to_col(tau)				 
	tempname T
	matrix `T' = r(Vec)
	return matrix tau = `T', copy
end 

capture program drop name_rows 
program name_rows, nclass 
	args tau ref 
	confirm matrix `tau' 
	row_to_col `tau'					
	matrix tau = r(tau)
	local cn = rowsof(tau)
	forvalues x = 1/`cn' {				
		if `x' == 1 {
			local nums `x'
		}
		else {
			local add `x'
			local nums: list nums | add
		}
	}
	local remove `ref'					 
	local nums2: list nums - remove

	// Treatment Effect Table			 
	global names 
	foreach x of local nums {
		local name   "Study `x'"
		global names `" $names "`name'" "'
	}
	// Difference, Equivalence, and Correspondence Tables 
	global names2 
	foreach x of local nums2 {
		local name2   "`x' vs. `ref'"
		global names2 `" $names2 "`name2'" "'
	}

end 

//---------------------------------------//
// Call Mata Function / Display Results
//---------------------------------------//


capture program drop correspondence
program correspondence, rclass 
	
	version 17 
	args tau se0 alpha e_th ref nothing
	
	confirm matrix `tau'
	confirm matrix `se0'
	if "`nothing'" != "" {
		display as err "`nothing' found where nothing expected"
		exit 198
	}

	name_rows `tau' `ref' 
	
	//-------------------
	// Call Mata Function 
	//-------------------
	mata: tau = st_matrix("`tau'")
	mata: se0 = st_matrix("`se0'")
	mata: wrapper(tau, se0, `alpha', `e_th', `ref') //  Main Function 
	
	tempname effect diff equiv corr 
	// Significance Tests
	matrix `effect' = r(effect)
	return matrix E = `effect', copy
	
	// Difference Tests
	matrix `diff' = r(diff)
	return matrix D = `diff', copy
	
	// Equivalence Test
	matrix `equiv' = r(equiv)
	return matrix EQ = `equiv', copy
	
	// Correspondence Tests
	matrix `corr' = r(corr)
	return matrix C = `corr', copy
	
	//-------------------
	// Display Results
	//-------------------
	matrix effect_stats     = r(effect)
	matrix difference_stats = r(diff)
	matrix equiv_stats      = r(equiv)
	matrix correspondence 	= r(corr)
	
	#delimit ; 
	matrix rownames effect_stats = 		 $names  ; 
	matrix colnames effect_stats = 		"Treatment Effect" 
										"Std Error" 
										"P" 
										"CI Lower" 
										"CI Upper" 
										"Significance" ;
								   
	matrix rownames difference_stats =  $names2  ; 
	matrix colnames difference_stats = "Difference" 
										"Std Error"
										"P" 
										"CI Lower" 
										"CI Upper" 
										"Significance" ;

	matrix rownames equiv_stats = 		 $names2  ; 
	matrix colnames equiv_stats = 		"Difference" 
										"P Below δE"
										"P Above δE" 
										"Significance" ;
	
	matrix rownames correspondence = 	 $names2  ;
	matrix colnames correspondence = 	"at `e_th'" ; 
	
	#delimit cr 
	display " "
	display " "
	display " "
	display "-----------------------------------------------------------------"
	display "---------------------  CORRESPONDENCE TEST  ---------------------"
	display "-----------------------------------------------------------------"
	display " "
	display "Reference Study = `ref'"
	display " "
	display "-------------------------------------------------------------------------------------------"
	display "**** Individual Treatment Effects ****" 
	esttab matrix(effect_stats) , title("Treatment Effects") nomtitles  
	display " "
	display "-------------------------------------------------------------------------------------------"
	display " "
	display "**** Difference Tests Between Study `ref' and Comparison Study ****" 
	display	"H0: τ1 − τ2 = 0" 
	display	"H1: τ1 − τ2 != 0" 
    esttab matrix(difference_stats) , title("Difference Tests") nomtitles  
	display " "
	display "-------------------------------------------------------------------------------------------"
	display "**** Equivalence Tests Between Study `ref' and Comparison Study **** " 
	display "H0:  τ1 − τ2 ≥ δE   | H0 : τ1 − τ2 ≤ −δE"
    display "H1 : τ1 − τ2 < δE   | H1 : τ1 − τ2 > −δE"
	display "Threshhold = `e_th'"
	esttab matrix(equiv_stats) , title("Equivalence Tests") nomtitles  
	display " "
	display "-------------------------------------------------------------------------------------------"
	display "**** Correspondence Tests Between Study `ref' and Comparison Study ****" 
	display " "
	esttab matrix(correspondence) , title("Correspondence Tests") nomtitles  
	#delimit ;
	display  
	   "---------------------------"   _newline 
		  "Key"    					   _newline 
	   "---------------------------"   _newline 
	   " 1 = Trivial Difference  "     _newline 
	   " 2 = Difference          "     _newline 
	   " 3 = Equivalence         "     _newline 
	   " 4 = Indeterminacy       "     _newline 
	   "---------------------------" ;
	display 
	"Trivial Difference (1) if DT(αR) = 1 & ET(δE, αR) = 1"    _newline
	"Difference (2) if DT(αR) = 1 & ET(δE, αR) = 0"            _newline
	"Equivalence (3) if DT(αR) = 0 & ET(δE, αR) = 1"           _newline
	"Indeterminacy (4) if DT(αR) = 0 & ET(δE, αR) = 0"         _newline ; 
	
	#delimit cr
	
end 

//---------------------------------------//
// Mata Functions 
//---------------------------------------//

version 17 
set matastrict on 
mata 
mata clear 

//---------------------------------------//
// Helper Functions
//---------------------------------------//


// Check column vector 
function is_colvec(z) return(anyof(("colvector","scalar"),orgtype(z)))

// Transpose to colvector 
void row_to_col(real vector v) {
	
	real matrix 			Vec
	
	if (is_colvec(v) == 0) {
		v = v'
	}
	else v = v 
	
	Vec = v
	st_matrix("r(Vec)", Vec)
}

//---------------------------------------//
// Define Structures 
//---------------------------------------//

// INPUTS
struct inputs 
{
	//----------------------- Input Variables ----------
	real scalar ref 					// reference Study 
	real scalar alpha 					// alpha level for significance tests
	real scalar e_th 					// std threshold for equiv test 
	real vector tau 					// vector of treatment effects
	real vector se0 					// vector of standard errors 
	//----------------------- Derived Variables  ------- 
	real scalar c      					// count of studies 
	real scalar z_crit 					// z_critical value 
	real vector tau1   					// reference treatment effect vector 
	real vector tau2   					// comparison treatment effect vector 
	real vector se1    					// reference se vector 
	real vector se2    					// comparison se vector 
	real vector se3    					// Pooled (pairwise) se vector 
	real vector delta  					// difference 	
}

// OUTPUTS 
struct effect_outputs
{
	real matrix effect       			// treatment effects & sig. tests
	real vector sig        				// significance indicator vector
}
struct diff_outputs
{
	real matrix diff 					// difference tests
	real vector sig_delta     			// difference tests sig. vector 
}
struct equiv_outputs
{
	real matrix equiv					// equivalence tests
	real vector sig_eq                  // equivalence tests sig. vector 
}
struct corr_outputs
{
	real vector corr 					// correspondence tests string vector 
	real vector corr_num			 	// correspondence tests numeric vector 
}
struct outputs 							// Contains sub-structures 
{
	struct effect_outputs scalar effect_o 
	struct diff_outputs scalar diff_o
	struct equiv_outputs scalar equiv_o
	struct corr_outputs scalar corr_o
}

//---------------------------------------//
// Call inputs Structure Function 
//---------------------------------------//

struct inputs scalar corr_inputs(real matrix tau, 
					             real matrix se0, 
					             real scalar alpha,
					             real scalar e_th, 
							     real scalar ref )
{

	//-------------------------------------------------------//
	struct inputs  				 	scalar in
	real scalar 				 	i , a , b 
	//----------------------- Input Variables ----------
	
	row_to_col(tau) 
	row_to_col(se0) 

	in.tau = tau
	in.se0 = se0
	in.ref = ref 
	in.alpha = alpha 
	in.e_th = e_th 
	
	//----------------------- Derived Variables  ------- 
	
	in.z_crit = abs(invnormal(.05/2))	
	in.c = rows(tau)
	in.tau1 = J(1,in.c-1,tau[ref,1])' 
	in.se1 = J(1,in.c-1,se0[ref,1])'	 
	in.tau2 = J(0,1,.) 
	in.se2 = J(0,1,.) 
	for (i=1;i<=in.c; i++) {
		if (i == ref) continue 
		else {
			a = tau[i,1]
			b = se0[i,1]
			in.tau2 = (in.tau2\a)
			in.se2 = (in.se2\b)
		}
	}
	in.se3 = sqrt(in.se1:^2 + in.se2:^2)   
	in.delta = abs(in.tau1 - in.tau2) 
	
	//-------------------------------------------------------//
	
	return(in)	
}
//---------------------------------------//
// Study Treatment Effects Function 
//---------------------------------------//

struct effect_outputs scalar effect_stat(struct inputs scalar x)
{
	//-------------------------------------------------------//
	struct effect_outputs           scalar r
	real scalar 					i, a 
	real vector 					z, p, u_tau, l_tau
	real matrix 					ci_tau
	r.effect = r.sig = 				.m
	//-------------------------------------------------------//

	// individual significance tests // 
	z = x.tau :/ x.se0
	p = J(0,1,.)
	for (i=1;i<=length(z); i++) {	
		a = (1 - normal(z[i,1])) 
		p = (p\a)
	}
	r.sig = p :<= x.alpha 
	
	// ci for individual studies 
	u_tau = x.tau +  x.z_crit:*x.se0 
	l_tau = x.tau + -x.z_crit:*x.se0 
	ci_tau = (l_tau,u_tau)
	
	r.effect = (x.tau,x.se0,p,ci_tau,r.sig) 
	return(r)
}

//---------------------------------------//
// Difference Test Function 
//---------------------------------------//
struct diff_outputs scalar diff_test(struct inputs scalar x)
{
	
	//-------------------------------------------------------//
	struct diff_outputs 			scalar r
	real scalar 					i, a 
	real vector 					z_delta, p_delta, u_delta, l_delta
	real matrix 					ci_delta
	r.diff = r.sig_delta = 			.m
	//-------------------------------------------------------//
	
	x.se3 = sqrt(x.se1:^2 + x.se2:^2) 
	x.delta = abs(x.tau1 - x.tau2) 
	z_delta = x.delta :/ x.se3  
	p_delta = J(0,1,.)
	
	for (i=1;i<=length(z_delta); i++) {	
		a = (1 - normal(abs(z_delta[i,1])))  
		p_delta = (p_delta\a)
	}
	
	r.sig_delta = p_delta :<= x.alpha / 2 
	
	// ci for the difference 
	u_delta  = x.delta + x.z_crit:*x.se3  
	l_delta = x.delta + -x.z_crit:*x.se3 
	ci_delta = (l_delta,u_delta)

	r.diff = (x.delta,x.se3,p_delta:*2,ci_delta,r.sig_delta)  
	return(r)
	
}
//---------------------------------------//
// Equivalence Test Function 
//---------------------------------------//
struct equiv_outputs scalar equiv_test(struct inputs scalar x)
{
	//-------------------------------------------------------//
	struct equiv_outputs 			scalar r
	real scalar 					i, a 
	real vector 					l_z_eq, u_z_eq, p_eqp, p_eqn 
	real matrix 					p_eq
	r.equiv = r.sig_eq = 			.m 
	//-------------------------------------------------------//
	
	x.se3 = sqrt(x.se1:^2 + x.se2:^2)  
	x.delta = abs(x.tau1 - x.tau2) 
	l_z_eq = (x.delta :+ -x.e_th):/x.se3 
	u_z_eq  = (x.delta :+ x.e_th):/x.se3 
	
	// pval below e_th
	p_eqp = J(0,1,.)
	for (i=1;i<=length(l_z_eq); i++) {	
		a = (normal(l_z_eq[i,1])) 
		p_eqp = (p_eqp\a)
	}
	
	// pval above e_th
	p_eqn = J(0,1,.)
	for (i=1;i<=length(u_z_eq); i++) {	
		a = ((1 - normal(u_z_eq[i,1]))) 
		p_eqn = (p_eqn\a)
	}
	
	// combine and check significance
	p_eq = (p_eqp, p_eqn)
	r.sig_eq = J(0,1,.)
	for (i=1;i<=rows(p_eq); i++) {  
		if (p_eq[i,1] <= x.alpha  & p_eq[i,2] <= x.alpha) {
			a = 1 
		}
		else a = 0 
		r.sig_eq = (r.sig_eq\a)
	}
	r.equiv = (x.delta, p_eq, r.sig_eq)
	return(r)
}
//---------------------------------------//
// Correspondence Test 
//---------------------------------------//
struct corr_outputs scalar corr_test(struct diff_outputs scalar x1,
									 struct equiv_outputs scalar x2)								
{
	//-------------------------------------------------------//
	struct corr_outputs 			scalar r
	real scalar 					i
	transmorphic matrix 			a  
	real matrix 					sig
	r.corr = r.corr_num = 			.m 
	//-------------------------------------------------------//
	
	sig = (x1.sig_delta, x2.sig_eq)
	
	r.corr = J(0,1,"")
	for (i=1;i<=rows(sig); i++) {	
		if (sig[i,1] == 1 & sig[i,2] == 1 ) {
		 a = "Trivial Difference"
				} 
		if (sig[i,1] == 1 & sig[i,2] == 0 ) {
		 a = "Difference"  
				} 
		if (sig[i,1] == 0 & sig[i,2] == 1) {
		 a = "Equivalence"   
			   } 
		if ( sig[i,1] == 0 & sig[i,2] == 0) {
           a = "Indeterminacy"     
				} 
		r.corr = (r.corr\a)
	}
	
	r.corr_num = J(0,1,.)
	for (i=1;i<=rows(sig); i++) {	
		if (sig[i,1] == 1 & sig[i,2] == 1 ) {
		 a = 1
                } 
		if (sig[i,1] == 1 & sig[i,2] == 0 ) {
		 a = 2
                } 
		if (sig[i,1] == 0 & sig[i,2] == 1) {
		 a = 3
                } 
		if ( sig[i,1] == 0 & sig[i,2] == 0) {
		 a = 4
                } 
		r.corr_num = (r.corr_num\a)
	}
	return(r)
}

void wrapper(real matrix tau, 
			 real matrix se0, 
		     real scalar alpha,
			 real scalar e_th,
			 real scalar ref )
{
	
	//-------------------------------------------------------//
	struct outputs 				scalar r
	struct inputs 				scalar a
	real matrix 				E, D, EQ, C
	//-------------------------------------------------------//
	
	// obtain inputs 
	a = corr_inputs(tau, se0, alpha, e_th, ref )
	
	// obtain outputs
	r.effect_o = effect_stat(a)
	r.diff_o = diff_test(a) 
	r.equiv_o = equiv_test(a)
	r.corr_o = corr_test(r.diff_o, r.equiv_o)

	E = r.effect_o.effect
	D = r.diff_o.diff
	EQ = r.equiv_o.equiv
	C = r.corr_o.corr_num

	st_rclear()
	st_matrix("r(effect)", E)
	st_matrix("r(diff)", D)
	st_matrix("r(equiv)", EQ)
	st_matrix("r(corr)", C)
}
end 






