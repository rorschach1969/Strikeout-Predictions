proc import datafile = 'S:\final\pitching.txt' out = pitching replace;
delimiter = '09'x;
PROC PRINT;
RUN;

Proc SGSCATTER;
MATRIX SO IP HR BB Balls Strikes Pitches K_9 BB_9 K_BB;
run;

proc univariate normal;
var SO;
histogram /normal(mu=est sigma=est);
run;


title "Residuals of Original Data Set";
PROC REG data=pitching;
MODEL SO = IP HR BB Balls Strikes Pitches K_9 BB_9 K_BB /vif;
plot student.*(IP HR BB Balls Strikes Pitches K_9 BB_9 K_BB);
plot student.*predicted.;
plot npp.*student.;
RUN;

data pitching;
set pitching;
drop Pitches;
run;

title "Residuals without Pitches";
PROC REG data=pitching;
MODEL SO = IP HR BB Balls Strikes K_9 BB_9 K_BB / vif;
plot student.*(IP HR BB Balls Strikes K_9 BB_9 K_BB);
plot student.*predicted.;
plot npp.*student.;
RUN;

data t1Pitching;
set pitching;
log_so = log(SO);
run;

title "Log Transformation SO";
PROC REG data=t1Pitching;
MODEL log_so = IP HR BB Balls Strikes K_9 BB_9 K_BB / vif;
plot student.*(IP HR BB Balls Strikes K_9 BB_9 K_BB);
plot student.*predicted.;
plot npp.*student.;
RUN;

data LogPitching;
set pitching;
log_so = log(SO);
log_ip = log(IP);
log_hr = log(HR);
log_bb = log(BB);
log_balls = log(Balls);
log_strikes = log(Strikes);
log_k_9 = log(K_9);
log_bb_9 = log(BB_9);
log_k_bb = log(K_BB);
run;

proc print data=LogPitching;
run;


title "Log Transformation SO";
PROC REG data=LogPitching;
MODEL log_so = IP HR BB Balls Strikes log_k_9 log_bb_9 log_k_bb /vif;
plot student.*(IP HR BB Balls Strikes log_k_9 log_bb_9 log_k_bb);
plot student.*predicted.;
plot npp.*student.;
RUN;

title "Log Transformation For all X and Y Variables";
PROC REG data=LogPitching;
MODEL log_so = log_ip log_hr log_bb log_balls log_strikes log_k_9 log_bb_9 log_k_bb / vif;
plot student.*(log_ip log_hr log_bb log_balls log_strikes log_k_9 log_bb_9 log_k_bb);
plot student.*predicted.;
plot npp.*student.;
RUN;

data sqrtPitching;
set pitching;
sqrt_so = sqrt(SO);
sqrt_ip = sqrt(IP);
sqrt_hr = sqrt(HR);
sqrt_bb = sqrt(BB);
sqrt_balls = sqrt(Balls);
sqrt_strikes = sqrt(Strikes);
sqrt_k_9 = sqrt(K_9);
sqrt_bb_9 = sqrt(BB_9);
sqrt_k_bb = sqrt(K_BB);
run;

title "Sqrt Transformation SO";
Proc SGSCATTER;
MATRIX sqrt_so IP HR BB Balls Strikes K_9 BB_9 K_BB;
run;

title "Square Root Transformation For all X and Y Variables";
PROC REG data=sqrtPitching;
MODEL sqrt_so = sqrt_ip sqrt_hr sqrt_bb sqrt_balls sqrt_strikes sqrt_k_9 sqrt_bb_9 sqrt_k_bb / vif;
plot student.*(sqrt_ip sqrt_hr sqrt_bb sqrt_balls sqrt_strikes sqrt_k_9 sqrt_bb_9 sqrt_k_bb);
plot student.*predicted.;
plot npp.*student.;
RUN;

/*Here is where we decided that something might be wrong with our code and decided to run proc reg again on the original data*/

title "Re-Run Original Data Set Residuals (Curved Predicted Residuals)";
/*This was our sign that we needed polynomial regression techniques to fix our dataset*/
PROC REG data=pitching;
MODEL SO = IP HR BB Balls Strikes K_9 BB_9 K_BB / vif;
plot student.*(IP HR BB Balls Strikes K_9 BB_9 K_BB);
plot student.*predicted.;
plot npp.*student.;
RUN;

/*Start of Polynomial Analysis*/
title "Polynomial Backward Selection";
proc glmselect data=pitching;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @2/selection = backwards (stop=adjrsq);
run;

data quadratic_backward;
set pitching;
ipbb = IP*BB;
ip_balls = IP*Balls;
bb_balls = BB*Balls;
ip_strikes = IP*Strikes;
bb_strikes = BB*Strikes;
bs = Balls*Strikes;
k9ip = (K_9)*(IP);
strikes_k9 = Strikes*K_9;
ip_bb9 = IP*BB_9;
hr_bb9 = HR*BB_9;
bb_bb9 = BB*BB_9;
k9_bb9 = K_9*BB_9;
ip_kbb = IP*K_BB;
hr_kbb = HR*K_BB;
bb_kbb = (BB)*(K_BB);
balls_kbb = Balls*K_BB;
strikes_kbb = Strikes*K_BB;
k9_kbb = K_9*K_BB;
bb9_kbb = (BB_9)*(K_BB);
run;

title "Polynomial Backward Selection Residuals";
PROC REG data = quadratic_backward;
MODEL SO = HR BB ipbb Balls ip_balls bb_balls Strikes ip_strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 hr_bb9 bb_bb9 k9_bb9 ip_kbb hr_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/ vif;
plot student.*(HR BB ipbb Balls ip_balls bb_balls Strikes ip_strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 hr_bb9 bb_bb9 k9_bb9 ip_kbb hr_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;

title "Means of variables for centering values";
proc means;
var IP HR BB Balls Strikes K_9 BB_9 K_BB;
run;

data center_values;
/*This is our dataset for centering values*/
set pitching;
ipc = IP- 107.3424476;
hrc = HR-12.2832691;
bbc = BB-34.7642276;
ballsc = Balls-627.8836115;
strikesc = Strikes-1114.70;
k9c = K_9-8.1100000;
bb9c = BB_9-3.0244929;
kbbc = K_BB-2.9804664;
run;

data quadratic_backward_center;
/*This one is good*/
set center_values;
ipbb = ipc*bbc;
ip_balls = ipc*ballsc;
bb_balls = bbc*ballsc;
ip_strikes = ipc*strikesc;
bb_strikes = bbc*strikesc;
bs = ballsc*strikesc;
k9ip = (k9c)*(ipc);
strikes_k9 = strikesc*k9c;
ip_bb9 = ipc*bb9c;
hr_bb9 = hrc*bb9c;
bb_bb9 = bbc*bb9c;
k9_bb9 = k9c*bb9c;
ip_kbb = ipc*kbbc;
hr_kbb = hrc*kbbc;
bb_kbb = (bbc)*(kbbc);
balls_kbb = ballsc*kbbc;
strikes_kbb = strikesc*kbbc;
k9_kbb = k9c*kbbc;
bb9_kbb = (bb9c)*(kbbc);
run;

title "Backward selection with Centering Values";
PROC REG data = quadratic_backward_center;
MODEL SO = HR BB ipbb Balls ip_balls bb_balls Strikes ip_strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 hr_bb9 bb_bb9 k9_bb9 ip_kbb hr_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/ selection=backward vif;
plot student.*(HR BB ipbb Balls ip_balls bb_balls Strikes ip_strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 hr_bb9 bb_bb9 k9_bb9 ip_kbb hr_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;


/*The following code is commented out to save you time processing the document. UnComment the code to run*/
/*each manual step for removing values with HIGH VIF values after backward selection.*/

/*
PROC REG data = quadratic_backward_center;
MODEL SO = BB ipbb ip_balls Strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 ip_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 ip_kbb bb_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bb_strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 ip_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bs K_9 k9ip strikes_k9 BB_9 ip_bb9 ip_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bs K_9 strikes_k9 BB_9 ip_bb9 ip_kbb balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bs K_9 strikes_k9 BB_9 ip_bb9 balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb ip_balls Strikes bs K_9 strikes_k9 BB_9 ip_bb9 balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb Strikes bs K_9 strikes_k9 BB_9 ip_bb9 balls_kbb strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

PROC REG data = quadratic_backward_center;
MODEL SO = ipbb Strikes bs K_9 strikes_k9 BB_9 ip_bb9 strikes_kbb k9_kbb bb9_kbb/vif;
RUN;
*/

title "Backward selection with Centering Values after high VIF removal";
PROC REG data = quadratic_backward_center;
/*This is the model for backward center data after high VIF values removed*/
MODEL SO = Strikes bs K_9 strikes_k9 BB_9 ip_bb9 strikes_kbb k9_kbb bb9_kbb/vif;
RUN;

/*Here we decided to check forward selection as well*/

title "Polynomial Forward Selection";
proc glmselect;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @2/selection = forward (stop=adjrsq);
run;

data quadratic_forward;
/*This one is good*/
set pitching;
hrbb = HR*BB;
hr_balls = HR*Balls;
ip_strikes = IP*Strikes;
balls_kbb = Balls*K_BB;
k9_kbb = K_9*K_BB;
strikes_k9 = Strikes*K_9;
k9ip = (K_9)*(IP);
bb_kbb = (BB)*(K_BB);
bb9_kbb = (BB_9)*(K_BB);
run;

title "Polynomial Forward Selection Residuals";
PROC REG data = quadratic_forward;
MODEL SO = IP hrbb hr_balls ip_strikes K_9 k9ip strikes_k9 bb_kbb balls_kbb k9_kbb bb9_kbb/ vif stb ;
plot student.*(IP hrbb hr_balls ip_strikes K_9 k9ip strikes_k9 bb_kbb balls_kbb k9_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;

data quadratic_forward_center;
set center_values;
hrbb = hrc*bbc;
hr_balls = hrc*ballsc;
ip_strikes = ipc*strikesc;
balls_kbb = ballsc*kbbc;
k9_kbb = k9c*kbbc;
strikes_k9 = strikesc*k9c;
k9ip = (k9c)*(ipc);
bb_kbb = (bbc)*(kbbc);
bb9_kbb = (bb9c)*(kbbc);
run;

title "Polynomial Forward w/ Centered Variables Selection Residuals";
PROC REG data = quadratic_forward;
MODEL SO = IP hrbb hr_balls ip_strikes K_9 k9ip strikes_k9 bb_kbb balls_kbb k9_kbb bb9_kbb/ vif stb ;
plot student.*(IP hrbb hr_balls ip_strikes K_9 k9ip strikes_k9 bb_kbb balls_kbb k9_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;


title "Polynomial Stepwise Selection";
proc glmselect;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @2/selection = stepwise (stop=adjrsq);
run;

data quadratic;
set pitching;
k9ip = (K_9)*(IP);
bb_kbb = (BB)*(K_BB);
bb9_kbb = (BB_9)*(K_BB);
run;

title "Polynomial Stepwise Selection with VIF";
PROC REG data = quadratic;
MODEL SO = K_9 k9ip bb_kbb bb9_kbb / vif stb ;
plot student.*(K_9 k9ip bb_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;

data quadratic_stepwise_center;
/*This one is good*/
set center_values;
k9ip = (k9c)*(ipc);
bb_kbb = (bbc)*(kbbc);
bb9_kbb = (bb9c)*(kbbc);
run;

PROC REG data = quadratic_stepwise_center;
MODEL SO = k9c k9ip bb_kbb bb9_kbb / vif stb ;
plot student.*(k9c k9ip bb_kbb bb9_kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;

/*Here is where we try different polynomial glmselect method*/

title "Polynomial @3/stepwise selection";
proc glmselect data = pitching;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @3/selection = stepwise (stop=adjr);
run;

data pitching_step;
set pitching;
ipk9 = IP*K_9;
bb_kbb = BB*K_BB;
bb_k9_kbb = BB*K_9*K_BB;
ip_bb9_kbb = IP*BB_9*K_BB;
run;

title "Checking VIF values from glmselect stepwise";
proc reg data = pitching_step;
model SO = ipk9 bb_kbb bb_k9_kbb ip_bb9_kbb/ vif;
run;

title "Removal of high VIF ipk9";
proc reg data = pitching_step;
model SO = bb_kbb bb_k9_kbb ip_bb9_kbb/vif;
run;

title "Removal of high VIF bb_kbb";
proc reg data = pitching_step;
model SO = bb_k9_kbb ip_bb9_kbb/vif;
run;

/*transformations*/
data sqrt;
set pitching_step;
sqrt_SO = sqrt(SO);
sqrt_bbk9kbb = sqrt(BB*K_9*K_BB);
sqrt_ipbb9kbb = sqrt(IP*BB_9*K_BB);
run;

title "Stepwise Polynomial Regression with SQRT transformations";
proc reg data=sqrt;
model sqrt_SO = sqrt_bbk9kbb sqrt_ipbb9kbb / vif;
plot student.*(sqrt_bbk9kbb sqrt_ipbb9kbb);
plot student.*predicted.;
plot npp.*student.;
run;

/*Forward didn’t give us any good, solid models, so we’ll ignore it. */
title "Polynomial regression with Forward Selection";
proc glmselect data = pitching;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @3/selection = forward (stop=adjr);
run;

data pitching_forward;
set pitching;
k9ip = (K_9)*(IP);
bb_kbb = (BB)*(K_BB);
bb9_kbb = (BB_9)*(K_BB);
bb_k9_kbb = (BB)*(K_9)*(K_BB);
ip_bb9_kbb = (IP)*(BB_9)*(K_BB);
ip_hr_bb = (IP)*(HR)*(BB);
hr_bb_balls = (HR)*(BB)*(Balls);
ip_strikes = (IP)*(Strikes);
ip_bb_strikes = (IP)*(BB)*(Strikes);
bb_balls_strikes = (BB)*(Balls)*(Strikes);
bb_k9 = (BB)*(K_9);
strikes_k9 = (Strikes)*(K_9);
ip_strikes_k9 = (IP)*(Strikes)*(K_9);
balls_strikes_k9 = (Balls)*(Strikes)*(K_9);
hr_bb_bb9 = HR*BB*BB_9;
ip_k9_bb9 = IP*K_9*BB_9;
bb_k9_bb9 = BB*K_9*BB_9;
balls_kbb = Balls*K_BB;
ip_balls_kbb = IP*Balls*K_BB;
hr_balls_kbb = HR*Balls*K_BB;
bb_strikes_kbb = BB*Strikes*K_BB;
ip_k9_kbb = IP*K_9*K_BB;
balls_k9_kbb = Balls*K_9*K_BB;
strikes_k9_kbb = Strikes*K_9*K_BB;
bb_bb9_kbb = BB*BB_9*K_BB;
run;

title "Checking VIF values from Forward Selection";
PROC REG data = pitching_forward;
MODEL SO = bb_balls_strikes hr_bb_bb9 bb_k9_bb9 hr_balls_kbb strikes_k9_kbb / vif;
RUN;

/*
The following are high VIF values that have been removed yet still keeping a good adjrsq
*remove k9ip;
*remove ip_strikes_k9;
*remove bb_k9;
*remove bb_kbb;
*remove ip_k9_bb9;
*remove strikes_k9;
*remove ip_bb_strikes;
*remove bb_strikes_kbb;
*remove ip_k9_kbb;
*remove balls_k9_kbb;
*remove ip_strikes;
*remove ip_bb9_kbb;
*remove ip_hr_bb;
*remove IP;
*remove bb_bb9_kbb;
*remove balls_strikes_k9;
*remove Balls;
*remove hr_bb_balls;
*remove balls_kbb;
*remove ip_balls_kbb;
*remove Strikes;
*remove bb_k9_kbb;
//*trying a log transformation; --- all transformations are ugly 
*/

data pitching_forward_log;
/*Transformation method invalid because of log(0)*/
set pitching_forward;
log_so = log(SO);
bb_balls_strikes2 = log(bb_balls_strikes);
hr_bb_bb92 = log(hr_bb_bb9);
bb_k9_bb92 = log(bb_k9_bb9);
hr_balls_kbb2 = log(hr_balls_kbb);
strikes_k9_kbb2 = log(strikes_k9_kbb);
run;

title "Checking log transformations for Forward Selection";
proc reg data = pitching_forward_log;
model log_so = bb_balls_strikes2 hr_bb_bb92 bb_k9_bb92 hr_balls_kbb2 strikes_k9_kbb2 / vif;
run;

title "Checking log transformations for Forward Selection";
proc reg data = pitching_forward_log;
model log_so = bb_balls_strikes2 bb_k9_bb92 hr_balls_kbb2 strikes_k9_kbb2 / vif;
run;

*trying a sqrt transformation;
data pitching_forward_sqrt;
set pitching_forward;
sqrt_so = sqrt(SO);
sqrt_bb_balls_strikes2 = sqrt(bb_balls_strikes);
sqrt_hr_bb_bb92 = sqrt(hr_bb_bb9);
sqrt_bb_k9_bb92 = sqrt(bb_k9_bb9);
sqrt_hr_balls_kbb2 = sqrt(hr_balls_kbb);
sqrt_strikes_k9_kbb2 = sqrt(strikes_k9_kbb);
run;

title "Residuals of Forward with Sqrt Transformation";
proc reg data = pitching_forward_sqrt;
model sqrt_so = sqrt_bb_balls_strikes2 sqrt_hr_bb_bb92 sqrt_bb_k9_bb92 sqrt_hr_balls_kbb2 sqrt_strikes_k9_kbb2 / vif;
run;

title "Residuals of Forward with Sqrt Transformation Final";
proc reg data = pitching_forward_sqrt;
model sqrt_so = sqrt_bb_balls_strikes2 sqrt_bb_k9_bb92 sqrt_hr_balls_kbb2 sqrt_strikes_k9_kbb2 / vif;
plot student.*(sqrt_bb_balls_strikes2 sqrt_bb_k9_bb92 sqrt_hr_balls_kbb2 sqrt_strikes_k9_kbb2);
plot student.*predicted.;
plot npp.*student.;
run; 

title "Polynomial with Backward Selection method";
proc glmselect data = pitching;
model SO = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @3/selection = backward (stop=adjrsq);
run;


data order3backward;
/*This is the start of backward selection*/
set pitching;
iphr = IP*HR;
ipbb = IP*BB;
hrbb = HR*BB;
ipb = IP*Balls;
iphrb = IP*HR*Balls;
bbb = BB*Balls;
ipbbb = IP*BB*Balls;
hrbbb = HR*BB*Balls;
ips = IP*Strikes;
hrs = HR*Strikes;
iphrs = IP*HR*Strikes;
bbs = BB*Strikes;
hrbbs = HR*BB*Strikes;
bs = Balls*Strikes;
ipk9 = IP*K_9;
hrk9 = HR*K_9;
iphrk9 = HR*K_9*IP;
bbk9 = BB*K_9;
ipbbk9 = IP*BB*K_9;
hrbbk9 = HR*BB*K_9;
bbsk9 = BB*Strikes*K_9;
bsk9 = Balls*Strikes*K_9;
ipbb9 = IP*BB_9;
hrbb9 = HR*BB_9;
iphrbb9 = IP*HR*BB_9;
bbbb9 = BB*BB_9;
ipbbbb9 = IP*BB*BB_9;
hrbbbb9 = HR*BB*BB_9;
bbb9 = Balls*BB_9;
hrbbb9 = HR*Balls*BB_9;
bbbbb9 = BB*Balls*BB_9;
ipsbb9 = IP*Strikes*BB_9;
hrsbb9 = HR*Strikes*BB_9;
bbsbb9 = BB*Strikes*BB_9;
ipk9bb9 = IP*K_9*BB_9;
hrk9bb9 = HR*K_9*BB_9;
bk9bb9 = Balls*K_9*BB_9;
ipkbb = IP*K_BB;
hrkbb = HR*K_BB;
iphrkbb = IP*HR*K_BB;
bbkbb = BB*K_BB;
hrbbkbb = HR*BB*K_BB;
bkbb = Balls*K_BB;
hrbkbb = HR*Balls*K_BB;
bbbkbb = BB*Balls*K_BB;
ipskbb = IP*Strikes*K_BB;
ipk9kbb = IP*K_9*K_BB;
bbk9kbb = BB*K_9*K_BB;
bk9kbb = Balls*K_9*K_BB;
ipbb9kbb = IP*BB_9*K_BB;
hrbb9kbb = HR*BB_9*K_BB;
bbbb9kbb = BB*BB_9*K_BB;
run;

/*Below is the code for removing the hgih VIF values from backward selection*/
/*to run this code, highlight the entire commented section and press ctrl +  */

/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs bbs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 iphrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot student.*(HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs bbs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 iphrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb);*/
/*plot student.*predicted.;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs bbs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot student.*(HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs bbs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb);*/
/*plot student.*predicted.;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot student.*(HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb);*/
/*plot student.*predicted.;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs ipk9 hrk9 iphrk9 bbk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs ipk9 hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb hrbb9kbb bbbb9kbb/vif;*/
/*plot npp.*student.;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*model SO = HR iphr bb ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/**/
/*proc reg data=order3backward;*/
/*/*removded bb*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 ipsbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipsbb9*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 bsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bsk9*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb ipbb9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbb9kbb*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 ipk9bb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipk9bb9*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ips hrs iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrs*/*/
/*model SO = HR iphr ipbb hrbb Balls ipb iphrb bbb ipbbb hrbbb ip iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipb*/*/
/*model SO = HR iphr ipbb hrbb Balls iphrb bbb ipbbb hrbbb ip iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbb*/*/
/*model SO = HR iphr ipbb Balls iphrb bbb ipbbb hrbbb ip iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbb*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb bkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bkbb*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs hrbbs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbbs*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs bs hrk9 iphrk9 ipbbk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbbk9*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs bs hrk9 iphrk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 ipbbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbbbb9*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs bs hrk9 iphrk9 hrbbk9 bbsk9 ipbb9 hrbb9 bbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbb9*/*/
/*model SO = HR iphr ipbb Balls iphrb ipbbb hrbbb ip iphrs bs hrk9 iphrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded iphrb*/*/
/*model SO = HR iphr ipbb Balls ipbbb hrbbb ip iphrs bs hrk9 iphrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded iphrk9*/*/
/*model SO = HR iphr ipbb Balls ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 bbsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbsbb9*/*/
/*model SO = HR iphr ipbb Balls ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 hrbbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbbb9*/*/
/*model SO = HR iphr ipbb Balls ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded Balls*/*/
/*model SO = HR iphr ipbb ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbb*/*/
/*model SO = HR iphr ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded iphr*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb iphrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded iphrkbb*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb bbkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbkbb*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrsbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrsbb9*/*/
/*/*adjr sq = 0.9994*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 hrbb9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbb9*/*/
/*/*adjr sq = 0.9994*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 hrbbk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbbk9*/*/
/*/*adjr sq = 0.9994*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb bbbb9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbbb9kbb*/*/
/*/*adjr sq = 0.9994*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb ipk9kbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipk9kbb*/*/
/*/*adjr sq = 0.9986*/*/
/*model SO = HR ipbbb hrbbb ip iphrs bs hrk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbbb*/*/
/*/*adjr sq = 0.9985*/*/
/*model SO = HR ipbbb ip iphrs bs hrk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb bbbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbbkbb*/*/
/*/*adjr sq = 0.9983*/*/
/*model SO = HR ipbbb ip iphrs bs hrk9 bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrk9*/*/
/*/*adjr sq = 0.9983*/*/
/*model SO = HR ipbbb ip iphrs bs bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bs*/*/
/*/*adjr sq = 0.9981*/*/
/*model SO = HR ipbbb ip iphrs bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 ipkbb hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipkbb*/*/
/*/*adjr sq = 0.9967*/*/
/*model SO = HR ipbbb ip iphrs bbsk9 bbbb9 hrbbbb9 bbb9 bbbbb9 hrk9bb9 bk9bb9 hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbb9*/*/
/*/*adjr sq = 0.9965*/*/
/*model SO = HR ipbbb ip iphrs bbsk9 bbbb9 hrbbbb9 bbbbb9 hrk9bb9 bk9bb9 hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbbbb9*/*/
/*/*adjr sq = 0.9965*/*/
/*model SO = HR ipbbb ip iphrs bbsk9 bbbb9 hrbbbb9 hrk9bb9 bk9bb9 hrkbb hrbkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbkbb*/*/
/*/*adjr sq = 0.9964*/*/
/*model SO = HR ipbbb ip iphrs bbsk9 bbbb9 hrbbbb9 hrk9bb9 bk9bb9 hrkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bbsk9*/*/
/*/*adjr sq = 0.9964*/*/
/*model SO = HR ipbbb ip iphrs bbbb9 hrbbbb9 hrk9bb9 bk9bb9 hrkbb ipskbb bbk9kbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bk9kbb*/*/
/*/*adjr sq = 0.9959*/*/
/*model SO = HR ipbbb ip iphrs bbbb9 hrbbbb9 hrk9bb9 bk9bb9 hrkbb ipskbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded HR*/*/
/*/*adjr sq = 0.9887*/*/
/*model SO = ipbbb ip iphrs bbbb9 hrbbbb9 hrk9bb9 bk9bb9 hrkbb ipskbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded bk9bb9*/*/
/*/*adjr sq = 0.9874*/*/
/*model SO = ipbbb ip iphrs bbbb9 hrbbbb9 hrk9bb9 hrkbb ipskbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipskbb*/*/
/*/*adjr sq = 0.9834*/*/
/*model SO = ipbbb ip iphrs bbbb9 hrbbbb9 hrk9bb9 hrkbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded ipbbb*/*/
/*/*adjr sq = 0.9822*/*/
/*model SO = ip iphrs bbbb9 hrbbbb9 hrk9bb9 hrkbb bk9kbb/vif;*/
/*RUN;*/
/**/
/*proc reg data=order3backward;*/
/*/*removded hrbbbb9*/*/
/*/*adjr sq = 0.9799*/*/
/*model SO = ip iphrs bbbb9 hrk9bb9 hrkbb bk9kbb/vif;*/
/*RUN;*/
;

title "Checking VIF values from Backward Selection";
proc reg data=order3backward;
/*removded iphrs*/
/*adjr sq = 0.9759*/
/*Fvalue = 18923.8*/
/*might need transformations because of high residual values (outside of +-3)*/
model SO = ip bbbb9 hrk9bb9 hrkbb bk9kbb/vif;
RUN;

/*transformations although might not be the best in comparison to other models*/
data t1_order3backward;
set order3backward;
log_so = log(SO);
inv_ip = 1/(IP); 
log_bbbb9 = log(bbbb9);
log_hrk9bb9 = log(hrk9bb9);
log_hrkbb = log(hrkbb);
log_bk9kbb = log(bk9kbb);
run;

title "Checking VIF values from Backward Selection";
proc reg data=t1_order3backward;
/*removded iphrs*/
/*adjr sq = 0.9759*/
/*Fvalue = 18923.8*/
/*might need transformations because of high residual values (outside of +-3)*/
model log_so = inv_ip log_bbbb9 log_hrk9bb9 log_hrkbb log_bk9kbb/vif;
plot student.*(inv_ip log_bbbb9 log_hrk9bb9 log_hrkbb log_bk9kbb);
plot student.*predicted.;
plot npp.*student.;
RUN;


/*Start of training and test set*/
proc surveyselect data=pitching out=pitching_train_forward seed = 54512 outall samprate = 0.75;
run;

data polynomial_new;
set pitching_train_forward;
if selected then train_y = SO;
run;

/*proc print data = pitching_train;*/
/*run;*/

proc glmselect data = polynomial_new ;
model train_y = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @3/selection = forward (stop=adjrsq);
run;

data forward_train_pitching;
set polynomial_new;
ip_balls = IP*Balls;
ip_balls_strikes = IP*Balls*Strikes;
ip_k9 = IP*K_9;
balls_k9 = Balls*K_9;
ip_strikes_k9 = IP*Strikes*K_9;
bb_kbb = BB*K_BB;
ip_bb_kbb = IP*BB*K_BB;
ip_balls_kbb = IP*Balls*K_BB;
bb_strikes_kbb = BB*Strikes*K_BB;
bb_k9_kbb = BB*K_9*K_BB;
bb9_kbb = BB_9*K_BB;
ip_bb9_kbb = IP*BB_9*K_BB;
balls_bb9_kbb = Balls*BB_9*K_BB;
strikes_bb9_kbb = Strikes*BB_9*K_BB;
run;

title "Removing predictors with high VIF Values";
proc reg data = forward_train_pitching;
model train_y = ip_balls ip_balls_strikes K_9 ip_k9 balls_k9 ip_strikes_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb bb9_kbb ip_bb9_kbb balls_bb9_kbb strikes_bb9_kbb/ vif;
run; 

/*please uncomment the code below to see the steps of removing each high VIF value manually*/

/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 ip_strikes_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb bb9_kbb ip_bb9_kbb balls_bb9_kbb strikes_bb9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb bb9_kbb ip_bb9_kbb balls_bb9_kbb strikes_bb9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb bb9_kbb ip_bb9_kbb strikes_bb9_kbb/ vif;*/
/*run;*/
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb bb9_kbb strikes_bb9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_strikes_kbb bb_k9_kbb strikes_bb9_kbb/ vif;*/
/*run;  */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 bb_kbb ip_bb_kbb ip_balls_kbb bb_k9_kbb strikes_bb9_kbb/ vif;*/
/*run;  */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 ip_bb_kbb ip_balls_kbb bb_k9_kbb strikes_bb9_kbb/ vif;*/
/*run;  */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls ip_balls_strikes K_9 balls_k9 ip_balls_kbb bb_k9_kbb strikes_bb9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls_strikes K_9 balls_k9 ip_balls_kbb bb_k9_kbb strikes_bb9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls_strikes K_9 balls_k9 ip_balls_kbb bb_k9_kbb/ vif;*/
/*run; */
/**/
/*title "Removing predictors with high VIF Values";*/
/*proc reg data = forward_train_pitching;*/
/*model train_y = ip_balls_strikes K_9 balls_k9 ip_balls_kbb/ vif;*/
/*run; */

title "Final step of removing high VIF values";
proc reg data = forward_train_pitching;
model train_y = K_9 balls_k9 ip_balls_kbb/ vif;
run; 

data pitching_train_forward_sqrt;
set forward_train_pitching;
sqrt_train_y = sqrt(train_y);
sqrt_k9 = sqrt(K_9);
sqrt_balls_k9 = sqrt(balls_k9);
sqrt_ip_balls_kbb = sqrt(ip_balls_kbb);
run;

title "Sqrt Transformation on Forward selection training model";
proc reg data = pitching_train_forward_sqrt;
model sqrt_train_y = sqrt_k9 sqrt_balls_k9 sqrt_ip_balls_kbb / vif;
plot student.*(sqrt_k9 sqrt_balls_k9 sqrt_ip_balls_kbb);
plot student.*predicted.;
plot npp.*student.;
run;

data pitching_train_forward_log;
set forward_train_pitching;
log_train_y = log(train_y);
log_k9 = log(K_9);
log_balls_k9 = log(balls_k9);
log_ip_balls_kbb = log(ip_balls_kbb);
run;

proc reg data = pitching_train_forward_log;
model log_train_y = log_k9 log_balls_k9 log_ip_balls_kbb/ vif;
plot student.*(log_k9 log_balls_k9 log_ip_balls_kbb);
plot student.*predicted.;
plot npp.*student.;
run;

proc reg data = pitching_train_forward_log;
model log_train_y = log_k9 log_balls_k9 log_ip_balls_kbb;
output out=outm2(where=(log_train_y=.)) p=yhat;
proc print;
run;

data log_test;
set outm2;
log_so = log(SO);
p_so = log_so - yhat;
ab_so = abs(p_so);
run;

title "Test set summary statistics for Forward Selection";
proc summary data = log_test;
var p_so ab_so;
output out = pred std(p_so)=rsme mean(ab_so)=mae;
run;

proc print;
run;

proc corr data = log_test;
var log_so yhat;
run;

title "Training set summary statistics for Forward";
proc reg data = test_transform_poly;
model sqrt_train_y = sqrt_bb9kbb sqrt_ipbb9kbb;
run;

/*This is the start of the stepwise training and testing sets*/

proc surveyselect data=pitching out= pitching_train seed=615824 samprate = 0.75 outall;
run;

data polynomial_new;
set pitching_train;
if selected then train_y = SO;
run;
proc print;
run;

proc glmselect data = polynomial_new ;
model train_y = IP|HR|BB|Balls|Strikes|K_9|BB_9|K_BB @3/selection = stepwise (stop=adjr);
run;

data test_polynomial;
set polynomial_new;
bbkbb = BB*K_BB;
ipk9 = IP*K_9;
bb9kbb = BB_9*K_BB;
ipbb9kbb = IP*BB_9*K_BB;
run;

proc reg data = test_polynomial ;
model train_y =  bbkbb ipk9 bb9kbb ipbb9kbb/vif;
run;

proc reg data = test_polynomial ;
model train_y =  bbkbb bb9kbb ipbb9kbb/vif;
run;

proc reg data = test_polynomial ;
model train_y = bb9kbb ipbb9kbb/vif;
run;

data test_transform_poly;
set test_polynomial;
sqrt_train_y = sqrt(train_y);
sqrt_bb9kbb = sqrt(bb9kbb);
sqrt_ipbb9kbb = sqrt(ipbb9kbb);
run;

proc reg data = test_transform_poly;
model sqrt_train_y = sqrt_bb9kbb sqrt_ipbb9kbb;
output out=outm1(where=(sqrt_train_y=.)) p=yhat;
proc print;
run;

data sqrt;
set outm1;
sqrt_so = sqrt(SO);
p_so = sqrt_so - yhat ;
ab_so = abs (p_so);
run;

title "Test set summary statistics for stepwise";
proc summary data = sqrt;
var p_so ab_so;
output out = pred std(p_so)=rsme mean(ab_so)=mae;
run;

proc print;
run;

proc corr data = sqrt;
var sqrt_so yhat;
run;

title "Training set summary statistics for stepwise";
proc reg data = test_transform_poly;
model sqrt_train_y = sqrt_bb9kbb sqrt_ipbb9kbb;
run;
