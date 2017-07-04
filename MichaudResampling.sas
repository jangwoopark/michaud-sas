libname cleandat "C:\SAS Data\Data";
libname result "C:\SAS Data\Output";

%let beginYear = 2007;
%let endYear = 2008;

data monthlyData(drop=date dlret dlretx);
		set cleandat.Monthlycrsp;
		rename shrout=sharesOut;
		rename prc=price;
		month=month(date);
		year=year(date);
		if missing(dlret)=0 then ret=dlret;
		if missing(dlretx)=0 then retx=dlretx;
		*if year>=&startyr-2 and year<=&endyr-2;
run;


proc import datafile="C:\SAS Data\Data\riskfree.xls" out=riskfree replace;
run;

data riskfree(drop=date dateString);
	set riskfree;
	dateString = PUT (date,6.0);
	year = substr(dateString,1,4);
	month = substr(dateString,5,2);
run;


proc sort data=monthlyData;
	by permno;
run;
%let yearlength = 5;

 
%macro outputret();
	%do yr=&beginYear %to &endYear;
	data crsp;
		set monthlyData;
		if (year >= &yr-&yearLength and year < &yr);
		*if year>=&startyr-2 and year<=&endyr-2;
	run;

	proc means data=crsp noprint;
			var permno ;
			by permno;
			output out=enoughData sum=month;
	run;

	data enoughData(keep = permno);
		set enoughData;
		if (_FREQ_ = 60);
	run;

	data top500;
		merge enoughData(in=k) crsp;
		by permno;
		if k;
	run;

	data top500;
		set top500;
		if (year+1 = &yr) and (month = 12);
		year = year + 1;
	run;

	proc sort data=	top500;
		by sharesOut price;
	run;

	data top500(keep=permno year mk_cap);
		set top500;
		by sharesOut price;
		mk_cap = sharesOut*price;
	run;


	proc sort data=crsp;
		by permno DESCENDING year descending month;
	run; 

	proc sort data=top500;
		by DESCENDING mk_cap ;
	run; 

	data top500;
		set top500;
	 	count + 1;
		by year;
		if first.year then count=1;
		if count <= 500;
	run;

	proc sort data=top500;
		by permno;
	run;

	proc sort data=crsp;
		by permno;
	run;

	data yearlyReturn(keep=permno year month ret count);
		merge top500(in=k) crsp;
		by permno;
		if k;
	run;

	proc sort data=Yearlyreturn;
		by year month;
	run;

	data Yearlyreturn(drop=month  year);
		set Yearlyreturn;
		by year month ;
		datedate = mdy(month, 1, year);
	run;

	proc sort data=Yearlyreturn;
		by datedate;
	run;


	PROC TRANSPOSE DATA=Yearlyreturn OUT=returns;
		BY datedate;
		id permno;
		var ret;
	run;

	proc sort data=returns;
		by descending datedate;
	run;

	proc corr data=returns cov out=cov_matrix noprint;
	run;


	data mean(drop = _TYPE_ _NAME_ datedate);  
		set cov_matrix;
		if _TYPE_ = "MEAN";
	run;

/**/


	data permno(keep = _NAME_);  
		set cov_matrix;
		if _NAME_ = "datedate" then delete;
		if _NAME_ = "" then delete;
		if _TYPE_ = "COV";
	run;

	data permno(keep = permnoStr);
		set permno;
		permnoStr = substr(_NAME_,2,5);
	run;

	data permno(keep = permno);
		set permno;
		permno = 1*permnoStr;
	run;

	data matrix(drop = _TYPE_ datedate _NAME_);  
		set cov_matrix;  
		if (_TYPE_ = "COV");
		if (_NAME_ = "datedate") then delete;
	run;

	proc sort data=riskfree;
		by year;
	run;

	data currentRiskFree;
		set riskfree;
		by year;
		rf = rf /100;
		if ((year>=&yr- &yearLength-1) and (year < &yr-1));
	run;

	data currentRiskFree(drop=month  year);
		set currentRiskFree;
		by year month ;
		datedate = mdy(month, 1, year);
	run;

	%do iter=1 %to 10;
		data randSelectedDates(keep=datedate rand);
			set Yearlyreturn;
			by datedate;
			rand = RANUNI(time());
			if first.datedate;
			if missing(rand) then rand=	RANUNI(time());
		run;

		proc sort data=	randSelectedDates;
			by descending rand;
		run;

		data randSelectedDates(keep=datedate);
			set randSelectedDates;
		 	count + 1;
			if first.year then count=1;
			if count <= 12;
		run;

		proc sort data=	randSelectedDates;
			by datedate;
		run;


		proc sort data =  currentRiskFree;
			by datedate;
		run;

		data randRiskFree;
			merge randSelectedDates(in=k) currentRiskFree;
			by datedate;
			if (k);
		run;

		proc means data=randRiskFree noprint;
			var rf ;
			output out=meanriskfree mean=meanReturn;
		run;

		data meanriskfree(keep=meanReturn);
			set meanriskfree;
		run;

		data meanReturn;
			merge randSelectedDates(in=k) Yearlyreturn;
			by datedate;
			if (k);
		run;

		proc sort data=meanReturn;
			by permno;
		run;

		proc means data=meanReturn noprint;
			var ret ;
			by permno;
			output out=meanReturn mean=meanReturn;
		run;

		data meanReturn(keep=meanReturn);
			set meanReturn;
		run;

		proc iml;
			reset noprint;
			use matrix;
			read All into cov_mat;
			use meanriskfree;
			read all into rfr;
			use Meanreturn;
			read all into Meanreturn;
			one = j(500,1);
			inverseSigma = ginv(cov_mat);
			tangencyWeights = inverseSigma*(Meanreturn-rfr*ONE);
			CREATE tangency FROM tangencyWeights;
			APPEND from tangencyWeights;
		quit;


		/*tangency*/

		data tangency;
			merge permno tangency;
		run;

		data tangency(keep=permno wt);
			set tangency;
			if COL1 < 0 then wt = 0;
			else wt = COL1;
		run;

		proc means data=tangency noprint;
			var wt;
			output out=tangencySum sum=allweight;
		run;

		data tangency;
			merge tangency tangencySum;
			retain sumwt;  
			if _n_ = 1 then do;
				sumwt = allweight; 
			end;
			drop _freq_ _type_ allweight;
		run;

		data tangency(keep=permno wt);
			set tangency;
			wt = wt/sumwt;
		run;


		proc sort data=tangency;
			by permno;
		run;

	
		%if &iter = 1 %then %do;
			data sampleWeight;
				set tangency;
			run; 
		%end;
		%else %do;
			data sampleWeight;
				set	sampleWeight tangency;
			run;
		%end;
	%end;
	
	proc sort data=sampleWeight;
		by permno;
	run;

	proc means data=sampleWeight noprint;
		var wt;
		by permno;
		output out=sampleWeight2 mean=mean_wt;
	run;

	data thisYearCrsp;
		set monthlydata;
		if (year = &yr);
	run;

	data tangency2(keep=permno mean_wt ret retx month year);
		merge sampleWeight2(in=k) thisYearCrsp;
		by permno;
		if (k);
	run;

	data tangency2;
		set tangency2;
		by permno year month;
		lagretx=lag(retx);
		if (first.permno) then lagretx=0;
	run;

	data tangency2;
		set tangency2;
		by permno year month;
		retain dyn_wt;
		if first.permno then do;
			dyn_wt = mean_wt;
		end;
		else dyn_wt = dyn_wt*(1+lagretx);
	run;

	proc sort data=tangency2;
		by year month;
	run;

	proc means data=tangency2 noprint;
		var ret retx;
		weight dyn_wt;
		by year month;
		output out=tangWeightedReturn mean=TR PR;
	run;

	proc sort data=tangWeightedReturn;
		by year month ;
	run;

	data  tangWeightedReturn;
		set tangWeightedReturn;
		if missing(year)=0;
	run;

	%if &yr = &beginYear %then %do;
		data tangPortfolio;
			set	tangWeightedReturn;
		run; 
	%end;
	%else %do;
		data tangPortfolio;
			set	tangPortfolio tangWeightedReturn;
		run;
	%end;

	data result.tangPortfolio;
		set tangPortfolio;
	run;


	%end;
%mend outputret;

	

%outputret()

		proc export data=result.tangPortfolio outfile="C:\SAS Data\Output\michaudBoot.xls" replace;
		run;


