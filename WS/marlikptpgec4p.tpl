DATA_SECTION
  init_int Syr;
  init_int Fyr;
  init_vector Land(Syr,Fyr);
  init_int NBiom;
  init_matrix CatDynBiom(1,NBiom,1,3);
  vector YrBiom(1,NBiom);
  vector MLBiom(1,NBiom);
  vector SDBiom(1,NBiom);
INITIALIZATION_SECTION
  logB0 10.6;
  logK 10.6;
  logp 0.693;
  logr -0.693;
PARAMETER_SECTION
  //init_bounded_number logB0(9.0,11.5,1);
  init_number logB0(1);
  init_bounded_number logK(9.5,11.5,1);
  //init_number logK(1);
  init_bounded_number logp(0.6,0.8,3);
  //init_number logp(3);
  init_bounded_number logr(-0.8,0.8,2);
  //init_number logr(2);
  vector Biom(Syr,Fyr+1);
  vector PredBiom(1,NBiom);
  vector SquDiff(1,NBiom);
  number Prod;
  sdreport_number B0;
  sdreport_number K;
  sdreport_number p;
  sdreport_number r;
  sdreport_vector Biom_sd(Syr,Fyr+1);
  objective_function_value ff;
PRELIMINARY_CALCS_SECTION
  YrBiom=column(CatDynBiom,1);
  MLBiom=column(CatDynBiom,2);
  SDBiom=column(CatDynBiom,3);
PROCEDURE_SECTION
  B0=mfexp(logB0);
  K=mfexp(logK);
  p=mfexp(logp);
  r=mfexp(logr);
  int yr;
  int yrBiom;
  Biom(Syr)=B0;
  for(yr=Syr+1;yr<=Fyr+1;yr++)
    {
     Prod=Biom(yr-1)+r*Biom(yr-1)*(1-pow(Biom(yr-1)/K,p-1));
     if(Prod<=Land(yr-1)) Prod=Land(yr-1)+1000;
     Biom(yr)=Prod-Land(yr-1);
    }
  dvariable loglikBiom;
  loglikBiom=0; 
  yrBiom=1;
  for(yr=Syr;yr<=Fyr;yr++)
    {
      if(yr==YrBiom(yrBiom))
        {
         PredBiom(yrBiom)=(Biom(yr)+Biom(yr+1))/2.;
         SquDiff(yrBiom)=square(MLBiom(yrBiom)-PredBiom(yrBiom));
         loglikBiom = loglikBiom - 0.5*(log(2*3.1416*square(SDBiom(yrBiom)))+SquDiff(yrBiom)/square(SDBiom(yrBiom)));
         yrBiom=yrBiom+1;
        }
    }
  //cout<<"loglikBiom"<<loglikBiom<<endl;
  ff = -(loglikBiom);
  Biom_sd=Biom;
