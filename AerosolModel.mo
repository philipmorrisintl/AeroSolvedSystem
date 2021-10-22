within AeroSolvedSystem;


model AerosolModel

  constant Real PI=3.1415;
  constant Real k=1.3806488e-23 "Boltzmann constant";
  constant Real NA=6.02214129e23 "Avogadro number";

  parameter Integer ns=3;
  parameter Real pc_s[ns]=1e6*{6.13,5.17,4.42};
  parameter Real alpha[ns]={-8.69,-8.54,-8.42};
  parameter Real beta[ns]={1.18,1.96,2.23};
  parameter Real gamma[ns]={-4.88,-7.69,-8.25};
  parameter Real delta[ns]={-1.59,-2.95,-0.71};
  parameter Real Tc_s[ns]={513.92,536.78,563.05};
  parameter Real A[ns]={1060.6,1050.1,1050.3};
  parameter Real B[ns]={0.96,0.85,0.88};
  parameter Real a[ns]=1e-3*{24.05,25.26,27.18};
  parameter Real b[ns]=1e-3*{0.083,0.078,0.090};
  parameter Real Tref=273.15;
  parameter Real ro=1.2;
  parameter Real mu=1e5;
  parameter Real sg=1.3;
  parameter Real M_s[ns]={0.04607,0.06009,0.07412};
  parameter Real D_s[ns]={1e-5,1e-5,1e-5} "to be dep on T, (46) in Winkelmann";
 
  // INPUTS
  Real T=273.15+20;
  Real p=101325;
  Real Y_s[ns]={0.1,0.2,0.6};
  Real Z_s[ns]={0.01,0.02,0.07};
  Real N=1e14;

  // INTERNALS
  Real Tr_s[ns];
  Real tau_s[ns];
  Real psat_s[ns];
  Real rol_s[ns]; 
  Real sigma_s[ns];
  Real rol;
  Real dm;
  Real dbar;
  Real m_s[ns];
  Real W_s[ns];
  Real sigma_d; // relative to droplet
  Real v_s[ns];
  Real E_s[ns];
  Real Xs_s[ns];
  Real mg;
  Real Ys_s[ns];
  Real lambda;
  Real f;
  Real pv_s[ns];
  Real S_s[ns];
  Real Sec_s[ns];
  // Nucleation
  Real w_s[ns];
  Real H_s[ns];
  Real alpha_n;
  Real v;
  Real sigma_c;
  Real r;
  Real DG;
  Real ceq;
  Real smon_s[ns];
  Real nnuc,Z,Ntot;
  Real m,N_s[ns],Ki_s[ns];
  Real Rav,JN,Snuc_s[ns];
  Real Kl,Ks,Kbar,Jc;
  
  Real sum_of_ws(start=1);
  
equation
  for i in 1:ns loop
    Tr_s[i] = T/Tc_s[i];
    tau_s[i] = 1-Tr_s[i];
    psat_s[i] = pc_s[i]
                *exp((alpha[i]*tau_s[i]+beta[i]*tau_s[i]^1.5
                +gamma[i]*tau_s[i]^2.5+delta[i]*tau_s[i]^5)/Tr_s[i]);
    rol_s[i] = A[i]-B[i]*T;
    sigma_s[i] = a[i]-b[i]*(T-Tref);
  end for;  
  
  rol = sum(Z_s)/sum(Z_s[i]/rol_s[i] for i in 1:ns);
  dm = (6*ro*sum(Z_s)/(PI*rol*N))^(1/3);
  dbar = dm*exp(-(log(sg))^2);
  
  for i in 1:ns loop
    m_s[i] = M_s[i]/NA; // Mm in kg/mol
    W_s[i] = Z_s[i]/m_s[i]/sum(Z_s[j]/m_s[j] for j in 1:ns);
  end for;
  
  sigma_d = W_s*sigma_s;
  
  for i in 1:ns loop
    v_s[i] = m_s[i]/rol_s[i];
    E_s[i] = exp(4*sigma_d*v_s[i]/(k*T*dbar));
    Xs_s[i] = W_s[i]*E_s[i]*psat_s[i]/p;
  end for;  
  
  mg = sum(Y_s)/sum(Y_s[i]/m_s[i] for i in 1:ns);  
  
  for i in 1:ns loop
    Ys_s[i] = Xs_s[i]*m_s[i]/(Xs_s[i]*m_s[i]+(1-Xs_s[i])*mg);
  end for;
  
  lambda = sqrt(8*k*T/(PI*mg))*(4*mu/5/p);
  
  f = (1+2*lambda/dbar)/(1+5.33*(lambda/dbar)^2+3.42*lambda/dbar);
  
  for i in 1:ns loop
    pv_s[i] = ro*k*T*Y_s[i]/m_s[i];
    S_s[i] = pv_s[i]/psat_s[i];
    Sec_s[i] = 2*PI*D_s[i]*dbar*ro*Ys_s[i]*f*(E_s[i]-S_s[i]/W_s[i])*N;
  end for;
  
  // Nucleation
  
  alpha_n = 2.6e28; /* FIXME */

  sum_of_ws = sum (w_s);  
 
  v = w_s*v_s;
  sigma_c = w_s*sigma_s; // relative to cluster
  r = 2*sigma_c/(k*T*alpha_n);
  DG = 4/3*PI*r^2*sigma_c;
  
  for i in 1:ns loop
    w_s[i] = S_s[i]*exp(-alpha_n*v_s[i]);
    H_s[i] = if w_s[i]>0 then 1 else 0; // Winkelmann p. 178, 1st line
  end for;  
  
  ceq = exp(-DG/k/T)*(H_s*pv_s)/k/T;
  
  for i in 1:ns loop
    smon_s[i] = (36*PI)^(1/3)*v_s[i]^(2/3);
  end for;
  
  nnuc = sum(H_s);
  Z = (sigma_c*v^2/(4*PI^2*k*T*r^4))^(1-nnuc/2);
  Ntot = 4/3*PI*r^3/v;
  
  m = N_s*m_s;
  
  for i in 1:ns loop
    N_s[i] = Ntot*w_s[i];
    Ki_s[i] = pv_s[i]/k/T*(3/4/PI)^(1/6)*sqrt(6*k*T)*sqrt(1/m_s[i]+1/m)
              *((m_s[i]/rol_s[i])^(1/3)+(4*PI/3)^(1/3)*r)^2;
  end for;
  
  Rav = w_s*w_s/sum(w_s[i]^2/Ki_s[i] for i in 1:ns);
  
  JN = Rav*Z*ceq;
  
  for i in 1:ns loop
    Snuc_s[i] = 2*JN*N_s[i]*m_s[i];
  end for;
  
  // Coagulation
  
  Kl = 2*k*T/3/mu*(1+exp((log(sg)^2))
       +2.49*lambda/dm*(exp(2+(log(sg))^2)+exp(4*(log(sg)^2))));
  
  Ks = sqrt(3*k*T*dm/rol/(1+1/sg))
       *(exp(19/8*(log(sg)^2))+2*exp(-1/8*(log(sg)^2))+exp(-5/8*(log(sg)^2)));
  
  Kbar = 1/(1/Kl^2+1/Ks^2)^2;
  
  Jc = Kbar*N^2;

annotation(
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002),
    Documentation(info = "<html><head></head><body><div><b>Model based on:</b></div><div><ul><li> Winkelmann, C., Kuczaj, A.K., Nordlund, M. et al. Simulation of aerosol formation due to rapid cooling of multispecies vapors. J Eng Math 108, 171–196 (2018). https://doi.org/10.1007/s10665-017-9918-6 . </body></html>"));
end AerosolModel;
