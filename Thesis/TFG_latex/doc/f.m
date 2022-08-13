function y = f(t,x)

% INTRINSIC PARAMETERS

    % Thalamic cells
    
    C_th = 1;
    I_th = 5;  
    Per_th = 50; 
    Dur_th = 5; 
    Shift_th = -80;
    Gna_th = 3;   
    Gk_th = 5;  
    Gl_th = 0.05; 
    Ena_th = 50;
    Ek_th = -90;
    El_th = -70; 
    Gt_th = 5;
    Et_th = 0;
    qht = 5.5;
    tadj = 1;
    apr = 4;
    apt = 0.3;
        
    % STN cells
    
    Cm_sn = 1;
    El_sn = -60; 
    Ena_sn = 55; 
    Ek_sn = -80; 
    thetam = 30;
    sm = 15;
    Gl_sn = 2.25; 
    Gna_sn = 37.5; 
    Gk_sn = 45; 
    Gahp_sn = 9;
    Gca_sn = 0.5; 
    Eca_sn = 140; 
    k1 = 15;
    eps = 5e-5;
    kca = 22.5;
    thetas = 39;
    ss = 8;
    thetah = -39;
    sh = 3.1;
    thetan = -32;
    sn = -8;
    taun0 = 1;
    taun1 = 100;
    thn = 80;
    sigman = 26;
    tauh0 = 1;
    tauh1 = 500;
    thh = 57;
    sigmah = 3;
    phi = 0.75;
    thetat = -63;
    kt = -7.8;
    Gt_sn = 0.5; 
    phir = 0.5;
    thetar = -67;
    kr = 2;
    taur0 = 7.1;
    taur1 = 17.5;
    thr = -68;
    sigmar = 2.2;
    alpha = 5;
    beta = 1;
    ab = -30;
    rth = 0.25;
    rsig = -0.07;
    Iapp_sn = 25; 
    alphai = 1;
    betai = 0.05;
        
    % GP cells
    
    Cm_g = 1; 
    Gna_g = 120; 
    Gk_g = 30;
    Gahp_g = 30; 
    Gt_g = 0.5; 
    Gca_g = 0.1;
    Gl_g = 0.1; 
    Ena_g = 55; 
    Ek_g = -80; 
    Eca_g = 120; 
    El_g = -55; 
    thetasg = -57;
    ksg = 2;
    thetas1g = -35;
    ks1g = 2;
    thetarg = -70;
    krg = -2;
    taurg = 30;
    thetamg = -37;
    sigmamg = 10;
    thetang = -50;
    sigmang = 14;
    taun0g = 0.05;
    taun1g = 0.27;
    thng = -40;
    sng = -12;
    thetahg = -58;
    sigmahg = -12;
    tauh0g = 0.05;
    tauh1g = 0.27;
    thhg = -40;
    shg = -12;
    k1g = 30;
    kcag = 15;
    epsg = 0.0001;
    phig = 1;
    deltang = 0.1;
    deltahg = 0.05;
    alphag = 2;
    abg = -20;
    
    % GPe
    
    Iapp_ge = -2.2; 
    betag = 0.04;
    
    % GPi
    
    Iapp_gi = 5;
    betagi = 0.08;
    kcagi = 15;
    alphaggi = 1;
    betaggi = 0.1;
      
 % CONNECTIVITY PARAMETERS
 
    % Thamalic cells
    
     Gsyn_gith = 0.06;  
     Esyn_gith = -85;   
     
     % GPi cells
     
     Gsyn_gegi = 1; 
     Esyn_gegi = -100; 
     
     Gsyn_sngi = 0.3; 
     Esyn_sngi = 0; 
     
     % STN cells
     
     Gsyn_gesn = 0.9; 
     Esyn_gesn = -100; 
     
     % GPe cells
     
     Gsyn_snge = 0.3; 
     Esyn_snge = 0; 
     
     Gsyn_gege = 0;  
     Esyn_gege = -80; 


% FUNCTIONS

    % Thalamic cells
    
    minf_th=@(V) 1./(1+exp(-(V+37)/7)); 
    hinf_th=@(V) 1./(1+exp((V+41)/4)); 
    rinf_th=@(V) 1./(1+exp((V+84)/4)); 
    pinf_th=@(V) 1./(1+exp(-(V+60)/6.2)); 
    tauh_th=@(V) 1./(0.128*exp(-(V+46)/18)+apr./(1+exp(-(V+23)/5))); 
    taur_th=@(V) (28+apt*exp(-(V+25)/10.5));
    
    Il_th=@(V) Gl_th*(V-El_th);  
    Ina_th=@(V,H) Gna_th*minf_th(V).^3.*H.*(V-Ena_th); 
    Ik_th=@(V,H) Gk_th*(0.75*(1-H)).^4.*(V-Ek_th); 
    It_th=@(V,R) Gt_th*pinf_th(V).^2.*R.*(V-Et_th); 
     
    % Thalamic sensorimotor input 
    
    hv=@(x) 1./(1+exp(-x/0.001));
    
    Finput_th=@(t,I_th,Per_th,Dur_th,Shift_th) I_th*hv(sin(2*pi*(t+Shift_th)/Per_th)).*(1-hv(sin(2*pi*(t+Shift_th+Dur_th)/Per_th)));   
    
    % STN cells
    
    ninf_sn=@(V) 1./(1+exp((V-thetan)/sn)); 
    minf_sn=@(V) 1./(1+exp(-(V+thetam)/sm)); 
    hinf_sn=@(V) 1./(1+exp((V-thetah)/sh)); 
    ainf_sn=@(V) 1./(1+exp((V-thetat)/kt));
    rinf_sn=@(V) 1./(1+exp((V-thetar)/kr)); 
    sinf_sn=@(V) 1./(1+exp(-(V+thetas)/ss)); 
    binf_sn=@(R) 1./(1+exp((R-rth)/rsig)) - 1./(1+exp(-rth/rsig));
    taun_sn=@(V) taun0 + taun1./(1+exp((V+thn)/sigman)); 
    tauh_sn=@(V) tauh0 + tauh1./(1+exp((V+thh)/sigmah)); 
    taur_sn=@(V) taur0 + taur1./(1+exp((V+thr)/sigmar)); 
    
    Il_sn=@(V) Gl_sn*(V-El_sn); 
    Ik_sn=@(V,N) Gk_sn*N.^4.*(V-Ek_sn); 
    Ina_sn=@(V,H) Gna_sn*minf_sn(V).^3.*H.*(V-Ena_sn); 
    It_sn=@(V,R) Gt_sn*ainf_sn(V).^3.*binf_sn(R).^2.*(V-Eca_sn); 
    Ica_sn=@(V) Gca_sn*sinf_sn(V).^2.*(V-Eca_sn); 
    Iahp_sn=@(V,CA) Gahp_sn*(V-Ek_sn).*CA./(CA+k1);
    
    % GP cells
    
    ninf_g=@(V) 1./(1+exp(-(V-thetang)/sigmang));
    minf_g=@(V) 1./(1+exp(-(V-thetamg)/sigmamg)); 
    hinf_g=@(V) 1./(1+exp(-(V-thetahg)/sigmahg));
    ainf_g=@(V) 1./(1+exp(-(V-thetasg)/ksg)); 
    rinf_g=@(V) 1./(1+exp(-(V-thetarg)/krg));
    sinf_g=@(V) 1./(1+exp(-(V-thetas1g)/ks1g));
    taun_g=@(V) taun0g+taun1g./(1+exp(-(V-thng)/sng)); 
    tauh_g=@(V) tauh0g+tauh1g./(1+exp(-(V-thhg)/shg)); 
    taur_g=@(V) taurg; 

    Il_g=@(V) Gl_g*(V-El_g); 
    Ik_g=@(V,N) Gk_g*N.^4.*(V-Ek_g); 
    Ina_g=@(V,H) Gna_g*minf_g(V).^3.*H.*(V-Ena_g); 
    It_g=@(V,R) Gt_g*ainf_g(V).^3.*R.*(V-Eca_g); 
    Ica_g=@(V) Gca_g*sinf_g(V).^2.*(V-Eca_g); 
    Iahp_g=@(V,CA) Gahp_g*(V-Ek_g).*(CA./(CA+k1g));
        
 % VARIABLES

     % Thalamic

     V_th = x(1:2);
     h_th = x(3:4);
     r_th = x(5:6);

     % GPi 

     V_gi = x(7:22);
     n_gi = x(23:38);
     h_gi = x(39:54);
     r_gi = x(55:70);
     ca_gi = x(71:86);
     s_gi = x(87:102);

     % GPe

     V_ge = x(103:118);
     n_ge = x(119:134);
     h_ge = x(135:150);
     r_ge = x(151:166);
     ca_ge = x(167:182);
     s_gegi = x(183:198);
     s_gesn = x(199:214);

     % STN 

     V_sn = x(215:230);
     n_sn = x(231:246);
     h_sn = x(247:262);
     r_sn = x(263:278);
     ca_sn = x(279:294);
     s_sngi = x(295:310);
     s_snge  = x(311:326);


 %CONNECTIVITY MATRICES

 Msyn_gith = [1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0; 
              0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1];

 ll = 0.2;

 Msyn_snge = [0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 ll;
              0 0 1 0 0 0 1 0 0 0 0 0 0 0 ll 0
              1 0 0 0 1 0 0 0 ll 0 0 0 0 0 0 0
              0 1 0 0 0 1 0 0 0 ll 0 0 0 0 0 0
              0 0 0 1 0 0 0 1 0 0 0 ll 0 0 0 0
              0 0 1 0 0 0 1 0 0 0 ll 0 0 0 0 0
              0 1 0 0 1 0 0 0 0 0 0 0 ll 0 0 0
              1 0 0 0 0 1 0 0 0 0 0 0 0 ll 0 0
              0 0 0 0 0 0 0 ll 0 0 0 1 0 0 0 1
              0 0 0 0 0 0 ll 0 0 0 1 0 0 0 1 0
              ll 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0
              0 ll 0 0 0 0 0 0 0 1 0 0 0 1 0 0
              0 0 0 ll 0 0 0 0 0 0 0 1 0 0 0 1
              0 0 ll 0 0 0 0 0 0 0 1 0 0 0 1 0
              0 0 0 0 ll 0 0 0 0 1 0 0 1 0 0 0
              0 0 0 0 0 ll 0 0 1 0 0 0 0 1 0 0];   

 Msyn_gege = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
              1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
              0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0
              1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0
              0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0
              0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0
              0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0
              0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0
              0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0
              0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1
              0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0
              0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0
              0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0
              0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0];

 Msyn_gesn = [0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0
              1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
              0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0
              0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0
              0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0
              1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
              0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0
              0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0
              0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0
              0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1
              0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0
              0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0
              0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0
              0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0];      

% Thalamic
Isyn_gith = Gsyn_gith*(V_th-Esyn_gith).*Msyn_gith*s_gi;

f1 = -(Il_th(V_th) + Ina_th(V_th,h_th) + Ik_th(V_th,h_th) + It_th(V_th,r_th))/C_th + (Finput_th(t,I_th,Per_th,Dur_th,Shift_th) - Isyn_gith)/C_th; 
f2 = tadj*(hinf_th(V_th)-h_th)./tauh_th(V_th);
f3 = qht*(rinf_th(V_th)-r_th)./taur_th(V_th);

% GPi

Isyn_gegi = Gsyn_gegi*s_gegi.*(V_gi-Esyn_gegi); 
Isyn_sngi = Gsyn_sngi*s_sngi.*(V_gi-Esyn_sngi); 

f4 = -(Il_g(V_gi) + Ik_g(V_gi,n_gi) + Ina_g(V_gi,h_gi) + It_g(V_gi,r_gi) + Ica_g(V_gi) + Iahp_g(V_gi,ca_gi))/Cm_g + (Iapp_gi - Isyn_gegi - Isyn_sngi)/Cm_g;
f5 = deltang*((ninf_g(V_gi)-n_gi)./taun_g(V_gi)); 
f6 = deltahg*((hinf_g(V_gi)-h_gi)./tauh_g(V_gi)); 
f7 = phig*(rinf_g(V_gi)-r_gi)./taur_g(V_gi); 
f8 = epsg*(-Ica_g(V_gi) - It_g(V_gi,r_gi) - kcagi*ca_gi);
f9 = alphag*(1-s_gi).*ainf_g(V_gi+abg) - betagi*s_gi; 

% GPe

Isyn_snge = Gsyn_snge*(V_ge-Esyn_snge).*Msyn_snge*s_snge;  
Isyn_gege = Gsyn_gege*(V_ge-Esyn_gege).*Msyn_gege*s_gesn; 

f10 = -(Il_g(V_ge) + Ik_g(V_ge,n_ge) + Ina_g(V_ge,h_ge) + It_g(V_ge,r_ge) + Ica_g(V_ge) + Iahp_g(V_ge,ca_ge))/Cm_g + (Iapp_ge - Isyn_snge - Isyn_gege)/Cm_g;  
f11 = deltang*((ninf_g(V_ge)-n_ge)./taun_g(V_ge));  
f12 = deltahg*((hinf_g(V_ge)-h_ge)./tauh_g(V_ge)); 
f13 = phig*(rinf_g(V_ge)-r_ge)./taur_g(V_ge); 
f14 = epsg*(-Ica_g(V_ge) - It_g(V_ge,r_ge) - kcag*ca_ge); 
f15 = alphaggi*(1-s_gegi).*ainf_g(V_ge+abg) - betaggi*s_gegi;
f16 = alphag*(1-s_gesn).*ainf_g(V_ge+abg) - betag*s_gesn; 

% STN

Isyn_gesn = Gsyn_gesn*(V_sn-Esyn_gesn).*Msyn_gesn*s_gesn; 

f17 = -(Il_sn(V_sn) + Ik_sn(V_sn,n_sn) + Ina_sn(V_sn,h_sn) + It_sn(V_sn,r_sn) + Ica_sn(V_sn) + Iahp_sn(V_sn,ca_sn))/Cm_sn + (Iapp_sn - Isyn_gesn)/Cm_sn; 
f18 = phi*((ninf_sn(V_sn)-n_sn)./taun_sn(V_sn)); 
f19 = phi*((hinf_sn(V_sn)-h_sn)./tauh_sn(V_sn)); 
f20 = phir*((rinf_sn(V_sn)-r_sn)./taur_sn(V_sn));  
f21 = phi*eps*(-Ica_sn(V_sn) - It_sn(V_sn,r_sn) - kca*ca_sn);  
f22 = alphai*(1-s_sngi).*sinf_sn(V_sn+ab) - betai*s_sngi;  
f23 = alpha*(1-s_snge).*sinf_sn(V_sn+ab) - beta*s_snge; 

% STN
y = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13; f14; f15; f16; f17; f18; f19; f20; f21; f22; f23];


