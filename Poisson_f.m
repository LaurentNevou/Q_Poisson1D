function[OUTPUT,ro3DEfn,ro3DEfp,ErrVec] = Poisson_f(Structure,En,Ep,T,EfL,EfR,Nloops,tau0,k,Video_convergence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 6.62606896E-34;           %% Planck constant J.s
hbar  = h/(2*pi);
e     = 1.602176487E-19;          %% charge de l electron Coulomb
m0    = 9.10938188E-31;           %% electron mass kg
%c     = 2.99792458e8;             %% speed of light (m/s)
Epsi0 = 8.854187817620E-12;       %% constant dielectric du vide F/m
%mu0   = 1/(Epsi0*c^2);            %% permeabiliy du vide
kB    = 1.3806488E-23;            %% Boltzmann's constant (J/K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z     = Structure(:,1)';
V0    = Structure(:,2)';
Eg    = Structure(:,3)';
Dop   = Structure(:,4)';
Epsi  = Structure(:,5)';
EfXX  = Structure(:,6)';
Mass_n= Structure(:,7)';
Mass_p= Structure(:,8)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Meshgrid of density matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ro3Dn_const = (1/(2*pi^2)) * ( (2*e*Mass_n*m0/(hbar^2)).^(3/2) );
ro3Dp_const = (1/(2*pi^2)) * ( (2*e*Mass_p*m0/(hbar^2)).^(3/2) );

[ro3Dn_const_M,EEn] = meshgrid(ro3Dn_const,En); % put the vector Mass_n in a matrix En-long
[ro3Dp_const_M,EEp] = meshgrid(ro3Dp_const,Ep); % put the vector Mass_p in a matrix Ep-long

ro3Dn = ro3Dn_const_M .* sqrt(  EEn );
ro3Dp = ro3Dp_const_M .* sqrt( -EEp );

[Eg_M]=meshgrid(Eg,Ep);       % put the vector Gap in a matrix E-long

[EfXX_Mn,EEn]=meshgrid(EfXX,En);
[EfXX_Mp,EEp]=meshgrid(EfXX,Ep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Voltage = EfXX(end)-EfXX(1);
%Vbi = EfL - EfR;
Fbi = (EfL-EfR)/(z(end)-z(1));
Vs  = -(EfR-EfL-Voltage)/(z(end)-z(1))*z;

ntot=0;
minErr = 1e-10;     % minimum error on the potential at which the program stop
ErrVec=[];
sumVbitotVec=1;
nloop=1;

if Video_convergence==1
    figure('position',[100 100 1000 800]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Start of the Poisson s loop %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (nloop<Nloops)
         
  %nloop;

  Vbending=Vs; 
  Vbitot=V0+Vbending;
  tau = tau0*(1 + 2^((nloop - Nloops*0.8 )/10)); % the tau will increase at each loop

  %%%%%%%%%%%%%%%%%%%%%%%%%%% matrix density calcul %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vbitot_Mn=meshgrid(Vbitot,En);
  Vbitot_Mp=meshgrid(Vbitot,Ep);

  %%%%%%%%%%%%%%%%%%% calcul of the electrons density %%%%%%%%%%%%%%%%%%%%%%%%%%

  FEc = 1./(1+exp((EEn +Vbitot_Mn -EfXX_Mn)/(kB*T/e))) ; 
  ro3DEfn = ro3Dn .* FEc  ;
  NtotX = trapz(En,ro3DEfn);

  %%%%%%%%%%%%%%%%%%%%%% calcul of the holes density %%%%%%%%%%%%%%%%%%%%%%%%%%%

  FEv = 1./(1+exp(-( EEp-Eg_M +Vbitot_Mp -EfXX_Mp )/(kB*T/e))) ; 
  ro3DEfp=ro3Dp.*FEv;   
  PtotX=trapz(Ep,ro3DEfp);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  NPtotX=NtotX-PtotX-Dop;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%% Damping injection methode %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ntot = ntot + (NPtotX-ntot)/tau;    %% It add slowly the total number of electrons in order to converge

  %%%%%%%%%%%%%%%%%%%%% Electrical Field calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%

  F   = e*cumtrapz(z,ntot)./(Epsi0*Epsi);     % integal on a nonlinear grid z
  MF  = trapz(z,F)/(z(end)-z(1));  % MF is the mean(F) function on a nonlinear grid z
  F   = F - MF - Fbi  - Voltage/(z(end)-z(1)) ;

  %%%%%%%%%%%%%%%%%%%%% New potentiel calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Vs  = -cumtrapz(z,F);      % integal on a nonlinear grid z
       
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if Video_convergence==1
    cla
    hold on
    plot(z*1e9,Vbitot   ,'b-')
    plot(z*1e9,Vbitot-Eg,'b-')
    plot(z*1e9,EfXX     ,'m-')
    xlim([0 z(end)*1e9])
    ylim([min([V0 Voltage])-2 max([V0 Voltage])+2])
    title(strcat('nloop=',num2str(nloop)))
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    pause(0.01)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%% for the plotting of the graph of the convergence %%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Err = abs(  1 - sumVbitotVec(end)/sum(Vs)  );
  sumVbitotVec(k) = sum(Vs);
  ErrVec = [ErrVec Err];
  nloop  = nloop+1;

  if Err<minErr
     Err;
     break 
  end

  OUTPUT = [Vbitot' F' NtotX' PtotX' ntot'];

end
