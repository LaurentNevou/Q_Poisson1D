%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update  29December2020, lne %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code solves the Poison equation in 1D on inhomogeneous grid in semiconductors
% heterostructures. As a results, it gives the band bending profile for any 
% heterostructures. It is specialy usefull for CMOS, diode and whatever npn & pnp junction.
% A strain model is included. It basically shifts the conduction and valence band edge
% The code does NOT take into account the non-parabolicity of the bands. The electron
% and the hole masses are loaded from the "material.csv" file for each material.
% Additionnal material can be added in the "materialDB_ZB.csv" file

% -> II-VI and cubic nitride material parameters are available but should
% be grabt in the "Library.m" file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the code doesn t converge:
% -> increase the damping, tau0
% -> increase the amount of loops, Nloops. Nloops should be more than 3 times higher than tau0
% -> increase the amount of points
% -> increase the resolution dE
% -> increase the temperature (T=0K is very bad while T=10K is already much better)
% -> decrease the doping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h     = 6.62606896E-34;           %% Planck constant J.s
hbar  = h/(2*pi);
e     = 1.602176487E-19;          %% charge de l electron Coulomb
m0    = 9.10938188E-31;           %% electron mass kg
c     = 2.99792458e8;             %% speed of light (m/s)
Epsi0 = 8.854187817620E-12;       %% constant dielectric du vide F/m
mu0   = 1/(Epsi0*c^2);            %% permeabiliy du vide
kB    = 1.3806488E-23;            %% Boltzmann's constant (J/K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StrainModel = 1;    % Activate Strain model
T=300;              % Temperature in Kelvin (never put zero)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Turm on Graph and Saving %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for turn off
% 1 for turn on
CBand1D=1;    Xfig=10;Yfig=100;Wfig=1000;Hfig=800;
CBand3D=1;

Convergence=1;
Video_convergence=0;
TotCharge_density=0;
Charge_density=0;
band3D=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Voltage sweep %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To make it easier, dV must be constant (homogeneous voltage grid)
dV=0.02;
Voltage=0;
%Voltage=-5:dV:5;
%Voltage=-2:dV:2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy
QuaternaryAlloy;          % compute the quaternary alloy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Energy grid definition: the grid is moving respect to the bending %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Electron Energy grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

En1 = linspace( 0 , 0.35, 31 );
En2 = linspace( En1(end)+0.01 , 1, 10 );
En  = [En1 En2]; En = sort(En);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holes Energy grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ep1 = linspace( 0 , -0.25, 20 );
Ep2 = linspace( Ep1(end)-0.01 , -1, 10 );
Ep  = [Ep1 Ep2]; Ep = sort(Ep);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt   = M(:,end-3)*1e-9;    % conversion of the length from Angstrom to meter

%%%%%%%%%%%%%%%%%%%%%%% Eg = Eg0 - (a*T.^2)./(T + b) %%%%%%%%%%%%%%%%%%%%%%%%%%%
EgG  = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   % Bandgap at Gamma point
EgX  = M(:,idx_EgX)  - (M(:,idx_alphaX)*T^2) ./ (T+M(:,idx_betaX));   % Bandgap at X point
EgL  = M(:,idx_EgL)  - (M(:,idx_alphaL)*T^2) ./ (T+M(:,idx_betaL));   % Bandgap at L point

Egt=min([EgG EgX EgL],[],2);

VBOt = M(:,idx_VBO);
CBOt = Egt + VBOt;         % CBO from band gap difference and temperature
Epsit= M(:,idx_Epsi);      % used for Poisson solver only

Doptn=M(:,end-2)*1e18*1e6;  % n doping conversion from cm-3 to m-3
Doptp=M(:,end-1)*1e18*1e6;  % p doping conversion from cm-3 to m-3
Dopt=Doptn-Doptp;

Masstn = M(:,idx_me);
Masstp = M(:,idx_mhh) ;
%Masstp = ( M(:,idx_mhh).^(3/2) + M(:,idx_mlh).^(3/2) ).^(2/3) ;

Pt=M(:,end);

Ntott=Dopt.*zt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at   = M(:,idx_a);           % lattice parameter
act  = M(:,idx_ac);          % Conduction band strain offset parameter
avt  = M(:,idx_av);          % Valence band strain offset parameter
bvt  = M(:,idx_bv);          % Valence band strain offset parameter
c11t = M(:,idx_c11);         % strain parameter
c12t = M(:,idx_c12);         % strain parameter

a0   = substrate(idx_a);

if StrainModel == 1
  exxt =  (a0-at)/a0; eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
else
  exxt =  (a0-at)/a0 * 0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z and a lot of other parameters

z(1)=0; Dop(1)=Dopt(1); Mass_n(1)=Masstn(1); Mass_p(1)= Masstp(1);Eg(1)=Egt(1);Epsi(1)=Epsit(1);
ac=act(1); av=avt(1); bv=bvt(1); exx=exxt(1); ezz=ezzt(1);
dzz=1E-12;

%V0(1)=CBOt(1);
CBO(1)=CBOt(1);
VBO(1)=VBOt(1);

for i=1:length(zt)
    t=zt(i);
    zv     = linspace( z(end)+dzz , z(end) + t , Pt(i) );
    z      = [ z                        zv        ];
    %V0     = [ V0      ones(size(zv)) * CBOt(i)   ];
    CBO    = [ CBO     ones(size(zv)) * CBOt(i)   ];
    VBO    = [ VBO     ones(size(zv)) * VBOt(i)   ];
    Dop    = [ Dop     ones(size(zv)) * Dopt(i)   ];
    Mass_n = [ Mass_n  ones(size(zv)) * Masstn(i) ];
    Mass_p = [ Mass_p  ones(size(zv)) * Masstp(i) ];
    Eg     = [ Eg      ones(size(zv)) * Egt(i)    ];
    Epsi   = [ Epsi    ones(size(zv)) * Epsit(i)  ];
    ac     = [ ac      ones(size(zv)) * act(i)    ];
    av     = [ av      ones(size(zv)) * avt(i)    ];
    bv     = [ bv      ones(size(zv)) * bvt(i)    ];
    exx    = [ exx     ones(size(zv)) * exxt(i)   ];
    ezz    = [ ezz     ones(size(zv)) * ezzt(i)   ];
end

[ZZn,EEn]=meshgrid(z,En);
[ZZp,EEp]=meshgrid(z,Ep);

eyy = exx;
DCBO   = -abs(ac).*(exx+eyy+ezz) ;                      % shift of the CB due to strain
DVBOHH = +abs(av).*(exx+eyy+ezz) - abs(bv).*(exx-ezz) ; % shift of the VB-HH due to strain
%DVBOLH = +abs(av).*(exx+eyy+ezz) + abs(bv).*(exx-ezz) ; % shift of the VB-LH due to strain
%DVBOSO = +abs(av).*(exx+eyy+ezz) ;                      % shift of the VB-SO due to strain

V0 = Eg + VBO + DCBO;
shift=V0(1);
V0=V0-shift;

Eg = Eg + DCBO - DVBOHH ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Finding the boundary conditions %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EfL = Fermi3D_np( Eg(1),T,Masstn(1),Masstp(1),Dop(1) );          % Fermi level on the left side
EfR = Fermi3D_np( Eg(end),T,Masstn(end),Masstp(end),Dop(end) );  % Fermi level on the right side

EfL = EfL+V0(1);
EfR = EfR+V0(end);

distance_L=0;
for i=1:Fermi_layerbreak_L
    distance_L=distance_L + zt(i);
end

distance_R=0;
for i=1:Fermi_layerbreak_R
    distance_R=distance_R + zt(i);
end

idxg = find( abs( z-distance_L ) <1e-10 , 1); % Here the index in the vector z where the Fermi level will be broken
idxd = find( abs( z-distance_R ) <1e-10 , 1); % Here the index in the vector z where the Fermi level will be broken

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% loop over the Voltage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ErrVec=[];

for k=1:length(Voltage)     % big loop over the applied voltage
  
  display(sprintf('Simulation at : %.2f Volt' ,Voltage(k)))
  
  %%%%%%%%%%%%%%%%%%% building of the quasi Fermi level %%%%%%%%%%%%%%%%%%%%%%%%
    
  for i=1:idxg
      EfXX(i)=EfL;
  end
  for i=idxg+1:idxd-1
      EfXX(i)= Voltage(k)*z(i)/(z(idxd)-z(idxg)) + EfL - Voltage(k)*z(idxg)/(z(idxd)-z(idxg));
  end
  for i=idxd:length(z)
      EfXX(i)=EfL+Voltage(k);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Structure=[z' V0' Eg' Dop' Epsi' EfXX' Mass_n' Mass_p'];

  [OUTPUT,ro3DEfn,ro3DEfp,Err] = Poisson_f(Structure,En,Ep,T,EfL,EfR,Nloops,tau0,k,Video_convergence);
  
  Vbitot = OUTPUT(:,1)';
  F      = OUTPUT(:,2)';
  NtotX  = OUTPUT(:,3)';
  PtotX  = OUTPUT(:,4)';
  ntot   = OUTPUT(:,5)';
  ErrVec = [ErrVec Err];
  
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

if CBand1D==1 || CBand3D==1;
    figure('position',[Xfig Yfig Wfig Hfig],'color','w')

    xscale =[z(1) z(end)]*1e9;
    yscale=[min(Vbitot-Eg) max(Vbitot)+0.1]*1.5;
    
    subplot(1,1,1,'fontsize',20)
    hold on; grid on;box on;
    colormap(jet)
    col=colormap;
    
    if CBand3D==1
        
        [Vbitot_Mn]=meshgrid(Vbitot,En);
        [Vbitot_Mp]=meshgrid(Vbitot,Ep);
        [Eg_M]=meshgrid(Eg,Ep);
        grid off
        
        pcolor(ZZn*1e9,EEn+Vbitot_Mn,ro3DEfn*1e-6)
        pcolor(ZZp*1e9,EEp+Vbitot_Mp-Eg_M,ro3DEfp*1e-6)
        plot(z*1e9,Vbitot,'w-','linewidth',2)
        plot(z*1e9,Vbitot-Eg,'w-','linewidth',2)
        
        set(gca,'color',col(1,:))
        shading flat
        
        hcb=colorbar;
        title(hcb,'\fontsize{8}cm-3')

    else
        plot(z*1e9,Vbitot   ,'b-' ,'linewidth',2)
        plot(z*1e9,Vbitot-Eg,'b-' ,'linewidth',2)
                
        if StrainModel==1
            plot(z*1e9, V0-DCBO,'k--','linewidth',1)
            plot(z*1e9,V0-DCBO-(Eg-DCBO+DVBOHH),'k--','linewidth',1)
        end
    end
    
    plot(z*1e9,V0,'b--','linewidth',1)
    plot(z*1e9,V0-Eg,'b--','linewidth',1)
    plot(z*1e9,EfXX,'g-','linewidth',1)
    text(z(1)*1e9,EfXX(1)-0.1,'\color{green}Fermi')
    text(z(end)*1e9*0.95,EfXX(end)+0.1,'\color{green}Fermi')

    xlim(xscale)
    ylim(yscale)
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    if StrainModel == 1
      title(strcat('\fontsize{15}T=',num2str(T),'K; Voltage=',num2str(Voltage(end)),'V ; with STRAIN'))
    else
      title(strcat('\fontsize{15}T=',num2str(T),'K; Voltage=',num2str(Voltage(end)),'V ; without STRAIN'))
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Convergence==1;
    
    figure
    semilogy(ErrVec,'b.-')  % convergence graph
    grid on;
    xlim([0 length(ErrVec)])
    xlabel('Convergence cycle');
    ylabel('Oscillation of the Build-in potential (%)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TotCharge_density==1;
    
    figure('position',[410 50 400 400]);

    subplot(3,1,1)
    plotyy(z*1e9,ntot*1e-6,z*1e9,Vbitot)
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    legend('ntot','Vbitot')
    xlim([z(1) z(end)]*1e9)
    
    % From the next graph, one can see if the convergence parameters are well 
    % set. In order to make the convergence, the charges are added 
    % progressively and not all are injected from the begining. As a 
    % results, it can be that the simulation finishes but not all the 
    % charges were injected.
    % NtotX is the full charge density while ntot is the one used during
    % the simulation. If both do not match, either nloop have to be
    % increased or tau0 has to be decreased.
    
    subplot(3,1,2)
    plot(z*1e9,(NtotX-PtotX-Dop)*1e-6,'r')
    hold on
    plot(z*1e9,ntot*1e-6,'b')
    plot(z*1e9,ntot*1e-6 + Dop*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    legend('NPtotX','ntot','ntot+Dop')
    xlim([z(1) z(end)]*1e9)
    title('Simulations good when both curves are on each other')
    
    subplot(3,1,3)
    hold on
    plot(z*1e9,F,'g')
    xlabel('z (nm)');
    ylabel('Electric field (V/m)');
    legend('F')
    xlim([z(1) z(end)]*1e9)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Charge_density==1;

    figure('position',[820 50 400 400]);

    subplot(2,1,1)
    yscale=[1e0 1e20];
    ylim(yscale)
    %hold on
    semilogy(z*1e9,NtotX*1e-6,'r')
    hold on;grid on;box on;
    semilogy(z*1e9,PtotX*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    xlim([z(1) z(end)]*1e9)
    ylim(yscale)
    
    subplot(2,1,2)
    plot(z*1e9,NtotX*1e-6,'r')
    hold on;grid on;box on;
    plot(z*1e9,PtotX*1e-6,'g')
    xlabel('z (nm)');
    ylabel('Charge density (cm-3)');
    xlim([z(1) z(end)]*1e9)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band3D==1

    [Vbitot_Mn]=meshgrid(Vbitot,En);
    [Vbitot_Mp]=meshgrid(Vbitot,Ep);
    [Eg_M]=meshgrid(Eg,Ep);
    
    figure('position',[1230 50 600 500]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,1,1)
    hold on
    pcolor(ZZn*1e9,EEn+Vbitot_Mn,(ro3DEfn))
    plot(z*1e9,EfXX,'g.-','linewidth',1)
    plot(z*1e9,Vbitot,'r.-')
    shading flat
    colormap(jet)
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    %zlabel('density of states * Fermi distribution')
    xlim([z(1) z(end)]*1e9)
    
    subplot(2,1,2)
    hold on
    pcolor(ZZp*1e9,EEp+Vbitot_Mp-Eg_M,(ro3DEfp))
    plot(z*1e9,EfXX,'g.-','linewidth',1)
    plot(z*1e9,Vbitot-Eg,'r.-')
    shading flat
    colormap(jet)
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    %zlabel('density of states * Fermi distribution')
    xlim([z(1) z(end)]*1e9)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%