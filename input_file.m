%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1E18 cm-3 
% fourth column is the p doping volumique of that layer in 1E18 cm-3 
% fifth column is the number of points (meshing) of that layer 
% You have to put a resonable amount of doping! Otherwise, it will diverge!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=InP;
% 
% M=[
% AlInAs        40   1       0   20
% AlInAs       100   0       0   20
% 
% AlInAs        20   0       0   20
% InGaAs30      30   0       0   20
% AlInAs        20   0       0   20
% 
% AlInAs       100   0       0   20
% AlInAs        40   1       0   20
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 1;   %% it chooses the layer number at which the Fermi level breaks on the right
% Nloops  = 300;            %% number of loops at which it stops
% tau0    = 100;            %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=InP;
% 
% M=[
% AlInAs        40   0      10   50
% AlInAs       100   0       0   100
% 
% AlInAs        20   0       0   10
% InGaAs30      30   0       0   10
% AlInAs        20   0       0   10
% InGaAs30      30   0       0   10
% AlInAs        20   0       0   10
% InGaAs30      30   0       0   10
% AlInAs        20   0       0   10
% InGaAs30      30   0       0   10
% AlInAs        20   0       0   10
% 
% AlInAs       100   0       0   100
% AlInAs        40  20       0   50
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 1;   %% it chooses the layer number at which the Fermi level breaks on the right
% Nloops  = 500;            %% number of loops at which it stops
% tau0    = 150;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=GaAs;
% 
% M=[
% GaAs        40   1       0   20
% GaAs        40   0       0   20
% 
% GaAs        20     0       0   20
% InGaAs30      30   0       0   120
% GaAs        20      0       0   20
% 
% GaAs        40   0       0   20
% GaAs        40   1       0   20
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 6;   %% it chooses the layer number at which the Fermi level breaks on the right
% Nloops  = 150;            %% number of loops at which it stops
% tau0    = 50;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=GaAs;
% 
% M=[
% GaAs        40   0       1   20
% GaAs        40   0       0   20
% 
% GaAs        20     0       0   20
% InGaAs30      30   0       0   120
% GaAs        20      0       0   20
% InGaAs30      30   0       0   120
% GaAs        20      0       0   20
% InGaAs30      30   0       0   120
% GaAs        20      0       0   20
% 
% GaAs        40   0       0   20
% GaAs        40   1       0   20
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left
% Fermi_layerbreak_R = 6;   %% it chooses the layer number at which the Fermi level breaks on the right
% Nloops  = 150;            %% number of loops at which it stops
% tau0    = 50;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NP junction on Si

% substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% M=[
% Si   40  1     0    30
% Si   30  0     1.3  50
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
% Fermi_layerbreak_R = 1;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
% Nloops  = 150;            %% number of loops at which it stops
% tau0    = 50;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PN junction on Si

% substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% M=[
% Si   50   0       1   30
% Si   100  0.2     0   50
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
% Fermi_layerbreak_R = 1;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
% Nloops  = 50;             %% number of loops at which it stops
% tau0    = 15;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIN junction on Si

% substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% M=[
% Si   50  0       1   30
% Si   50  0       0   10
% Si   50  0.8     0   50
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
% Fermi_layerbreak_R = 2;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
% Nloops  = 50;            %% number of loops at which it stops
% tau0    = 15;            %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NPN junction on Si

% substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% M=[
% Si    50   1       0     20
% Si   200   0       0.5   60
% Si   150   0.2     0     50
% ];
% 
% Fermi_layerbreak_L = 1;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
% Fermi_layerbreak_R = 2;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
% Nloops  = 500;            %% number of loops at which it stops
% tau0    = 150;            %% Damping value

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5V N-MOS with 12nm oxide gate

substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)

M = [
Poly    3     50       0      3
Poly    2     50       0     15
Oxide  12      0       0      3
Si      0.5    0       0.2   15
Si      1.0    0       0.2   10
Si      2.5    0       0.2    5
Si      100    0       0.2   30
];

Fermi_layerbreak_L = 2;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
Fermi_layerbreak_R = 3;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
Nloops  = 50;             %% number of loops at which it stops
tau0    = 15;             %% Damping value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.8V P-MOS with 3.5nm oxide gate

% substrate=Si;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% 
% M = [
% Poly   3     0       50     5
% Poly   2     0       50     5
% Oxide  3.5   0        0     3
% Si     0.5   0.5      0    15
% Si     1.0   0.5      0    10
% Si     2.5   0.5      0     5
% Si    70     0.5      0    30
% ];
% 
% Fermi_layerbreak_L = 2;   %% it chooses the layer number at which the Fermi level breaks on the left in case of applied voltage
% Fermi_layerbreak_R = 3;   %% it chooses the layer number at which the Fermi level breaks on the right in case of applied voltage
% Nloops  = 50;             %% number of loops at which it stops
% tau0    = 15;             %% Damping value
