%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testSynapticInput.m --- for simulating synaptic potentials 
% written by Go Ashida, December 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time parameters
DT = 0.002; % [ms] time step
T0 = 5.0; % [ms] time before stimulus
T1 = 15.0; % [ms] time after stimulus
N0 = round(T0/DT); % steps 
N1 = round(T1/DT); % steps 
Ntot = N0 + N1; 
tv = (0:Ntot-1)*DT - T0; % time vector [ms]

%% stimulus vector
spEx = zeros(1,Ntot);
spIn = zeros(1,Ntot);

%% data vector
Nmax = 8; % max number of inputs 
vPLke = zeros(Nmax+1,Ntot); 
vPLki = zeros(Nmax+1,Ntot); 
vPSpe = zeros(Nmax+1,Ntot); 
vPSpi = zeros(Nmax+1,Ntot); 
vPExe = zeros(Nmax+1,Ntot); 
vPExi = zeros(Nmax+1,Ntot); 
vALke = zeros(Nmax+1,Ntot); 
vALki = zeros(Nmax+1,Ntot); 
vASpe = zeros(Nmax+1,Ntot); 
vASpi = zeros(Nmax+1,Ntot); 
vAExe = zeros(Nmax+1,Ntot); 
vAExi = zeros(Nmax+1,Ntot); 

%% main loop 

for i = 0:Nmax
  % excitatory input
  spEx(N0) = i; 
  spIn(N0) = 0;
  [~,vPLke(i+1,:)] = LSOmodelPLkIF(spEx,spIn,DT);
  [~,vPSpe(i+1,:)] = LSOmodelPSpIF(spEx,spIn,DT);
  [~,vPExe(i+1,:)] = LSOmodelPExIF(spEx,spIn,DT);
  [~,vALke(i+1,:)] = LSOmodelALkIF(spEx,spIn,DT);
  [~,vASpe(i+1,:)] = LSOmodelASpIF(spEx,spIn,DT);
  [~,vAExe(i+1,:)] = LSOmodelAExIF(spEx,spIn,DT);

  % inhibitory input
  spEx(N0) = 0; 
  spIn(N0) = i;
  [~,vPLki(i+1,:)] = LSOmodelPLkIF(spEx,spIn,DT);
  [~,vPSpi(i+1,:)] = LSOmodelPSpIF(spEx,spIn,DT);
  [~,vPExi(i+1,:)] = LSOmodelPExIF(spEx,spIn,DT);
  [~,vALki(i+1,:)] = LSOmodelALkIF(spEx,spIn,DT);
  [~,vASpi(i+1,:)] = LSOmodelASpIF(spEx,spIn,DT);
  [~,vAExi(i+1,:)] = LSOmodelAExIF(spEx,spIn,DT);
  
end


%% plotting
figure(231); clf;
set(gcf,'Position',[100,100,800,800]);

% color vectors
cPLk = [0.0,0.3,0.0]; cPSp = [0.3,0.6,0.1]; cPEx = [0.6,0.9,0.2];
cALk = [0.0,0.0,0.3]; cASp = [0.3,0.1,0.6]; cAEx = [0.6,0.2,0.9];

% PLk
subplot(4,3,1); cla; hold on; 
for i = 0:Nmax; plot(tv,vPLke(i+1,:),'-','color',cPLk); end
xlim([-1,9]); ylim([-80,0]); 
title('Passive Leaky IF Model');
subplot(4,3,4); cla; hold on; 
for i = 0:Nmax; plot(tv,vPLki(i+1,:),'-','color',cPLk); end
xlim([-1,9]); ylim([-75,-55]);

% PSp
subplot(4,3,2); cla; hold on; 
for i = 0:Nmax; plot(tv,vPSpe(i+1,:),'-','color',cPSp); end
xlim([-1,9]); ylim([-80,0]);
title('Passive IF with Spike Current');
subplot(4,3,5); cla; hold on; 
for i = 0:Nmax; plot(tv,vPSpi(i+1,:),'-','color',cPSp); end
xlim([-1,9]); ylim([-75,-55]);

% PEx
subplot(4,3,3); cla; hold on; 
for i = 0:Nmax; plot(tv,vPExe(i+1,:),'-','color',cPEx); end
xlim([-1,9]); ylim([-80,0]);
title('Passive Exponential IF');
subplot(4,3,6); cla; hold on; 
for i = 0:Nmax; plot(tv,vPExi(i+1,:),'-','color',cPEx); end
xlim([-1,9]); ylim([-75,-55]);

% ALk
subplot(4,3,7); cla; hold on; 
for i = 0:Nmax; plot(tv,vALke(i+1,:),'-','color',cALk); end
xlim([-1,9]); ylim([-80,0]);
title('Active Leaky IF Model');
subplot(4,3,10); cla; hold on; 
for i = 0:Nmax; plot(tv,vALki(i+1,:),'-','color',cALk); end
xlim([-1,9]); ylim([-75,-55]);
xlabel('time [ms]'); ylabel('potential [mV]');

% ASp
subplot(4,3,8); cla; hold on; 
for i = 0:Nmax; plot(tv,vASpe(i+1,:),'-','color',cASp); end
xlim([-1,9]); ylim([-80,0]);
title('Active IF with Spike Current');
subplot(4,3,11); cla; hold on; 
for i = 0:Nmax; plot(tv,vASpi(i+1,:),'-','color',cASp); end
xlim([-1,9]); ylim([-75,-55]);

% AEx
subplot(4,3,9); cla; hold on; 
for i = 0:Nmax; plot(tv,vAExe(i+1,:),'-','color',cAEx); end
xlim([-1,9]); ylim([-80,0]);
title('Active Exponential IF');
subplot(4,3,12); cla; hold on; 
for i = 0:Nmax; plot(tv,vAExi(i+1,:),'-','color',cAEx); end
xlim([-1,9]); ylim([-75,-55]);

