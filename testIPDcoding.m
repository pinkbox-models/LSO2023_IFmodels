%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testIPDcoding.m --- for simulating IPD/ITD coding of LSO models  
% written by Go Ashida, December 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: Running this script may take several minutes. Be patient. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% time parameters
DT = 0.002; % [ms] 
Tinit = 80.0; % [ms]
Tmain = 2000.0; % [ms]
Tlast = 20.0; % [ms]
Ninit = round(Tinit/DT); % steps 
Nmain = round(Tmain/DT); % steps 
Nlast = round(Tlast/DT); % steps 
Ntot = Ninit+Nmain+Nlast; 
tv = (0:Ntot)*DT; % time vector [ms]; caution:length=Ntot+1
lmain = logical( [zeros(1,Ninit),ones(1,Nmain),zeros(1,Nlast+1)] );

%% input parameters
FQ = 300; % [Hz] modulation frequency
RT = 180-0.03*FQ; % frequency-dependent input rate 
VS = 0.65 * (1-exp( (FQ-2000)/500) ) ./ (1+exp( (FQ-2000)/500) ); % VS
P0 = 0; % zero initial phase 
Mex = 20; % number of excitatory inputs 
Min = 8;  % number of inhibitory inputs 
phases = -180:15:180; % [deg] phase difference between ipsi and contra inputs 

%% data array for output rates
RToutPLk = zeros(1,length(phases));
RToutPSp = zeros(1,length(phases));
RToutPEx = zeros(1,length(phases));
RToutALk = zeros(1,length(phases));
RToutASp = zeros(1,length(phases));
RToutAEx = zeros(1,length(phases));

%% spike input vectors
spEx  = sum( PhaseLock(Mex,length(tv),FQ,VS,RT,P0,DT), 1 );
spTmp = sum( PhaseLock(Min,length(tv),FQ,VS,RT,P0,DT), 1 );

%% main loop 
for i = 1:length(phases)
 
 % time lag between excitatory and inhibitory inputs 
 pp = phases(i); % [deg] 
 Tlag = 1000.0/FQ*pp/360; % [ms]
 Nlag = round(Tlag/DT); 

 % assign inhibitory vector according to phase difference  
 if(Nlag<=0) 
   spIn = [ zeros(1,-Nlag), spTmp(1:end+Nlag) ];    
 else  
   spIn = [ spTmp(Nlag:end), zeros(1,Nlag) ]; 
 end

 % calling LSO models 
 [spPLk,vPLk] = LSOmodelPLkIF(spEx,spIn,DT); 
 [spPSp,vPSp] = LSOmodelPSpIF(spEx,spIn,DT); 
 [spPEx,vPEx] = LSOmodelPExIF(spEx,spIn,DT); 
 [spALk,vALk] = LSOmodelALkIF(spEx,spIn,DT); 
 [spASp,vASp] = LSOmodelASpIF(spEx,spIn,DT); 
 [spAEx,vAEx] = LSOmodelAExIF(spEx,spIn,DT); 

 % getting the main part of the response
 dataPLk = spPLk(lmain); 
 dataPSp = spPSp(lmain); 
 dataPEx = spPEx(lmain); 
 dataALk = spALk(lmain); 
 dataASp = spASp(lmain); 
 dataAEx = spAEx(lmain); 
 
 % calculating output spike rates
 RToutPLk(i) = sum(dataPLk)*1000/Tmain;
 RToutPSp(i) = sum(dataPSp)*1000/Tmain;
 RToutPEx(i) = sum(dataPEx)*1000/Tmain;
 RToutALk(i) = sum(dataALk)*1000/Tmain;
 RToutASp(i) = sum(dataASp)*1000/Tmain;
 RToutAEx(i) = sum(dataAEx)*1000/Tmain;

end

%% display results
sprintf('PLk model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutPLk),min(RToutPLk),max(RToutPLk)-min(RToutPLk))
sprintf('PSp model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutPSp),min(RToutPSp),max(RToutPSp)-min(RToutPSp))
sprintf('PEx model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutPEx),min(RToutPEx),max(RToutPEx)-min(RToutPEx))
sprintf('ALk model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutALk),min(RToutALk),max(RToutALk)-min(RToutALk))
sprintf('ASp model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutASp),min(RToutASp),max(RToutASp)-min(RToutASp))
sprintf('AEx model: max=%.2f; min=%.2f; dep=%.2f',...
 max(RToutAEx),min(RToutAEx),max(RToutAEx)-min(RToutAEx))

%% plotting 
figure(232); clf;
set(gcf,'Position',[100,100,800,600]);

% color vectors
cPLk = [0.0,0.3,0.0]; cPSp = [0.3,0.6,0.1]; cPEx = [0.6,0.9,0.2];
cALk = [0.0,0.0,0.3]; cASp = [0.3,0.1,0.6]; cAEx = [0.6,0.2,0.9];

subplot(2,3,1); cla; hold on; 
plot(phases,RToutPLk,'o-','color',cPLk);
title('Passive Leaky IF Model');
xlim([-180,180]); ylim([0,140]);

subplot(2,3,2); cla; hold on; 
plot(phases,RToutPSp,'o-','color',cPSp);
title('Passive IF with Spike Current');
xlim([-180,180]); ylim([0,140]);

subplot(2,3,3); cla; hold on; 
plot(phases,RToutPEx,'o-','color',cPEx);
title('Passive Exponential IF');
xlim([-180,180]); ylim([0,140]);

subplot(2,3,4); cla; hold on; 
plot(phases,RToutALk,'o-','color',cALk);
title('Active Leaky IF Model');
xlim([-180,180]); ylim([0,140]);
xlabel('phase difference [deg]');
ylabel('rate [spikes/sec]');

subplot(2,3,5); cla; hold on; 
plot(phases,RToutASp,'o-','color',cASp);
title('Active IF with Spike Current');
xlim([-180,180]); ylim([0,140]);

subplot(2,3,6); cla; hold on; 
plot(phases,RToutAEx,'o-','color',cAEx);
title('Active Exponential IF');
xlim([-180,180]); ylim([0,140]);

