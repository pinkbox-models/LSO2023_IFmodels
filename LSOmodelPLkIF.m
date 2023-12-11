function [spOut, vOut] = LSOmodelPLkIF(spEx, spIn, DT, Iext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passive Leaky Integrate-and-fire Model of LSO (with Voltage Reset at Spiking)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulation parameters
% If Iext is not given, then create it as an empty vector
if ~exist('Iext','var')
 Iext = []; 
end

% Determine the simulation time length 
if length(Iext)>1
 Nsteps = min([length(spEx), length(spIn), length(Iext)]);
else
 Nsteps = min([length(spEx), length(spIn)]);
end

% Make external input vector depending on the input argument 
if length(Iext)>1 % if the Iext is a vector, convert it from [pA] to [uA]
 jext = Iext * 1e-6; 
elseif length(Iext)==1 % if Iext is a scalar, then set a constant input 
 jext = ones(1,Nsteps) * Iext(1) * 1e-6; 
else % if Iext is an empty vector, then use zero input
 jext = zeros(1,Nsteps); 
end

% Create output vectors 
spOut = zeros(1,Nsteps);
vOut = zeros(1,Nsteps);

%% model parameters

% membrane parameters
Cs  = 24.0e-6; % 24.0*1e-6[uF] = 24.0[pF] 
gLL = 26.4e-6; % 26.4*1e-6[mS] = 26.4[nS] 
EL = -60.0; % [mV] leak reversal potential 
Vth = -45.1; % [mV] threshold 

% refractory period 
Tref = 2.0; % [ms] 
Nref = round(Tref/DT); % steps
refCount = 0; % refractory counter
Vref = -60.0; % [mV] reset (holding) potential for refractory period

% synaptic input parameters
Aex =  3.5; % [nS] excitatory synaptic input amplitude 
Ain = 12.0; % [nS] inhibitory synaptic input amplitude 
Tex = 0.16; % [ms] excitatory synaptic input time constant 
Tin = 0.32; % [ms] inhibitory synaptic input time constant 
Eex =  0.0; % [mV] excitatory synaptic input reversal potential
Ein =-75.0; % [mV] inhibitory synaptic input reversal potential  

% initialization
Vinit = EL; 
vOut(1) = Vinit; 

%% synaptic inputs 
gEx = AlphaSynapse(spEx, DT, Aex, Tex); % [nS] 
gIn = AlphaSynapse(spIn, DT, Ain, Tin); % [nS] 

%% main loop 
for t = 1:Nsteps-1

 % calculating all values at time t
 u = vOut(t); % potential 
 iLL = gLL * (EL-u); % leak current
 iEx = gEx(t) * (Eex-u) * 1e-6; % [uA] excitatory synaptic current 
 iIn = gIn(t) * (Ein-u) * 1e-6; % [uA] inhibitory synaptic current
 itot = iLL + iEx + iIn + jext(t); % total current 
  
 % step forward 
 v = vOut(t) + (itot/Cs) * DT; 

 % check for threshold crossing 
 % (1) if in refractory period, then hold voltage at reset and decrement counter
 if(refCount>0) 
  spOut(t) = 0; v = Vref; refCount = refCount-1;  
 % (2) if threshold is reached, count a spike and set the refractory counter
 elseif(v>=Vth) 
  spOut(t) = 1; refCount = Nref; 
 % (3) if no threshold crossing happened, then no spike output 
 else 
  spOut(t) = 0; 
 end 

 % save membrane potential 
 vOut(t+1) = v; 

end
