function [spOut, vOut] = LSOmodelASpIF(spEx, spIn, DT, Iext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Active Leaky Integrate-and-fire Model with Spike-mimicking current of LSO
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
gLL = 14.4e-6; % 14.4*1e-6[mS] = 14.4[nS] leak conductance 
gKL = 21.6e-6; % 21.6*1e-6[mS] = 21.6[nS] KLVA conductance
EL = -56.0; % [mV] leak reversal potential 
EK = -75.0; % [mV] potassium reversal potential
Vth = -45.8; % [mV] threshold

% refractory period 
Tref = 2.0; % [ms] 
Nref = round(Tref/DT); % steps
refCount = 0; % refractory counter

% spike-associated current
Asp1 = 24.0e-3; % 24.0*1e-6[uA] = 24.0[nA] amplitude (for depolarization)
Asp2 = 12.0e-3; % 12.0*1e-6[uA] = 12.0[nA] amplitude (for repolarization)
Tsp1 = 0.15; % [ms] time constant (for depolarization)
Tsp2 = 0.30; % [ms] time constant (for repolarization)
Dsp1 = exp(-DT/Tsp1); % decay factor (for depolarization)
Dsp2 = exp(-DT/Tsp2); % decay factor (for repolarization)

% synaptic input parameters
Aex =  3.5; % [nS] excitatory synaptic input amplitude 
Ain = 12.0; % [nS] inhibitory synaptic input amplitude 
Tex = 0.16; % [ms] excitatory synaptic input time constant 
Tin = 0.32; % [ms] inhibitory synaptic input time constant 
Eex =  0.0; % [mV] excitatory synaptic input reversal potential
Ein =-75.0; % [mV] inhibitory synaptic input reversal potential  

% initialization 
Vinit = -60.6; 
vOut(1) = Vinit; 
xKL = infKL(Vinit); % KLVA activation variable
Isp1 = 0.0; % spike associated current (for depolarization)
Isp2 = 0.0; % spike associated current (for repolarization)


%% synaptic inputs 
gEx = AlphaSynapse(spEx, DT, Aex, Tex); % [nS] 
gIn = AlphaSynapse(spIn, DT, Ain, Tin); % [nS] 

%% main loop 
for t = 1:Nsteps-1

 % calculating all values at time t
 u = vOut(t); % potential 
 iLL = gLL * (EL-u); % leak current
 iKL = gKL * xKL * (EK-u); % KLVA current 
 iEx = gEx(t) * (Eex-u) * 1e-6; % [uA] excitatory synaptic current 
 iIn = gIn(t) * (Ein-u) * 1e-6; % [uA] inhibitory synaptic current
 iSp = Isp1 - Isp2; % [uA] spike-associated current 
 itot = iLL + iKL + iEx + iIn + iSp + jext(t); % total current 

 % step forward 
 v = vOut(t) + (itot/Cs) * DT; 
 xKL = xKL + ( infKL(u) - xKL )/ tauKL(u) * DT; 
 Isp1 = Isp1 * Dsp1; 
 Isp2 = Isp2 * Dsp2; 

 % check for threshold crossing 
 % (1) if in refractory period, then count no spike and decrement counter 
 if(refCount>0) 
  spOut(t) = 0; refCount = refCount-1;  
 % (2) if threshold is reached, count a spike, set the refractory counter
 %     and initiate a spike-associated current 
 elseif(v>=Vth) 
  spOut(t) = 1; refCount = Nref; 
  Isp1 = Isp1 + Asp1; 
  Isp2 = Isp2 + Asp2; 
 % (3) if no threshold crossing happened, then no spike output 
 else 
  spOut(t) = 0; 
 end 

 % save membrane potential 
 vOut(t+1) = v; 

end

%% end of main function 
end 

%% Internal Functions for KLVA Kinetics 
function x = alphaKL(v) % KLVA activation 
  x = 0.5 * exp(   ( v + 50.0 ) / 16.0 );  
end
function x = betaKL(v) % KLVA activation 
  x = 0.5 * exp( - ( v + 50.0 ) / 16.0 );   
end
function x = tauKL(v) % KLVA time constant 
  x = 1.0 ./ ( alphaKL(v) + betaKL(v) ); 
end
function x = infKL(v) % KLVA steady state value 
  a = alphaKL(v);
  x = a ./ ( a + betaKL(v) ); 
end
