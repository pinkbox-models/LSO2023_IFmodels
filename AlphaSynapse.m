function Aout = AlphaSynapse(Sp,DT,Amp,Tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alpha function for simulating synaptic input waveforms 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%   Sp  : spike vector (to show the number of inputs at each step)
%   DT  : size of time step [ms]
%   Amp : input amplitude [arbitrary unit]
%   Tau : time constant [ms]
% Output
%   Aout : simulated synaptic input waveform (same unit as Amp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by GA, 2017-2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating output vector
Aout = zeros(1,length(Sp));

% initial values 
x = 0;
y = 0; 

% factors used for calculation
aa = Amp * exp(1-DT/Tau) / Tau; % amplitude factor for normalizing the input 
ee = exp(-DT/Tau); % decrease factor used at each time step

% step-by-step calculation 
for t = 1:length(Sp)
 y = ee*y + aa*Sp(t); 
 x = ee*x + DT*y;
 Aout(t) = x; % store data
end
