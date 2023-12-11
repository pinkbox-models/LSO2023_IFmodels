---------------------------------------------------------------------------------------
 Integrate-and-fire-type models of the lateral superior olive -- Matlab Implementation
---------------------------------------------------------------------------------------

%%% Versions %%% 
+ Ver. 0.9 (12.12.2023): Initial release of the code on GitHub. 

%%% Author %%% 
Go Ashida (University of Oldenburg) go.ashida@uni-oldenburg.de

%%% Contents %%% 
LSOmodelPLkIF : Passive Leaky IF 
LSOmodelPSpIF : Passive IF with Spike Current
LSOmodelPExIf : Passive Exponential IF
LSOmodelALkIF : Active Leaky IF 
LSOmodelASpIF : Active IF with Spike Current 
LSOmodelAExIF : Active Exponential IF
PhaseLock     : Code for generating phase-locked input sequences 
AlphaSynapse  : Code for calculating synaptic input modeled as an alpha function 
testSynapticInput : Sample code for plotting simulated synaptic inputs of each model 
testIPDcoding     : Sample code for plotting binaural phase-tuning curve of each model 
testILDcoding     : Sample code for plotting binaural intensity-tuning curve of each model 

%%% Notes %%%
The LSOmodelPLkIF and LSOmodelASpIF models are respectively the same as the 
LSOmodelPIF and LSOmodelAIF models in the 2017 package but with a revised 
refractory period of 2.0 ms. For more information of each model, see the reference below. 

%%% Reference %%% 
Ashida G, Wang T, Kretzberg J (2023/2024) To be submitted 
"Integrate-and-fire-type models of the lateral superior olive" 


%%%%%% Copyright 2023 Go Ashida (go.ashida@uni-oldenburg.de) %%%%%%%%%%%%%
% Permission is hereby granted under the Apache License, Version 2.0; 
% Users of this package must be in compliance with this license, a copy  
% of which may be obtained at http://www.apache.org/licenses/LICENSE-2.0
% This package is provided on an "AS IS" basis, WITHOUT WARRANTIES OR 
% CONDITIONS OF ANY KIND, either express or implied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

