%%This is the wrapper for running simulations of the Modified Attention
%%Model.

%% Set condition to run

cond = 2; %Condition to simulate (seventeen conditions are available here)
%1. Dichoptic gratings / attended   (attended rivalry)
%2. Dichoptic gratings / unattended (unattended rivalry)
%3. Monocular plaids / attended
%4. Monocular plaids / unattended
%5. Swapping with a blank gap before swap (blank-before-swap condition) 
%6. Swapping with no blank (static image condition) 
%7. Swapping with flicker  (flicker-and-swap condition)
%8. Binocular plaid to rivalry
%9. Rivalry to binocular plaid
%10. Binocular fusion (same orientation in both eyes)
%11. Fusion to rivalry
%12. Rivalry to fusion
%13. Tristable zone - 'near fusion'
%14. Static disparity simulation using orientation tuning simulation
%15. Rotating out simulation using orientation tuning simulation
%16. Rotating in simulation using orientation tuning simulation
%17. Monocular plaid with orientation tuning simulation


%% 
%Setup the stimulus sequence and parameters
p      = setParameters(cond); %set parameters
p      = setStim(p);          %draw stimuli
p      = initTimeSeries(p);   %preallocate data matrices
p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

%Run the model
fprintf('%s / input strength: %1.2f %1.2f \n', p.condnames{p.cond}, p.input(1), p.input(2));
p = n_model_tuned(p)


%Plot results
plotTimeSeries(p);
