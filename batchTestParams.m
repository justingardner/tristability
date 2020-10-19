function batchTestParams(noises,fattcs)

%Conditions to simulate (seventeen conditions are available here)
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
conds = [1 2 3 4 5 14 15 16 15 16];

nNoises = size(noises,2);
nFattcs = size(fattcs,2);
nConds = size(conds,2);

binSumLayer = cell(nNoises,nFattcs,10);

condNames = {'Rivalrous attended','Rivalrous unattended','Plaid attended','Plaid unattended','Swaps with blanks','Tristable','Fast rot out','Fast rot in','Slow rot out','Slow rot in'};
durs = [23000 23000 23000 23000 23000 23000 20000 20000 160000 160000];

for n=1:nNoises
    noise = noises(n)
    for f=1:nFattcs
        fattc = fattcs(f)
        for cond=1:nConds

            %% 
            %Setup the stimulus sequence and parameters
            p      = setParameters(conds(cond),noise,fattc,durs(cond)); %set parameters
            p      = setStim(p);          %draw stimuli
            p      = initTimeSeries(p);   %preallocate data matrices
            p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
            p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

            %Run the model
            fprintf('%s / input strength: %1.2f %1.2f \n', condNames{cond}, p.input(1), p.input(2));
            p = n_model_tuned(p)
            binSumLayer{n,f,cond} = p.r{3};
        end
    end
end
FileName=['~/proj/rivalry/tunedThreeNodeAttentionModel/paramTests/paramTest',datestr(now, 'dd-mmm-yyyy-HH:MM:SS')];
save(FileName,'binSumLayer','condNames','noises','fattcs','conds','durs');
