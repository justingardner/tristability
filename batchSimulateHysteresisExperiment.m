function [ ] = batchSimulateHysteresisExperiment(noises,fattcs)
%Wrapper for simulating hysteresis experiments using the modified attention
%model that has three orientations and noise. Will save resulting data file
%with the given subjID.

%% Set condition to run


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

nNoises = size(noises,2);
nFattcs = size(fattcs,2);
subjID = 900;
subjIDs = [subjID:1:subjID+nNoises*nFattcs-1];
subjIDs = reshape(subjIDs,[nNoises,nFattcs]);

hystBatchResults = cell(nNoises,nFattcs);
%Fix SubjID
parfor noiseInd=1:nNoises
    for fattcInd = 1:nFattcs

        nFastOut = 10;%10;
        nFastIn = 10;%10;
        nSlowOut = 5;%5;
        nSlowIn = 5;%5;
        allResults = [];

        for i=1:nSlowOut
            %% 
            cond = 15; %Condition to simulate - rotating out
            %Setup the stimulus sequence and parameters
            p      = setParameters(cond,noises(noiseInd),fattcs(fattcInd),160000); %set parameters

            %Override duration to be slow
            p.T      = 160000;%total duration to simulate (ms)
            p.nt     = p.T/p.dt+1; %number of time point
            p.tlist  = 0:p.dt:p.T; %time vector (ms)

            p      = setStim(p);          %draw stimuli
            p      = initTimeSeries(p);   %preallocate data matrices
            p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
            p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

            %Run the model
            %p = n_model_tuned(p);

            simRecord = p.r{1,6};%pulls out the attn layer for this subject
            [excVals,resps] = max(simRecord(:,:));
            buttons = zeros(size(resps));
            buttons(resps==1)=124;
            buttons(resps==2)=127;
            buttons(resps==3)=125;
            specs = {4,1,1,1};%slow out
            results = {specs, buttons};

            allResults = [allResults;results];
        end

        for i=1:nSlowIn
            %% 
            cond = 16; %Condition to simulate - rotating in
            %Setup the stimulus sequence and parameters
            p      = setParameters(cond,noises(noiseInd),fattcs(fattcInd),160000); %set parameters

            %Override duration to be slow
            p.T      = 160000;%total duration to simulate (ms)
            p.nt     = p.T/p.dt+1; %number of time point
            p.tlist  = 0:p.dt:p.T; %time vector (ms)

            p      = setStim(p);          %draw stimuli
            p      = initTimeSeries(p);   %preallocate data matrices
            p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
            p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

            %Run the model
            %p = n_model_tuned(p)


            simRecord = p.r{1,6};%pulls out the attn layer for this subject
            [excVals,resps] = max(simRecord(:,:));
            buttons = zeros(size(resps));
            buttons(resps==1)=124;
            buttons(resps==2)=127;
            buttons(resps==3)=125;
            specs = {4,-1,1,1};%slow in
            results = {specs, buttons};

            allResults = [allResults;results];
        end

        for i=1:nFastOut
            %% 
            cond = 15; %Condition to simulate - rotating out
            %Setup the stimulus sequence and parameters
            p      = setParameters(cond,noises(noiseInd),fattcs(fattcInd),20000); %set parameters

            %Override duration to be fast
            p.T      = 20000;      %total duration to simulate (ms)
            p.nt     = p.T/p.dt+1; %number of time point
            p.tlist  = 0:p.dt:p.T; %time vector (ms)

            p      = setStim(p);          %draw stimuli
            p      = initTimeSeries(p);   %preallocate data matrices
            p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
            p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

            %Run the model
            %p = n_model_tuned(p)

            simRecord = p.r{1,6};%pulls out the attn layer for this subject
            [excVals,resps] = max(simRecord(:,:));
            buttons = zeros(size(resps));
            buttons(resps==1)=124;
            buttons(resps==2)=127;
            buttons(resps==3)=125;
            specs = {.5,1,1,1};%fast out
            results = {specs, buttons};

            allResults = [allResults;results];
        end

        for i=1:nFastIn
            %% 
            cond = 16; %Condition to simulate - rotating in
            %Setup the stimulus sequence and parameters
            p      = setParameters(cond,noises(noiseInd),fattcs(fattcInd),20000); %set parameters

            %Override duration to be fast
            p.T      = 20000;      %total duration to simulate (ms)
            p.nt     = p.T/p.dt+1; %number of time point
            p.tlist  = 0:p.dt:p.T; %time vector (ms)

            p      = setStim(p);          %draw stimuli
            p      = initTimeSeries(p);   %preallocate data matrices
            p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
            p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

            %Run the model
            %p = n_model_tuned(p)


            simRecord = p.r{1,6};%pulls out the attn layer for this subject
            [excVals,resps] = max(simRecord(:,:));
            buttons = zeros(size(resps));
            buttons(resps==1)=124;
            buttons(resps==2)=127;
            buttons(resps==3)=125;
            specs = {.5,-1,1,1};%fast in
            results = {specs, buttons};

            allResults = [allResults;results];
        end


        results = allResults;
        currNoise = noises(noiseInd);
        currFattc = fattcs(fattcInd);
        subjID = subjIDs(noiseInd,fattcInd);
        %save(strcat(num2str(subjID),'_hysteresis.mat'),'results','currNoise','currFattc');
        %subjID=subjID+1;
        hystBatchResults{noiseInd,fattcInd} = results;

    end
end

FileName=['~/proj/rivalry/tunedThreeNodeAttentionModel/paramTests/hystParamTest',datestr(now, 'dd-mmm-yyyy-HH:MM:SS')];
save(FileName,'hystBatchResults','noises','fattcs','-v7.3');
end

