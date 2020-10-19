function [] = testHysteresisConditions(noise,fattc)


cond = 15; %Condition to simulate - rotating out
%Setup the stimulus sequence and parameters
dur=20000;%total duration to simulate (ms)
p      = setParameters(cond,noise,fattc,dur); %set parameters

p.nt     = p.T/p.dt+1; %number of time point
p.tlist  = 0:p.dt:p.T; %time vector (ms)

p      = setStim(p);          %draw stimuli
p      = initTimeSeries(p);   %preallocate data matrices
p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

%Run the model
p = n_model_tuned(p)

 %Plot results
simRecord = p.r{1,6};%pulls out the attention layer for this subject
figure();
subplot(2,1,1);
plot(simRecord(1,:));
hold on;
plot(simRecord(2,:));
plot(simRecord(3,:));
title('Fast Out');

cond = 16; %Condition to simulate - rotating in
%Setup the stimulus sequence and parameters
dur=20000;%total duration to simulate (ms)
p      = setParameters(cond,noise,fattc,dur); %set parameters

p.nt     = p.T/p.dt+1; %number of time point
p.tlist  = 0:p.dt:p.T; %time vector (ms)

p      = setStim(p);          %draw stimuli
p      = initTimeSeries(p);   %preallocate data matrices
p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
p.i{2} = p.stimR;             %assign stimulus to the inputs of monocular layers

%Run the model
p = n_model_tuned(p)

 %Plot results
simRecord = p.r{1,6};%pulls out the attention layer for this subject
subplot(2,1,2);
plot(fliplr(simRecord(1,:)));
hold on;
plot(fliplr(simRecord(2,:)));
plot(fliplr(simRecord(3,:)));
title('Fast In (Flipped)');

end