function [] = testRivalryConditions(noise,fattc)

dur = 20000;

for condition = 1:4
    p      = setParameters(condition,noise,fattc,dur); %set parameters
    p      = setStim(p);          %draw stimuli
    p      = initTimeSeries(p);   %preallocate data matrices
    p.i{1} = p.stimL;             %assign stimulus to the inputs of monocular layers
    p.i{2} = p.stimR;
    p = n_model_tuned(p)
    figure();
    subplot(2,1,1);
    attn = p.r{6};
    plot(attn(1,:));
    title(strcat(p.condnames{condition},' attention layer'));
    hold on;
    plot(attn(2,:));
    plot(attn(3,:));


    subplot(2,1,2);
    binSum = p.r{3};
    plot(binSum(1,:));
    title(strcat(p.condnames{condition},' binSum layer'));
    hold on;
    plot(binSum(2,:));
    plot(binSum(3,:));
end
