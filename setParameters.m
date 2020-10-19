function p = setParameters(cond,noise,fattc,dur)

if(nargin<2)
    noise=1.0;
    fattc=.05;
    dur=20000;
end

p.cond      = cond;
p.condnames =  {'Attended rivalry','Unattended rivalry','Attended plaid','Unattended plaid','Blank before swap','Swap with no blank','Flicker-and-Swap','Plaid to Rivalry','Rivalry to Plaid','Binocular Fusion','Fusion to Rivalry','Rivalry to Fusion','Tristable Zone','Static orientation with tuning functions','Rotating out with tuning curves','Rotating in with tuning curves','Monocular plaid'};
                               
% experiment parameters
p.input  = [1 1];%[.5 .5];  %input strength at [left eye, right eye]; for monocular pliad, only left-eye value will be used. Monocular unit tuning curves are set to work with [1 1]
p.dt     = .5;         %time-step in forward Euler method (ms)
p.T      = dur;%20000;%1800000;%20000;%15000;      %total duration to simulate (ms)
p.nt     = p.T/p.dt+1; %number of time point
p.tlist  = 0:p.dt:p.T; %time vector (ms)
p.nx     = 1;          %simulate one location
p.ntheta = 3;          %simulate three orientations
p.noise = noise;%1.0;%20;%.004;%0;.005;        %add this much noise to the adaptation updates at monocular and third layer of model

p.nodeCenters = [-30 0 30];
%p.nodeCenters = [-45 0 45];

% some properties of simulated neurons
p.nLayers       = 6; %2 monocular + 1 binocular-summation + 2 opponency + 1 attention layer
p.rectSmoothFlag= 1; %options for half-wave rectification function. 0 is non-smoothed version; 1 is the smoothed version; both work here.
                     %if you want to implement the model using ODE solver in MATLAB (or in Mathematica), or to run the model in AUTOp07, 
                     %the smoothed version is preferred.

% model parameters
p.n_m   = 1;   %exponent for monocular neurons
p.n     = 2;   %exponent for all other neurons
p.sigma = .5;  %suppression constant (for all layers, except attention layer, see p.sigma_a below)
p.m     = 2;   %gain (scaling factor) of monocular neurons
p.wh    = 2;   %weights of self-adaptation
p.wah   = .47;   %weights of self-adaptation for attention neurons (my addition)
p.wo    = .65; %weights of mutual inhibition
p.ao    = .45; %weights of attentional competition (my addition)
p.wa    = .6;  %weights of attentional modulation
p.tau_s = 5;   %time constant for monocular and binocular-summation neurons
               %The range 1-10 ms has been tested, and generated similar results
p.tau_o = 20;  %time constant for opponency neurons
p.tau_a = 150; %time constant for attention
p.tau_h = 2000;%time constant for adaptation
p.mfac = 0;%Amount of facilitation for monocular layers
p.fattc = fattc;%.05;%.22;%Amount of facilitation for attention layer
p.fac = 0;%Amount of facilitation for all other layers

% some stuff for attention layer
p.sigma_a = .2;          %suppression constant for attention layer
%p.aKernel = [1 -1;-1 1]; %weight from binocular summation neurons to attention neurons (subtraction in Eq.3 in the paper)
%p.aKernel = [1 -1 -1;-1 1 -1;-1 -1 1]; %weights from binocular summation neurons to attention neurons (subtraction in Eq.3 in the paper)
%p.aKernel = [1 -.5 -.5;-.5 1 -.5;-.5 -.5 1]; %weights from binocular summation neurons to attention neurons (subtraction in Eq.3 in the paper)


%% Set up transient/decay at stimulus onset and offset
p.alpha_t  = 3;   %parameter that influences onset transient duration (shape)
p.alphaAmp = .5;  %onset transient amplitude (relative to constant input)
p.tan_t    = 30;  %parameter that influences offset decay duration

%% condition specific settings
switch cond
    case 1 %attended binocualr rivalry
        
    case 2 %unattended binocualr rivalry
        p.wa = 0;
    case 3 %attended plaid 
        
    case 4 %un attended plaid
        p.wa = 0;
    
    %Check the sitmulus sequence and make sure that it is correct, if you change any stimulus setting for condition 5-7 here.
    case 5 %blank before swap
        p.ISP   = 333; %interswap period (ms)
        p.blank = 150; %blank duration inserted before swap (ms); has to be smaller than ISP
    case 6 %swap with no blank (static image).
        p.ISP   = 333; %interswap period (ms)
        p.blank = 0;
    case 7 %flicker and swap
        p.fHz = 18;    %flicker rate
        p.sHz = 3;     %swap rate
end

end