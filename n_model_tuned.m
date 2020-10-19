function [p,results] = n_model_tuned(p)
%This is the model code for the Modified Attention model.
%p is set by setStim and setParameters


%rng(54); %For debugging

%Layers are defined as:
%Layer 1 = Left-eye monocular neurons
%Layer 2 = Right-eye monocular neurons
%Layer 3 = Binocular-summation neurons
%Layer 4 = Left-minus-right opponency neurons
%Layer 5 = Right-minus-left opponency neurons
%Layer 6 = Attention

%Old code for attention kernel
%Get cross-connectivity matrix based on tuning functions and
%centers of each node
% p.crossConns = getCrossConns(p);
% zeroedKernel = p.crossConns-(1./p.ntheta);
% newKernel = zeros(3,3);
% for i=1:3
%     newKernel(i,:) = (zeroedKernel(i,:)./max(zeroedKernel(i,:)));
% end
%p.aKernel = newKernel';

p.aKernel = [0 -1 -1; -1 0 -1; -1 -1 0];%Attention kernel with mutual inhibition

% Define the rectification functions
switch p.rectSmoothFlag
    case 0
        h = @(x)halfExp(x,1);
    case 1
        h = @(x)halfExp_smooth(x);
end

idx = 1;
for t = p.dt:p.dt:p.T
    
    idx = idx+1;
    
    % Print progress (msec) in command window
    if mod(p.tlist(idx),5000) <= p.dt/2
        counterdisp(p.tlist(idx));
    end
    
    %% Monocular Layers
    for lay = [1 2]
        % Input: for monocular layer the input is the stimulus input strength
        nNodes = size(p.nodeCenters,2);
        inp = p.i{lay}(:,idx);
        
        %Save input value
        p.inp{lay}(:,idx) = inp;
        
        % Updating excitatory drive: E
        p.d{lay}(:,idx) =h(inp.^p.n_m - p.o{lay}(:,idx-1)*p.wo) .* h(1 + p.r{6}(:,idx-1)*p.wa);%Original
        %p.d{lay}(:,idx) =h(inp.^p.n_m - p.r{6-lay}(:,idx-1)*p.wo) .* h(1 + p.r{6}(:,idx-1)*p.wa);%Node by node opponency
        %p.d{lay}(:,idx) =h(inp.^p.n_m - p.o{lay}(:,idx-1)*p.wo) .* h(1 + p.aKernel*p.r{6}(:,idx-1)*p.wa);%attention feedback with [inhibitory] cross connections
        %p.d{lay}(:,idx) =h(inp.^p.n_m * 1.22 - p.o{lay}(:,idx-1)*p.wo) .* h(1 + p.r{6}(:,idx-1)*p.wa);%With collinear facilitation
        %p.d{lay}(:,idx) =h(inp.^p.n_m - p.o{lay}(:,idx-1)*p.wo) + h(p.r{6}(:,idx-1)*p.wa);%Additive attention feedback
    end
    for lay = [1 2]
        
        % Defining normalization pool (six monocular neurons)
        pool = [p.d{1}(:,idx) p.d{2}(:,idx)];
        
        % Compute suppressive drive: S
        p.s{lay}(:,idx) = sum(pool(:)); % suppressive drive: sum over all the units in normalization pool
        
        % Asymptotic firing rate (normalization equation)
        p.f{lay}(:,idx) = p.m*p.d{lay}(:,idx) ./ (p.s{lay}(:,idx) + p.sigma.^p.n_m + p.h{lay}(:,idx-1).^p.n_m);
        
        % Update response: R
        p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_s)*((-1+p.mfac).*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));
        
        % Update adaptation: H
        p.h{lay}(:,idx) = p.h{lay}(:,idx-1) + (p.dt/p.tau_h)*(-p.h{lay}(:,idx-1) + p.r{lay}(:,idx-1)*p.wh + [randn;randn;randn].*p.noise);
    end
    
    %% Binocular-summation and Opponency Layers
    for lay = 3:5
        % Input
        switch lay
            case 3 %Binocular-summation neurons
                inp = p.r{1}(:,idx-1) + p.r{2}(:,idx-1);%Original
                %inp = (1+p.r{1}(:,idx-1)) .* (1+p.r{2}(:,idx-1)).^2.-1;%Attempt at making a sort of squaring nonlinearity
                %inp = p.crossConns*p.r{1}(:,idx-1) + p.crossConns*p.r{2}(:,idx-1);%With cross connections included
                
            case 4 %LE-RE opponency neurons
                inp = h(p.r{1}(:,idx-1) - p.r{2}(:,idx-1)); %Original
                %inp = h(p.r{1}(:,idx-1) - p.r{2}(:,idx-1)./(sum(p.r{1}(:,idx-1)))); %Normalized without magnitude                
                %inp = h(p.crossConns*p.r{1}(:,idx-1) - p.crossConns*p.r{2}(:,idx-1)); % Smooths before differencing
                %inp = h(p.r{1}(:,idx-1) - p.r{2}(:,idx-1)).*h(1 - .1*h(p.r{6}(:,idx-1))); %Non-normalized, with attention inhibition
                %inp = h((p.r{1}(:,idx-1) - p.r{2}(:,idx-1)).*(p.r{1}(:,idx-1) - p.r{2}(:,idx-1))./(sum(p.r{1}(:,idx-1)))); %Normalized keeping magnitude
                %inp = h(p.crossConns*(p.r{1}(:,idx-1) - p.r{2}(:,idx-1))); %Smooths after differencing
            case 5 %RE-LE opponency neurons
                inp = h(p.r{2}(:,idx-1) - p.r{1}(:,idx-1)); %Original
                %inp = h(p.r{2}(:,idx-1) - p.r{1}(:,idx-1)./(sum(p.r{2}(:,idx-1)))); %Normalized without magnitude                
                %inp = h(p.crossConns*p.r{2}(:,idx-1) - p.crossConns*p.r{1}(:,idx-1)); %Smooths before differencing
                %inp = h(p.r{2}(:,idx-1) - p.r{1}(:,idx-1)).*h(1 - .1*h(p.r{6}(:,idx-1))); %Non-normalized, with attention inhibition
                %inp = h((p.r{2}(:,idx-1) - p.r{1}(:,idx-1)).*(p.r{2}(:,idx-1) - p.r{1}(:,idx-1))./(sum(p.r{2}(:,idx-1)))); %Normalized keeping magnitude
                %inp = h(p.crossConns*(p.r{2}(:,idx-1) - p.r{1}(:,idx-1))); %Smooths after differencing
        end
        %Save input value
        p.inp{lay}(:,idx) = inp;
        
        % Updating excitatory drive: E
        p.d{lay}(:,idx) = inp.^p.n;
    end
    for lay = 3:5
        % Compute suppressive drive (S), asymptotic response (normalization equation) (f), update response (R)
        switch lay
            case 3
                p.s{lay}(:,idx) = p.d{lay}(:,idx);
                p.f{lay}(:,idx) = p.d{lay}(:,idx) ./ (p.s{lay}(:,idx) + p.sigma.^p.n + p.h{lay}(:,idx-1).^p.n);
                p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_s)*((-1+p.fac).*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));
            case {4,5}
                pool = p.d{lay}(:,idx);
                p.s{lay}(:,idx) = sum(pool(:));
                p.f{lay}(:,idx) = p.d{lay}(:,idx) ./ (p.s{lay}(:,idx) + p.sigma.^p.n);
                p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_o)*((-1+p.fac).*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));
        end
        
        % Update adaptation: H
        if lay == 3
            p.h{lay}(:,idx) = p.h{lay}(:,idx-1) + (p.dt/p.tau_h)*(-p.h{lay}(:,idx-1) + p.r{lay}(:,idx)*p.wh + [randn;randn;randn].*p.noise);
        end
        
        % Pooling the responses of opponency neurons for each eye, for eliciting inhibition
        if lay == 4
            p.o{2}(:,idx) = sum(p.r{lay}(:,idx)); % Inhibition sent to lay 2
        elseif lay == 5
            p.o{1}(:,idx) = sum(p.r{lay}(:,idx)); % Inhibition sent to lay 1
        end
    end
    
    %% Update attention
    for lay=6
        % This is written in an alternative form of Eq.3 in the paper.
        % This is a short-hand for combining on- and off- channels,
        % aSign keeps track of the sign of attention (attentional enhacement vs. attentional suppression)
        % so this will work regardless whether the exponent term is an an odd or even number

        inp     = p.r{3}(:,idx);%Original
        %inp = [1;1;1];%stand-in to test the dynamics of the attention layer alone
        %inp     = p.aKernel*p.r{3}(:,idx);%With inhibitory cross connections
        
        %Save input value
        p.inp{lay}(:,idx) = inp;
        
        %aDrive  = abs(p.aKernel*inp);
        %aSign   = sign(p.aKernel*inp);
        
        % Excitatory drive
        %p.d{lay}(:,idx) = aSign.*(aDrive.^p.n);%Original
        %p.d{lay}(:,idx) = inp.^p.n;
        p.d{lay}(:,idx) = h(inp.^p.n + p.ao.*p.aKernel*p.r{lay}(:,idx-1))+p.fattc.*p.r{lay}(:,idx-1);%Paperdraft
        %p.d{lay}(:,idx) = h(inp + p.ao.*p.aKernel*p.r{lay}(:,idx-1)).^p.n+p.fattc.*p.r{lay}(:,idx-1);%Intermediate
        %p.d{lay}(:,idx) = (inp + p.ao.*p.aKernel*p.r{lay}(:,idx-1)).^p.n;%Modified
        
        % Suppressive drive (S) + suppression constant (sigma)
        %p.s{lay}(:,idx) = repmat((sum(aDrive.^p.n) + p.sigma_a^p.n),p.ntheta,1);%Original
        p.s{lay}(:,idx) = repmat((sum(inp.^p.n) + p.sigma_a^p.n),p.ntheta,1);%Paperdraft
        %p.s{lay}(:,idx) = repmat((sum(h(inp))),p.ntheta,1);%Modified
        
        % Asymptotic firing rate (normalization equation)
        p.f{lay}(:,idx) = p.d{lay}(:,idx) ./(p.s{lay}(:,idx) + p.h{lay}(:,idx-1));%Original
        
        % Update responses: R
        p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_a)*(-1.*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));%Paperdraft
        %p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_a)*((-1+p.fac).*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));%Intermediate
        %p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_a).*((-1+p.fattc).*p.r{lay}(:,idx-1) + p.f{lay}(:,idx));%Original
        %p.r{lay}(:,idx) = p.r{lay}(:,idx-1) + (p.dt/p.tau_a)*(-1.*p.r{lay}(:,idx-1) + p.f{lay}(:,idx) - h(p.h{lay}(:,idx-1)));%with adaptation
        
        %Update adaptation (my addition)
        p.h{lay}(:,idx) = p.h{lay}(:,idx-1) + (p.dt/p.tau_h)*(-1.*p.h{lay}(:,idx-1) + p.r{lay}(:,idx)*p.wah + [randn;randn;randn].*p.noise);
    end
end

    function crossConns = getCrossConns(p)%returns a matrix of cross-connections between nodes based on their centers
        nNodes = size(p.nodeCenters,2);
        crossConns = zeros(nNodes,nNodes);
        for node=1:nNodes
            for otherNode = 1:nNodes
                crossConns(node,otherNode)=tuningResponse(p.nodeCenters(node),p.nodeCenters(otherNode));
            end
        end
        for row=1:nNodes
            crossConns(row,:) = crossConns(row,:).*(1/sum(crossConns(row,:)));
        end
    end

    function out = tuningResponse(mu,ori)%Returns responses of a particular gaussian centered on mu, for input of ori
        FWHM = 30.; %Set desired full-width half max
        sd = FWHM/2.355; %Set SD based on desired full width half max
        out = 1/(2*pi*sd)*exp(-(ori-mu).^2/(2*sd^2));
    end
    
    % Function for displaying progress
    function counterdisp(i)
        fprintf('%d msec \r', i);
    end

    % Rectification (non-smoothed version)
    function [x] = halfExp(base,n)
        if nargin == 1
            n=1;
        end
        x = (max(0,base)).^n;
    end

    % Rectification (smoothed version)
    function y=halfExp_smooth(x)
        thresh=0.05;
        slope=30;
        
        idx_0 = x<0;
        idx_1 = x>0;
        
        y = zeros(size(x));
        y(idx_0) = 0;
        y(idx_1) = x(idx_1).*1./(1+exp(-slope*(x(idx_1)-thresh)));
    end
end