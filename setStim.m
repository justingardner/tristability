function p = setStim(p)

switch p.cond
    case 13 %Tristable zone
        modulator = ones([1 p.nt]);%p.nt is n steps
        %onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        %modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient

        OriL = [.75 .85 0]'; %Orientation in left  eye: one on and one off
        OriR = [0 .85 .75]'; %Orientation in right eye: one off and one on
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,3,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,3,1);

    case 14 %Tristable zone with orientation tuning curves for inputs
        modulator = ones([1 p.nt]);%p.nt is n steps
        disp = 33;
        %80x should make it so they respond with activity 1 to their ideal stimulus orientation when the input strength is [1 1]
        OriL = 80.*[tuningResponse(p.nodeCenters(1),-.5*disp) tuningResponse(p.nodeCenters(2),-.5*disp) tuningResponse(p.nodeCenters(3),-.5*disp)]';
        OriR = 80.*[tuningResponse(p.nodeCenters(1),.5*disp) tuningResponse(p.nodeCenters(2),.5*disp) tuningResponse(p.nodeCenters(3),.5*disp)]';
        
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,3,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,3,1);
        
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);

    case 17 % Monocular plaid with orientation tuning simulation
        modulator = ones([1 p.nt]);%p.nt is n steps
        disp = 90;
        OriL = 80.*[tuningResponse(p.nodeCenters(1),-.5*disp) tuningResponse(p.nodeCenters(2),-.5*disp) tuningResponse(p.nodeCenters(3),.5*disp)];
        OriR = 80.*[tuningResponse(p.nodeCenters(1),-20.*disp) tuningResponse(p.nodeCenters(2),-20.*disp) tuningResponse(p.nodeCenters(3),-20.*disp)];
        
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,3,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,3,1);
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);
        
    case 11 %'Rotation' from binocular fusion to rivalry
        modulator = ones([1 p.nt]);%p.nt is n stepsp.r{lay}(:,idx-1)
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        
        OriL = [0 1 1]'; %Orientation in left  eye: one on and one off
        OriR = [1 1 0]'; %Orientation in right eye: one off and one on
        rampDown = 1-1./p.nt:-1./p.nt:0;
        rampUp = 1./p.nt:1./p.nt:1;
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,3,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,3,1);
        
        
        %Each eye will start at the central orientation, and switch to
        %opposing ones by the end.
        p.stimL(2,:) = p.stimL(2,:).*rampDown;
        p.stimL(3,:) = p.stimL(3,:).*rampUp;
        p.stimR(2,:) = p.stimR(2,:).*rampDown;
        p.stimR(1,:) = p.stimR(1,:).*rampUp;

    case 15 %'Rotation from fusion to rivalry, with tuning curves simulated
       
        maxDisp = 50;
        minDisp = 10;
        dDisp = maxDisp-minDisp;
        oriDisp = minDisp+dDisp/p.nt:dDisp/p.nt:maxDisp;
        OriL = [.5*oriDisp;.5*oriDisp;.5*oriDisp];
        OriR = [-.5*oriDisp;-.5*oriDisp;-.5*oriDisp];
        
        for i=1:p.nt
            for n=1:p.ntheta
                OriL(n,i) = tuningResponse(p.nodeCenters(n),OriL(n,i));
                OriR(n,i) = tuningResponse(p.nodeCenters(n),OriR(n,i));
            end
        end
        
        p.stimL = 80.*OriL;
        p.stimR = 80.*OriR;
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);

    case 16 %'Rotation from rivalry to fusion, with tuning curves simulated

        maxDisp = 50;
        minDisp = 10;
        dDisp = maxDisp-minDisp;
        oriDisp = maxDisp-dDisp/p.nt:-1.*dDisp/p.nt:minDisp;
        OriL = [.5*oriDisp;.5*oriDisp;.5*oriDisp];
        OriR = [-.5*oriDisp;-.5*oriDisp;-.5*oriDisp];
        
        for i=1:p.nt
            for n=1:p.ntheta
                OriL(n,i) = tuningResponse(p.nodeCenters(n),OriL(n,i));
                OriR(n,i) = tuningResponse(p.nodeCenters(n),OriR(n,i));
            end
        end
        
        p.stimL = 80.*OriL;
        p.stimR = 80.*OriR;
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);
        
    case 12 %'Rotation' from rivalry to binocular fusion
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        
        OriL = [0 1 1]'; %Orientation in left  eye: one on and one off
        OriR = [1 1 0]'; %Orientation in right eye: one off and one on
        rampDown = 1-1./p.nt:-1./p.nt:0;
        rampUp = 1./p.nt:1./p.nt:1;
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,3,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,3,1);
        
        %Each eye will start with opposite orientations, and shift to the
        %central one by the end
        p.stimL(2,:) = p.stimL(2,:).*rampUp;
        p.stimL(3,:) = p.stimL(3,:).*rampDown;
        p.stimR(2,:) = p.stimR(2,:).*rampUp;
        p.stimR(1,:) = p.stimR(1,:).*rampDown;
    case 10 %Binocular fusion
        modulator = ones([1 p.nt]);
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        
        OriL = [1 0]'; %Same orientation in both eyes
        OriR = [1 0]'; 
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,2,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,2,1);
    case 8 %'Rotation' from plaid to rivalry
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        
        OriL = [1 1]'; %Orientation in left  eye: one on and one off
        OriR = [1 1]'; %Orientation in right eye: one off and one on
        rampDown = 1-1./p.nt:-1./p.nt:0;
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,2,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,2,1);
        
        %Now decrease input to one direction in each eye over time, so that
        %we go from both orientations in both eyes fully on (binocular
        %plaid) to just one competing one in each eye (rivalry)
        p.stimL(2,:) = p.stimL(2,:).*rampDown;
        p.stimR(1,:) = p.stimR(1,:).*rampDown;
    case 9 %'Rotation' from rivalry to plaid
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        
        OriL = [1 1]'; %Orientation in left  eye: one on and one off
        OriR = [1 1]'; %Orientation in right eye: one off and one on
        rampUp = 1./p.nt:1./p.nt:1;
        p.stimL = kron(ones(1,p.nt), OriL) .* repmat(modulator,2,1);
        p.stimR = kron(ones(1,p.nt), OriR) .* repmat(modulator,2,1);
        
        %Now increase input to one direction in each eye over time, so that
        %we go from just one competing one in each eye (rivalry) to both 
        %orientations in both eyes fully on (binocular plaid)
        p.stimL(2,:) = p.stimL(2,:).*rampUp;
        p.stimR(1,:) = p.stimR(1,:).*rampUp;
    case {1,2} %Binocular Rivalry 
        OriL = [45 45 45]'; %Orientation in left  eye: one on and one off
        OriR = [-45 -45 -45]'; %Orientation in right eye: one off and one on
        
        OriL = repmat(OriL,[1,p.nt]);
        OriR = repmat(OriR,[1,p.nt]);
        
        for i=1:p.nt
            for n=1:p.ntheta
                OriL(n,i) = tuningResponse(p.nodeCenters(n),OriL(n,i));
                OriR(n,i) = tuningResponse(p.nodeCenters(n),OriR(n,i));
            end
        end
        
        p.stimL = 80.*OriL;
        p.stimR = 80.*OriR;
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);
        
    case {3,4} %Monocular Plaid
        OriL = [-45 45 45]'; %Monocular plaid in left eye
        OriR = [-100 -100 -100]'; %Nothing in right eye
        
        OriL = repmat(OriL,[1,p.nt]);
        OriR = repmat(OriR,[1,p.nt]);
        
        for i=1:p.nt
            for n=1:p.ntheta
                OriL(n,i) = tuningResponse(p.nodeCenters(n),OriL(n,i));
                OriR(n,i) = tuningResponse(p.nodeCenters(n),OriR(n,i));
            end
        end
        
        p.stimL = 80.*OriL;
        p.stimR = 80.*OriR;
        
        modulator = ones([1 p.nt]);%p.nt is n steps
        onset     = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4); %skip this line if ignore onset transient
        modulator(1:length(onset)) = modulator(1:length(onset))+onset*p.alphaAmp; %skip this line if ignore onset transient
        p.stimL = p.stimL .* repmat(modulator,3,1);
        p.stimR = p.stimR .* repmat(modulator,3,1);
        
    case {5,6} %Blank
        ts_L     = zeros(p.ntheta,p.nt); %preallocate stimulus for left eye
        ts_R     = zeros(p.ntheta,p.nt); %preallocate stimulus for right eye
        disp = 90;
        %80x should make it so they respond with activity 1 to their ideal stimulus orientation when the input strength is [1 1]
        state_1 = 80.*[tuningResponse(p.nodeCenters(1),-.5*disp) tuningResponse(p.nodeCenters(2),-.5*disp) tuningResponse(p.nodeCenters(3),-.5*disp)]';
        state_2 = 80.*[tuningResponse(p.nodeCenters(1),.5*disp) tuningResponse(p.nodeCenters(2),.5*disp) tuningResponse(p.nodeCenters(3),.5*disp)]';
        
        ts_state = mod(floor(p.tlist/p.ISP),2)+1; %time series of the state of the stimulus (1 or 2, determine which stimulus is in which eye)
        
        %draw onset transient
        onsetIdx  = abs([1 diff(ts_state)]); %time series indexing the time point where stimulus swap
        ts_alpha  = conv(onsetIdx,makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4)); %convolve onsetIdx with the shape of onset transient
        
        %add blank before swap
        modulator = ones(p.ntheta,p.nt);           %preallocate modulator of the stimulus (onset- offset- transient, and blank)
        changeIdx = find(abs([0 diff(ts_state)])); %time point where stimulus swap
        nblank    = round(p.blank/p.dt);           %length of blank in unit of dt
        for i = 1:nblank
            modulator(:,changeIdx-i)= 0; %insert blank before swap
        end
        
        %draw offset transient
        offsetIdx = double(([1 diff(modulator(1,:))])==-1);   %time series indexing the time point of stimulus offset
        ts_offset = conv(offsetIdx,makeoffset(p.dt,p.tan_t)); %convolve onsetIdx with the shape of offset decay
        
        modulator = modulator + repmat(ts_alpha(1:p.nt),p.ntheta,1)*p.alphaAmp + repmat(ts_offset(1:p.nt),p.ntheta,1);
        ts_L(:,ts_state==1) = repmat(state_1,1,sum(ts_state==1));
        ts_L(:,ts_state==2) = repmat(state_2,1,sum(ts_state==2));
        ts_R(:,ts_state==1) = repmat(state_2,1,sum(ts_state==1));
        ts_R(:,ts_state==2) = repmat(state_1,1,sum(ts_state==2));
        
        p.stimL = ts_L .* modulator * p.input(1);
        p.stimR = ts_R .* modulator * p.input(2);
        
    case 7 %Flicker and swap
        ts_L = zeros(p.ntheta,p.nt); %preallocate stimulus for left eye
        ts_R = zeros(p.ntheta,p.nt); %preallocate stimulus for right eye
        disp = 90;
        %80x should make it so they respond with activity 1 to their ideal stimulus orientation when the input strength is [1 1]
        state_1 = 80.*[tuningResponse(p.nodeCenters(1),-.5*disp) tuningResponse(p.nodeCenters(2),-.5*disp) tuningResponse(p.nodeCenters(3),-.5*disp)]';
        state_2 = 80.*[tuningResponse(p.nodeCenters(1),.5*disp) tuningResponse(p.nodeCenters(2),.5*disp) tuningResponse(p.nodeCenters(3),.5*disp)]';
        
        %draw flicker
        nframe_on  = round(1000/p.fHz/p.dt/2); %number of dt per on-cycle in the flikcer
        nrep       = ceil(p.nt / (nframe_on*2));
        flickerIdx = kron(ones(1,nrep),[ones(1,nframe_on) zeros(1, nframe_on)]); %we use even duty cycle
        flickerIdx = flickerIdx(1:p.nt);
        
        %add onset transient and offset decay to the flicker
        alpha = makealpha(p.dt,p.alpha_t*100,p.alpha_t,10e-4);
        alpha = alpha(1:min(nframe_on,length(alpha)));
        onsetIdx  = double(([1 diff(flickerIdx)])==1);
        offsetIdx = double(([1 diff(flickerIdx)])==-1);
        ts_alpha  = conv(onsetIdx,alpha);
        ts_offset = conv(offsetIdx,makeoffset(p.dt,p.tan_t));
        flickerIdx = repmat(flickerIdx+ts_alpha(1:p.nt)*p.alphaAmp+ts_offset(1:p.nt), p.ntheta,1);
        
        swapCycle = p.fHz/p.sHz;
        ts_state  = mod(floor(((1:p.nt)-1)/(swapCycle*nframe_on*2)),2)+1;
        ts_L(:,ts_state==1) = repmat(state_1,1,sum(ts_state==1));
        ts_L(:,ts_state==2) = repmat(state_2,1,sum(ts_state==2));
        ts_R(:,ts_state==1) = repmat(state_2,1,sum(ts_state==1));
        ts_R(:,ts_state==2) = repmat(state_1,1,sum(ts_state==2));
        p.stimL = ts_L.*flickerIdx * p.input(1);
        p.stimR = ts_R.*flickerIdx * p.input(2);
end

    function out = tuningResponse(mu,ori)%Returns responses of a particular gaussian centered on mu, for input of ori
        FWHM = 30.; %Set desired full-width half max
        sd = FWHM/2.355; %Set SD based on desired full width half max
        out = 1/(2*pi*sd)*exp(-(ori-mu).^2/(2*sd^2));
    end
    function alpha = makealpha(dt,T,tau,bound)
        if ~exist('bound','var')
            bound = 10e-4;
        end
        tlist = 0:dt:T;
        alpha = tlist./tau.*exp(1-tlist/tau);
        alpha(tlist > tau & alpha<bound) = [];
    end
    function offset = makeoffset(dt,duration)
        xtanh  = linspace(pi-1,-pi+0.4,(duration+2)/dt);
        offset = 1/2*tanh(xtanh)+1/2;
        offset = offset - min(offset);
        offset = offset / max(offset);
        offset = offset(2:end-1);
    end
end