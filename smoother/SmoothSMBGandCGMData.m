function output=SmoothSMBGandCGMData(t_in,y_in_fp,y_in_cgm,varargin)

    %Parse the variable arguments
    parsedArgs = parseInputVarArgs(varargin);
    
    %Check units
    unit1 = autoDetectGlucoseUnit(y_in_fp);
    unit2 = autoDetectGlucoseUnit(y_in_cgm);
    if strcmp(unit1,unit2) == 0
       error('Supplied units are not consistent for fingerpricks and CGM')
    end

    if strcmp(unit1,'mmol_L')==0
        error('TODO unsupported unit')
    end
    
    %Handle time
    if isdatetime(t_in)
        %convert to relative time
        startDateTime = t_in(1);
        t_in=convertToRelativeTime(t_in, startDateTime);
    else
        %Assume relative time is passed in
        startDateTime = NaN;
    end

    if parsedArgs.presmooth
        %Run smoothing on each set individually
        output.smoothed_fp = SmoothSMBGData(t_in,y_in_fp,'outlierRemoval',1,'dynamicModel',2);
        output.smoothed_cgm = SmoothSMBGData(t_in,y_in_cgm,'outlierRemoval',1,'dynamicModel',2);
        t_i = output.smoothed_fp.t_i;
        y_fp = output.smoothed_fp.y_smoothed;
        var_fp = output.smoothed_fp.y_smoothed_sd.^2;
        y_cgm = output.smoothed_cgm.y_smoothed;
        var_cgm = output.smoothed_cgm.y_smoothed_sd.^2;
        %Set dynamic model to use for CGM SMBG fusion
        dynModel = augmentDynamicModelBiasAndLag(output.smoothed_fp.dynModel);
        Nstates = size(dynModel.Phi,1);
    else
        error('TODO implement no presmooth')
    end
    
    useBiasStartGuess = true;
    
    %%% Storage
    x_hat_f = nan(Nstates,length(t_i));         % A priori state vector storage, forward pass
    x_bar_f = nan(Nstates,length(t_i));         % A posteriori state vector storage, forward pass
    P_hat_f = nan(Nstates,Nstates,length(t_i));       % A priori covariance matrix storage, forward pass
    P_bar_f = nan(Nstates,Nstates,length(t_i));       % A posteriori covariance matrix storage, forward pass
    x_smoothed = nan(Nstates,length(t_i));            % State vector storage, backward pass
    P_smoothed = nan(Nstates,Nstates,length(t_i));    % Covariance matrix storage, backward pass
    Phis = nan(Nstates,Nstates,length(t_i));    
    PhiEKFs = nan(Nstates,Nstates,length(t_i));    
    
    init = false;
    endk = nan;
    for k = 1:length(t_i)
        if ~init
            if ~isnan(y_fp(k)) && ~isnan(y_cgm(k)) %We can initialise
                %%% Initialization
                xBar = zeros(Nstates,1);
                xBar(1)=y_fp(k);
                xBar(dynModel.Nin+1)=y_fp(k); %NB - the state is not set equal to the first CGM measurement, due to bias
                if(useBiasStartGuess)
                    %xBar(dynModel.Nin+2)=y_cgm(k)-y_fp(k);%Start guess on the bias - first sample
                    %xBar(dynModel.Nin+2)=mean(y_cgm-y_fp,'omitnan');%Start guess on the bias - mean difference
                    xBar(dynModel.Nin+2)=(y_cgm(k)-y_fp(k)+ mean(y_cgm-y_fp,'omitnan') )/2;%Start guess on the bias - mean between 
                   
                    %disp(sprintf('Used guess for bias: %f', xBar(dynModel.Nin+2)));%Start guess on the bias
                    
                    %The guess needs to be computed in a different way if
                    %fingerpricks are sparse (only look at the CGM
                    %measurements close to the FPs
                end
                
                xBar(dynModel.Nin+3)=dynModel.k_default;

                xHat=xBar;
                PBar=dynModel.initCov;
                PHat=PBar;
                init = true;
                startk = k;
            end
        elseif init && (isnan(y_fp(k)) || isnan(y_cgm(k))) %time to end
            endk = k-1;
            break
        else %We have initialised and it is not time to end yet, do filtering
            %%% Kalman filter forward pass
            %TU - Time update
            if dynModel.nonLinear
                [dynModel.Phi dynModel.PhiEKF ]= computePhiEKF(dynModel,xHat);
                Phis(:,:,k)=dynModel.Phi;
                PhiEKFs(:,:,k)=dynModel.PhiEKF;
                
            end
            xBar = dynModel.Phi*xHat;
            PBar = dynModel.PhiEKF*PHat*dynModel.PhiEKF' + dynModel.Q;
            %Store
            x_bar_f(:,k)=xBar;
            P_bar_f(:,:,k)=PBar;

            %MU - Measurement Update
            %Switch between measurements, both are not always available
            if ~isnan(y_fp(k)) & ~isnan(y_cgm(k)) 
                H = dynModel.H_both;
                y = [y_fp(k);y_cgm(k)];
                R = [var_fp(k) 0 ;0 var_cgm(k)];

            elseif ~isnan(y_cgm(k))
                H = dynModel.H_cgm;
                y = y_cgm(k);
                R = var_cgm(k);

            elseif ~isnan(y_fp(k))
                H = dynModel.H_fp;
                y = y_fp(k);
                R = var_fp(k);
            else
                error('Encountered empty measurement, should have been filtered out, something is wrong')
            end
            dz = y-H*xBar;
            Pz = (H*PBar*H'+R);

            %Measurement update
            K=PBar*H'/Pz;
            xHat = xBar + K*dz;
            PHat = (eye(size(PBar))-K*H)*PBar;
            [~,p] = chol(PHat);
            if p~=0
                disp(sprintf('Warning - filtering aborted at t = %f due to non-posdef covmatrix',t_i(k)))
                break
            end
            if xHat(end)>10
                %disp(sprintf('Warning - k estimated to %f in time step %d',xHat(end),k))
                xHat(end)=10;
            end
            if xHat(end)<0.01
                %disp(sprintf('Warning - k estimated to %f in time step %d',xHat(end),k))
                xHat(end)=0.01;
            end
            % Store
            x_hat_f(:,k)=xHat;
            P_hat_f(:,:,k)=PHat;
        end
    end % for k
    if(isnan(endk))
        endk=length(t_i);
    end

    %%% Rauch-Tung-Striebel backward pass
    x_smoothed(:,endk)=xHat;
    P_smoothed(:,:,endk)=PHat;
    for k = endk-1:-1:startk
        if dynModel.nonLinear
            %[dynModel.Phi dynModel.PhiEKF ] = computePhiEKF(dynModel,x_smoothed(:,k+1));
            %PhiForRTS = dynModel.Phi;
            %PhiForRTS = dynModel.PhiEKF;
            %PhiForRTS = Phis(:,:,k);
            PhiForRTS = PhiEKFs(:,:,k);
        end
        C=(P_hat_f(:,:,k)*PhiForRTS')/P_bar_f(:,:,k+1);
        x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
        if x_smoothed(end,k)>10
            %disp(sprintf('Warning - k smoothed to %f in time step %d',x_smoothed(end,k),k))
            x_smoothed(end,k)=10;
        end
        if x_smoothed(end,k)<0.01
            %disp(sprintf('Warning - k smoothed to %f in time step %d',x_smoothed(end,k),k))
            x_smoothed(end,k)=0.01;
        end
        P_hat_s=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';  
        [~,p] = chol(P_hat_s);
        if p~=0
            disp(sprintf('Warning - smoothing aborted at t = %f due to non-posdef covmatrix',t_i(k)))
            break
        end
        P_smoothed(:,:,k) = P_hat_s;
    end

    %Generate output structs
    output.y_fp_smoothed = nan(size(t_i));
    output.y_fp_smoothed_sd = nan(size(t_i));
    output.y_cgm_smoothed = nan(size(t_i));
    output.y_cgm_smoothed_sd = nan(size(t_i));
    output.bias_smoothed = nan(size(t_i));
    output.bias_smoothed_sd = nan(size(t_i));
    output.lagk_smoothed = nan(size(t_i));
    output.lagk_smoothed_sd = nan(size(t_i));
    output.bias_filtered = nan(size(t_i));
    output.bias_filtered_sd = nan(size(t_i));
    output.lagk_filtered = nan(size(t_i));
    output.lagk_filtered_sd = nan(size(t_i));
    
    output.y_fp_smoothed(startk:endk) = x_smoothed(1,startk:endk);
    output.y_cgm_smoothed(startk:endk) = x_smoothed(dynModel.Nin+1,startk:endk);
    
    output.bias_smoothed(startk:endk) = x_smoothed(dynModel.Nin+2,startk:endk);
    output.lagk_smoothed(startk:endk) = x_smoothed(dynModel.Nin+3,startk:endk);
    
    output.bias_filtered(startk:endk) = x_hat_f(dynModel.Nin+2,startk:endk);
    output.lagk_filtered(startk:endk) = x_hat_f(dynModel.Nin+3,startk:endk);
    
    for k = startk:endk
        output.y_fp_smoothed_sd(k) = sqrt(P_smoothed(1,1,k));
        output.y_cgm_smoothed_sd(k) = sqrt(P_smoothed(dynModel.Nin+1,dynModel.Nin+1,k));
        
        output.bias_smoothed_sd(k) = sqrt(P_smoothed(dynModel.Nin+2,dynModel.Nin+2,k));
        output.bias_filtered_sd(k) = sqrt(P_hat_f(dynModel.Nin+2,dynModel.Nin+2,k));
        
        output.lagk_smoothed_sd(k) = sqrt(P_smoothed(dynModel.Nin+3,dynModel.Nin+3,k));
        output.lagk_filtered_sd(k) = sqrt(P_hat_f(dynModel.Nin+3,dynModel.Nin+3,k));
    end 

    output.t_i = t_i;
    if isdatetime(startDateTime)
        output.t_i_relative = t_i;
        output.t_i = convertToAbsoluteTime(t_i, startDateTime);
    end

    %Add internal states
    output.x_filtered = x_hat_f;
    output.x_smoothed = x_smoothed;

    %Add user supplied wanted time
    output.y_cgm_corrected_at_tout = closestValues(parsedArgs.tout,t_i,output.y_cgm_smoothed, startDateTime);
    output.y_cgm_corrected_sd_at_tout = closestValues(parsedArgs.tout,t_i,output.y_cgm_smoothed_sd, startDateTime);

     %debug plotting
    if parsedArgs.debugPlot
        subplot(3,1,1)
        plot(t_in,y_in_fp,'r.','DisplayName','Fingerpricks')
        hold on;
        plot(t_in,y_in_cgm,'b.','DisplayName','CGM measurements')
        plot(t_i,y_fp,'r-','DisplayName','Smoothed fingerpricks')
        plot(t_i,y_cgm,'b-','DisplayName','Smoothed CGM')
        %plot(t_i,output.y_fp_smoothed,'g-','DisplayName','Merge-smoothed FPs')
        %plot(t_i,output.y_cgm_smoothed,'k-','DisplayName','Merge-smoothed CGM')
        
        hold off;
        xlabel('Time [min]','FontWeight','bold','FontSize',12);
        ylabel('Glucose [mmol/L]','FontWeight','bold','FontSize',12);
        title('Smoothing CGM and SMBG')
        legend('show')
        subplot(3,1,2)
        plot(t_i,output.bias_smoothed,'r-','DisplayName','Bias - smoothed')
        hold on;
        plot(t_i,output.bias_smoothed-2.5*output.bias_smoothed_sd,'r--','DisplayName','Bias variance')
        plot(t_i,output.bias_smoothed+2.5*output.bias_smoothed_sd,'r--')
        plot(t_i,output.bias_filtered,'g-','DisplayName','Bias - filtered')
        plot(t_i,output.bias_filtered-2.5*output.bias_filtered_sd,'g--','DisplayName','Bias variance')
        plot(t_i,output.bias_filtered+2.5*output.bias_filtered_sd,'g--')
        ylabel('Bias','FontWeight','bold','FontSize',12);
        ylim([-3 3])
        %legend('show')
        
        subplot(3,1,3)
        plot(t_i,output.lagk_smoothed,'r-','DisplayName','Lag k - smoothed')
        hold on;
        plot(t_i,output.lagk_smoothed-2.5*output.lagk_smoothed_sd,'r--','DisplayName','Lag k variance')
        plot(t_i,output.lagk_smoothed+2.5*output.lagk_smoothed_sd,'r--')
        plot(t_i,output.lagk_filtered,'g-','DisplayName','Lag k - filtered')
        plot(t_i,output.lagk_filtered-2.5*output.lagk_filtered_sd,'g--','DisplayName','Lag k variance')
        plot(t_i,output.lagk_filtered+2.5*output.lagk_filtered_sd,'g--')
        ylabel('Lag','FontWeight','bold','FontSize',12);
        ylim([0 1])
        %legend('show')
        
        
    end   



    %Check if we need to convert back to original unit
    if strcmp(unit1,'mg_dL')==1
        error('Todo Unit change not implemented yet')
    end
end%function


%Helper method to parse the arguments
function parsedArgs = parseInputVarArgs(varargs)
    %Set defaults
    parsedArgs.tout = [];
    parsedArgs.presmooth = true;
    parsedArgs.debugPlot = false;
    if mod(size(varargs,2),2)~=0
        error(['key/value arguments are not even:' num2str(size(varargs))])
    end
    for i = 1:2:size(varargs,2)
        if isfield(parsedArgs,varargs{i})
            %disp(['Parsing argument:' varargs{i}])
            parsedArgs.(varargs{i})=varargs{i+1};
        end
    end
end

function [Phi PhiEKF] = computePhiEKF(dynModel,x)
    ik = dynModel.Nin+3;%Index of the k_isf state
    k = x(ik);
    Gp = x(1);
    iGisf = dynModel.Nin+1;%Index of the Gisf state
    Gisf = x(iGisf);
    
    Phi = dynModel.Phi;
    Phi(iGisf,1)=k*dynModel.delta_t;
    Phi(iGisf,iGisf)=1-k*dynModel.delta_t;
    PhiEKF = Phi;
    PhiEKF(iGisf,end)=(Gp-Gisf)*dynModel.delta_t;
    
    
end

%Helper method to augment dynamic model with states needed for CGM dynamics
%estimation (add bias and lag)
function dynModel=augmentDynamicModelBiasAndLag(dynModelIn)
    dynModel.k_default = 1/5; %k = 1/T_isf (inverse lag for numerical stability);
    dynModel.Nin = size(dynModelIn.F,1);
    dynModel.Naug = 3;
    
    dynModel.stateNames = dynModelIn.stateNames;
    dynModel.stateNames{dynModel.Nin+1} = 'G_{isf}';
    dynModel.stateNames{dynModel.Nin+2} = 'b_cgm';
    dynModel.stateNames{dynModel.Nin+3} = 'k_isf';
    
    Ntot = dynModel.Naug+dynModel.Nin;
    dynModel.delta_t = dynModelIn.delta_t;
    dynModel.id = 2;
    dynModel.nonLinear = true;
    dynModel.F = zeros(Ntot,Ntot);               % System matrix (continuous)
    dynModel.F(1:dynModel.Nin,1:dynModel.Nin) = dynModelIn.F;
    %T_isf and bias terms are not added in system matrix, assumed to be non-varying
    dynModel.F(dynModel.Nin+1,1) = dynModel.k_default;
    dynModel.F(dynModel.Nin+1,dynModel.Nin+1) = -dynModel.k_default;
    
    dynModel.Q = zeros(Ntot,Ntot);               % Process noise covariance matrix.
    dynModel.Q(1:dynModel.Nin,1:dynModel.Nin) = dynModelIn.Q;     
    %No more noise terms added for the augmentation, this assumes bias and lag is
    %modeled as an unknown constant,  and the G_isf is not affected by random noise directly.
    
    dynModel.H_both=[dynModelIn.H 0 0 0 ;
                     zeros(1,dynModel.Nin) 1 1 0];      % Measurement matrix for simultaneous fingerpricks and cgm measurements
    dynModel.H_fp=[dynModelIn.H 0 0 0];        % Measurement matrix for fingerprick only
    dynModel.H_cgm=[zeros(1,dynModel.Nin) 1 1 0];      % Measurement matrix for cgm measurements only
    
    if ~rank(obsv(dynModel.F,dynModel.H_both))==size(dynModel.F,1)
        error('Augmented system is not observable with both measurments')
    end
    if ~rank(obsv(dynModel.F,dynModel.H_fp))==size(dynModel.F,1)
        disp('Augmented system is not observable with only fp measurements')
    end
    if ~rank(obsv(dynModel.F,dynModel.H_cgm))==size(dynModel.F,1)
        disp('Augmented system is not observable with only cgm measurements')
    end
    dynModel.initCov = dynModelIn.initCov;   % Initial covariance
    dynModel.initCov(dynModel.Nin+1,dynModel.Nin+1) = 0.25;   % Initial covariance for the G_isf state
    dynModel.initCov(dynModel.Nin+2,dynModel.Nin+2) = 1;   % Initial covariance for the bias. Biases seem to range from -2 to 2 in my data. Initial guessing based on the data could 
    dynModel.initCov(dynModel.Nin+3,dynModel.Nin+3) = 0.01;    % Initial covariance for the inverse lag. Rougly based on an a priori assumed range of 1 to 30 minutes value for T_isf (=1/k_isf).
    dynModel.initCov(1,1) = 0.25; %Since we use init we assume low glucose 
    
    %%% Discretization
    dynModel.Phi=eye(size(dynModel.F))+dynModel.F*dynModel.delta_t; % Discrete state transition matrix
    

    
    
end








