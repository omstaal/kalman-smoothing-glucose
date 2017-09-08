function output=SmoothSMBGandCGMData(t_in,y_in_fp,y_in_cgm,varargin)

    %Parse the variable arguments
    parsedArgs = parseInputVarArgs(varargin);
    
    %Check units
    unit1 = autoDetectGlucoseUnit(y_in_fp);
    unit2 = autoDetectGlucoseUnit(y_in_cgm);
    if strcmp(unit1,unit2) == 0
       error('Supplied units are not consistent for fingerpricks and CGM')
    end

    %TODO handle other units than mmol/L

    %Handle time
    if isdatetime(t_in)
        %convert to relative time
        startDateTime = t_in(1);
        t_in=convertToRelativeTime(t_in, startDateTime);
    else
        %Assume relative time is passed in
        startDateTime = NaN;
    end

    %Run smoothing on each set individually
    output.smoothed_fp = SmoothSMBGData(t_in,y_in_fp);
    output.smoothed_cgm = SmoothSMBGData(t_in,y_in_cgm);

    t_i = output.smoothed_fp.t_i;
    y_fp = output.smoothed_fp.y_smoothed;
    var_fp = output.smoothed_fp.y_smoothed_sd.^2;
    y_cgm = output.smoothed_cgm.y_smoothed;
    var_cgm = output.smoothed_cgm.y_smoothed_sd.^2;
    
    %Set dynamic model to use for CGM SMBG fusion
    dynModel = augmentDynamicModel(output.smoothed_fp.dynModel);
    Nstates = size(dynModel.Phi,1);

    %%% Storage
    x_hat_f = nan(Nstates,length(t_i));         % A priori state vector storage, forward pass
    x_bar_f = nan(Nstates,length(t_i));         % A posteriori state vector storage, forward pass
    P_hat_f = nan(Nstates,Nstates,length(t_i));       % A priori covariance matrix storage, forward pass
    P_bar_f = nan(Nstates,Nstates,length(t_i));       % A posteriori covariance matrix storage, forward pass
    x_smoothed = nan(Nstates,length(t_i));            % State vector storage, backward pass
    P_smoothed = nan(Nstates,Nstates,length(t_i));    % Covariance matrix storage, backward pass
    init = false;
    endk = nan;
    for k = 1:length(t_i)
        if ~init
            if ~isnan(y_fp(k)) || ~isnan(y_cgm(k)) %We can initialise
                %%% Initialization
                xBar = zeros(Nstates,1);
                if ~isnan(y_fp(k)) && ~isnan(y_cgm(k))
                    xBar(1)=y_fp(k);
                    xBar(end)=y_cgm(k);
                elseif ~isnan(y_fp(k))
                    xBar(1)=y_fp(k);
                    xBar(end)=y_fp(k);
                elseif ~isnan(y_cgm(k))
                    xBar(1)=y_cgm(k);
                    xBar(end)=y_cgm(k);
                end

                xHat=xBar;
                PBar=dynModel.initCov;
                PHat=PBar;
                init = true;
                startk = k;
            end
        elseif init &&  isnan(y_fp(k)) && isnan(y_cgm(k)) %time to end
            endk = k-1;
            break
        else %We have initialised and it is not time to end yet, do filtering
            %%% Kalman filter forward pass
            %TU - Time update
            xBar = dynModel.Phi*xHat;
            PBar = dynModel.Phi*PHat*dynModel.Phi' + dynModel.Q;
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
        C=P_hat_f(:,:,k)*dynModel.Phi'*inv(P_bar_f(:,:,k+1));
        x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
        P_smoothed(:,:,k)=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';  
    end

    %Generate output struct
    output.y_fp_smoothed = nan(size(t_i));
    output.y_fp_smoothed_sd = nan(size(t_i));
    output.y_cgm_smoothed = nan(size(t_i));
    output.y_cgm_smoothed_sd = nan(size(t_i));
    output.bias_smoothed = nan(size(t_i));
    output.bias_smoothed_sd = nan(size(t_i));
    
    output.bias_filtered = nan(size(t_i));
    output.bias_filtered_sd = nan(size(t_i));

    output.y_fp_smoothed(startk:endk) = x_smoothed(1,startk:endk);
    output.y_cgm_smoothed(startk:endk) = x_smoothed(Nstates,startk:endk);
    output.bias_smoothed(startk:endk) = x_smoothed(Nstates-1,startk:endk);
    output.bias_filtered(startk:endk) = x_hat_f(Nstates-1,startk:endk);

    for k = startk:endk
        output.y_fp_smoothed_sd(k) = sqrt(P_smoothed(1,1,k));
        output.y_cgm_smoothed_sd(k) = sqrt(P_smoothed(Nstates,Nstates,k));
        output.bias_smoothed_sd(k) = sqrt(P_smoothed(Nstates-1,Nstates-1,k));
        output.bias_filtered_sd(k) = sqrt(P_hat_f(Nstates-1,Nstates-1,k));
        
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
    if false
        subplot(2,1,1)
        plot(t_in,y_in_fp,'r.','DisplayName','Fingerpricks')
        hold on;
        plot(t_in,y_in_cgm,'b.','DisplayName','CGM measurements')
        plot(t_i,y_fp,'r-','DisplayName','Smoothed fingerpricks')
        plot(t_i,y_cgm,'b-','DisplayName','Smoothed CGM')
        plot(t_i,output.y_fp_smoothed,'g-','DisplayName','Merge-smoothed FPs')
        plot(t_i,output.y_cgm_smoothed,'k-','DisplayName','Merge-smoothed CGM')
        
        hold off;
        xlabel('Time [min]','FontWeight','bold','FontSize',12);
        ylabel('Glucose [mmol/L]','FontWeight','bold','FontSize',12);
        title('Smoothing CGM and SMBG')
        subplot(2,1,2)
        plot(t_i,output.bias_smoothed,'r-','DisplayName','Bias - smoothed')
        hold on;
        plot(t_i,output.bias_smoothed-2.5*output.bias_smoothed_sd,'r--','DisplayName','Bias variance')
        plot(t_i,output.bias_smoothed+2.5*output.bias_smoothed_sd,'r--')
        plot(t_i,output.bias_filtered,'g-','DisplayName','Bias - filtered')
        plot(t_i,output.bias_filtered-2.5*output.bias_filtered_sd,'g--','DisplayName','Bias variance')
        plot(t_i,output.bias_filtered+2.5*output.bias_filtered_sd,'g--')
        ylabel('Bias','FontWeight','bold','FontSize',12);
       
        ylim([-3 3])
    end   



    %Check if we need to convert back to original unit
    if strcmp(unit1,'mg_dL')==1
        error('Todo Not implemented yet')
    end
end%function


%Helper method to parse the arguments
function parsedArgs = parseInputVarArgs(varargs)
    %Set defaults
    parsedArgs.tout = [];
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

%Helper method to augment dynamic model with states needed for CGM dynamics
%estimation (add bias and 
function dynModel=augmentDynamicModel(dynModelIn)
    T_isf = 20;
    Nin = size(dynModelIn.F,1);
    Naug = 2;
    Ntot = Naug+Nin;
    
    dynModel.F = zeros(Ntot,Ntot);               % System matrix (continuous)
    dynModel.F(1:Nin,1:Nin) = dynModelIn.F;
    %bias term is not added in system matric, assumed non-varying
    dynModel.F(Nin+2,1) = 1/T_isf;
    dynModel.F(Nin+2,end) = -1/T_isf;
    
    dynModel.Q = zeros(Ntot,Ntot);               % Process noise covariance matrix.
    dynModel.Q(1:Nin,1:Nin) = dynModelIn.Q;     
    %No more noise terms added for the augmentation, this assumes bias is
    %modeled as an unknown constant,  and the G_isf is not affected by ransom noise directly.
    
    dynModel.H_both=[dynModelIn.H 0 0;
                     zeros(1,Nin) 1 1];      % Measurement matrix for simultaneous fingerpricks and cgm measurements
    dynModel.H_fp=[dynModelIn.H 0 0];        % Measurement matrix for fingerprick only
    dynModel.H_cgm=[zeros(1,Nin) 1 1];      % Measurement matrix for cgm measurements only
    
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
    dynModel.initCov(Nin+1,Nin+1) = 1;   % Initial covariance for the bias. Based on 39 sets of freestyle libre data
    dynModel.initCov(Nin+2,Nin+2) = 0.25;   % Initial covariance for the G_isf state
    
    %%% Discretization
    dynModel.Phi=expm(dynModel.F*dynModelIn.delta_t); % Discrete state transition matrix
    
    dynModel.stateNames = dynModelIn.stateNames;
    dynModel.stateNames{Nin+1} = 'cgm bias';
    dynModel.stateNames{Nin+2} = 'G_{isf}';
    
end










