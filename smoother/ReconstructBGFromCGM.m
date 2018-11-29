function output = ReconstructBGFromCGM( t_in,y_in_cgm, bias, lag,varargin )
%RECONSTRUCTBGFROMCGM Creates a smoothed SMBG blood glucose curve from CGM glucose readings,
% and an estimated bias and lag constant
% Usage:
% output=ReconstructBGFromCGM(t,y,t_i,...)
%   Generates a smoothed estimate from the input data t (array of datetime, 
%   or time in minutes as doubles) and y (CGM glucose values)
%
%   output is a struct with fields:
%       y_smoothed : smoother estimates, mean 
%       y_smoothed_sd : smoother estimates, standard deviation
%       y_filtered : forward pass KF estimates, mean 
%       y_filtered_sd : forward pass KF estimates, standard deviation
%       t_i : interpolated time corresponding to the above vectors
%       The above are all from first to last measurement with 10 sec
%       resolution
%       If a tout argument has been supplied, the interpolated values at
%       those times can be retrieved in the fields
%       y_smoothed_at_tout and y_smoothed_sd_at_tout
%       
%
%   
%   The supported variable arguments are as follows:
%   'y_error' : [] or an array of same length as y
%   'outlierRemoval' : 0,1 or 2
%   'outlierSDlimit : a number > 0
%   'plotResult' : 0, 1 or 2
%   'plotInternalStates' : 0 or 1
%   'dynamicModel' : 1,2 or 3
%   'tout' : user-supplied vector of times that an estimate is wanted for, 
%            either relative time or array of datetimes
%   'unit'  : string decribing which glucose unit the y data is given in
%
%   y_error can be set to:
%    - [], an empty array, to signify that limits from ISO 15197 should be used
%    - a user supplied array of same length as y, that has the errors for
%    individual measurements in y
%   Default if not supplied is []
%   
%   The outlierRemoval parameter controls outlier removal:
%   Measurements are considered to be outliers if the innovation is outside
%   X std devs of the innovation variance in the forward pass Kalmanfilter
%   If outlierRemoval==1, outlier removal is performed based on the smoothed estimate.
%   If outlierRemoval==2, outlier removal is performed based on the filtered estimate.
%   If outlierRemoval is any other value, outlier removal is not done
%   Note that the smoothing approach inherently suppresses outliers if
%   surrounding data allows
%   Default if not supplied is 0 (no outlier removal, only suppression)
%
%   outlierSDlimit is a double specifying how many standard deviations of the estimate to use
%   use when removing outliers. 2 would be a conservative value while 1
%   would lead to quite aggressive outlier removal
%   Default if not supplied is 2 
%
%   If plotResult==1, a plot will be produced showing the result of the
%   smoothing in a new figure. If plotResult==2, the estimates from the forward pass KF
%   will be added to the same plot.
%   Default if not supplied is 0
%
%   if plotInternalStates==1, the internal states of the dynamic model in
%   use will be plotted in a new figure.
%   Default if not supplied is 0
%
%   startDateTime is a datetime defining the start 
%   of the experiment. Allows for plotting the experiment with the datetime format.
%   If not defined, the plot will start at 0 [min].
%   Default if not supplied is 0
%   
%   if dynamicModel == 1 a simple 2-state model is used, that only
%   describes the glucose rate of change
%   if dynamicModel == 2 a 3-state model is used, that
%   describes the glucose rate of change being affected by a remote
%   compartment, which is fed by a central compartment (insulin and glucose
%   is lumped together to one state in this model)
%   Default if not supplied is 1
%
%   If unit is set to 'mg_dl' the input data (y, y_error) is assumed to be in mg/dL
%   If set to 'mmol_l' the input data is assumed to be in mmol/L
%   If set to 'auto', a autodetection routine is run on y that guesses
%   which unit is used based on the values found
%   Default if not supplied is 'auto'
% This code has been tested on Matlab R2016b. There are issues with using
% this code on earlier versions.

%Parse the variable arguments
    parsedArgs = parseInputVarArgs(varargin);
    
    %Check units
    unit = autoDetectGlucoseUnit(y_in_cgm);
    if strcmp(unit,'mmol_L')==0
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

    %Run smoothing on each set individually
    output.smoothed_cgm = SmoothSMBGData(t_in,y_in_cgm,'outlierRemoval',parsedArgs.outlierRemoval,'dynamicModel',parsedArgs.dynamicModel);
    t_i = output.smoothed_cgm.t_i;
    output.t_i = t_i;
    y_cgm = output.smoothed_cgm.y_smoothed;
    var_cgm = output.smoothed_cgm.y_smoothed_sd.^2;
    %Set dynamic model to use for CGM SMBG fusion
    dynModel = augmentDynamicModelKnownLag(output.smoothed_cgm.dynModel, lag);
    Nstates = size(dynModel.Phi,1);
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
            if ~isnan(y_cgm(k)) %We can initialise
                %%% Initialization
                xBar = zeros(Nstates,1);
                xBar(1)=y_cgm(k)-bias;
                xBar(dynModel.Nin+1)=y_cgm(k)-bias; 
    
                    
                xHat=xBar;
                PBar=dynModel.initCov;
                
                PHat=PBar;
                init = true;
                startk = k;

                H = dynModel.H;
                R = var_cgm(k);

            end
        elseif init && isnan(y_cgm(k)) %time to end
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
            y = y_cgm(k)-bias;
             
            dz = y-H*xBar;
            Pz = (H*PBar*H'+R);

            %Measurement update
            K=PBar*H'/Pz;
            xHat = xBar + K*dz;
            %PHat = (eye(size(PBar))-K*H)*PBar;
            hlp = (eye(size(PBar))-K*H);
            PHat = hlp*PBar*hlp' + K*R*K';
            
            
             [~,p] = chol(PHat);
             if p~=0
                 error(sprintf('Phat not pos def at k=%d',k))
             end
            if ~issymmetric(PHat)
                PHat = (PHat+PHat')/2; %Make symmetric
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
    output.abortedSmoothing = false;
    PhiForRTS = dynModel.Phi;
    for k = endk-1:-1:startk+1
        if dynModel.nonLinear
            PhiForRTS = PhiEKFs(:,:,k);
        end
        C=(P_hat_f(:,:,k)*PhiForRTS')/P_bar_f(:,:,k+1);
        x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
        
        P_hat_s=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';  
        [~,p] = chol(P_hat_s);
        if p~=0
            [l m] = eig(P_hat_s)
            disp(sprintf('Warning - smoothing aborted at t = %f due to non-posdef covmatrix. time remaining to smooth: %f',t_i(k),t_i(k)-t_i(startk)))
            output.abortedSmoothing = true;
            break
        end
        if ~issymmetric(P_hat_s)
            P_hat_sym = (P_hat_s+P_hat_s')/2; %Make symmetric
            err = max(P_hat_s - P_hat_sym);
            %disp(['Unsymmetric P_hat_s ,k=' num2str(k) ', err=' err])
            P_hat_s = P_hat_sym;
        end
        P_smoothed(:,:,k) = P_hat_s;
        lastSmoothK = k;
    end

    %Generate output structs
    output.y_fprec_smoothed = nan(size(t_i));
    output.y_fprec_smoothed_sd = nan(size(t_i));
    output.y_cgm_smoothed = nan(size(t_i));
    output.y_cgm_smoothed_sd = nan(size(t_i));
    
    output.y_fprec_smoothed(startk:endk) = x_smoothed(1,startk:endk);
    output.y_cgm_smoothed(startk:endk) = x_smoothed(dynModel.Nin+1,startk:endk);
    
    for k = startk:endk
        output.y_fprec_smoothed_sd(k) = sqrt(P_smoothed(1,1,k));
        output.y_cgm_smoothed_sd(k) = sqrt(P_smoothed(dynModel.Nin+1,dynModel.Nin+1,k));
    end 
    
    %Add internal states
    output.x_filtered = x_hat_f;
    output.x_smoothed = x_smoothed;

    %Add user supplied wanted time
    output.y_fprec_at_tout = closestValues(parsedArgs.tout,t_i,output.y_fprec_smoothed, startDateTime);
    output.y_fprec_sd_at_tout = closestValues(parsedArgs.tout,t_i,output.y_fprec_smoothed_sd, startDateTime);

     %debug plotting
    if parsedArgs.debugPlot
        sdsInConfInterval=2;
        %figure()
        plot(t_in,y_in_cgm,'b.','DisplayName','FGM measurements', 'MarkerSize',15)
        hold on;
        plot(t_i,y_cgm,'b-','DisplayName','Smoothed FGM','LineWidth',2)
        plot(t_i, output.y_fprec_smoothed,'k-','LineWidth',2,'DisplayName','Bias + lag corrected FGM')
        if length(parsedArgs.debugPlotTrueYfp)>0
            plot(t_in, parsedArgs.debugPlotTrueYfp,'r.','DisplayName','True SMBG','MarkerSize',15)
        end 
        hold off;
        xlabel('Time [min]');
        ylabel('Glucose [mmol/L]');
        title(['FGM corrected for bias ' num2str(bias,'%.1f') ' and lag ' num2str(lag,'%.1f')])
       % legend('Location','northwest')
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
   
    end   



    %Check if we need to convert back to original unit
    if strcmp(unit,'mg_dL')==1
        error('Todo Unit change not implemented yet')
    end
end%function


%Helper method to parse the arguments
function parsedArgs = parseInputVarArgs(varargs)
    %Set defaults
    parsedArgs.tout = [];
    parsedArgs.debugPlot = false;
    parsedArgs.debugPlotTrueYfp = [];
    parsedArgs.dynamicModel = 2;
    parsedArgs.outlierRemoval = 1;
    if mod(size(varargs,2),2)~=0
        error(['key/value arguments are not even:' num2str(size(varargs))])
    end
    for i = 1:2:size(varargs,2)
        if isfield(parsedArgs,varargs{i})
            %disp(['Parsing argument:' varargs{i}])
            parsedArgs.(varargs{i})=varargs{i+1};
        else
            disp(['Parsed unknown argument:' varargs{i}])
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
function dynModel=augmentDynamicModelKnownLag(dynModelIn, lag)
    dynModel.k_default = 1/5; %k = 1/T_isf (inverse lag for numerical stability);
    dynModel.Nin = size(dynModelIn.F,1);
    dynModel.Naug = 1;
    
    dynModel.stateNames = dynModelIn.stateNames;
    dynModel.stateNames{dynModel.Nin+1} = 'G_{isf}';
    
    Ntot = dynModel.Naug+dynModel.Nin;
    dynModel.delta_t = dynModelIn.delta_t;
    dynModel.id = 2;
    dynModel.nonLinear = false;
    dynModel.F = zeros(Ntot,Ntot);               % System matrix (continuous)
    dynModel.F(1:dynModel.Nin,1:dynModel.Nin) = dynModelIn.F;
    %T_isf and bias terms are not added in system matrix, assumed to be non-varying
    dynModel.F(dynModel.Nin+1,1) = 1/lag;
    dynModel.F(dynModel.Nin+1,dynModel.Nin+1) = -1/lag;
    
    dynModel.Q = zeros(Ntot,Ntot);               % Process noise covariance matrix.
    dynModel.Q(1:dynModel.Nin,1:dynModel.Nin) = dynModelIn.Q;     
    useParamProcessNoise=true;
    
    dynModel.H=[zeros(1,dynModel.Nin) 1];      % Measurement matrix for cgm measurements only
    
    if ~rank(obsv(dynModel.F,dynModel.H))==size(dynModel.F,1)
        disp('Augmented system is not observable with only cgm measurements')
    end
    dynModel.initCov = dynModelIn.initCov;   % Initial covariance
    dynModel.initCov(dynModel.Nin+1,dynModel.Nin+1) = 0.25;   % Initial covariance for the G_isf state
    dynModel.initCov(1,1) = 0.25; %Since we use init we assume low glucose covariance
    
    %%% Discretization
    dynModel.Phi=eye(size(dynModel.F))+dynModel.F*dynModel.delta_t; % Discrete state transition matrix
    dynModel.PhiEKF = dynModel.Phi;
end

