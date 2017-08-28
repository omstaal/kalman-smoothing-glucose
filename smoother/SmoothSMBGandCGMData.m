function output=SmoothSMBGandCGMData(t_in,y_in_fp,y_in_cgm,varargin)
%Parse the variable arguments
parsedArgs = parseInputVarArgs(varargin);

%Handle unit
if strcmp(parsedArgs.unit,'auto')==1
   parsedArgs.unit = autoDetectUnit(y_in_fp);
   unit2 = autoDetectUnit(y_in_cgm);
   if strcmp(parsedArgs.unit,unit2) == 0
       error('Supplied units are not consistent for fingerpricks and CGM')
   end
end
if strcmp(parsedArgs.unit,'mg_dL')==1 %This code assumes mmol/L, so convert to that and convert back at the end
    y_in_fp = convertTo_mmol_L(y_in_fp);
    y_in_cgm = convertTo_mmol_L(y_in_cgm);
    
    parsedArgs.y_error_fp = convertTo_mmol_L(parsedArgs.y_error_fp);
    parsedArgs.y_error_cgm = convertTo_mmol_L(parsedArgs.y_error_cgm);
    
    %TODO need to handle user supplied fp error and cgm
end

%Handle time
if isdatetime(t_in)
    %convert to relative time
    t=convertToRelativeTime(t_in, t_in(1));
    %override startDateTime
    if(parsedArgs.startDateTime~=0)
        disp('Warning! Overriding startDateTime since t_in was a datetime array')
    end
    parsedArgs.startDateTime = t_in(1);
else
    %Assume relative time is passed in
    t=t_in;
end

%Run smoothing on each set individually
smoothed_fp = SmoothSMGBData(t,y_fp,'y_error',parsedArgs.y_error_fp)
smoothed_cgm = SmoothSMGBData(t,y_cgm,'y_error',parsedArgs.y_error_cgm)

%Set dynamic model to use
dynModel = setDynamicModel(parsedArgs.dynamicModel, delta_t);
Nstates = size(dynModel.Phi,1);

%Set up error to variance computation
sdsInConfInterval = 2.5; %2 for 95% CI, 2.5 for 99% CI
error2var = @(error) (error/sdsInConfInterval)^2; %Assumes normal distribution, and  error is given as a 99% confidence interval

output.outliers = false(size(y));
doneFindingOutliers=false;

while ~doneFindingOutliers;
    %%% Storage
    x_hat_f = zeros(Nstates,length(t_i));         % A priori state vector storage, forward pass
    x_bar_f = zeros(Nstates,length(t_i));         % A posteriori state vector storage, forward pass
    P_hat_f = zeros(Nstates,Nstates,length(t_i));       % A priori covariance matrix storage, forward pass
    P_bar_f = zeros(Nstates,Nstates,length(t_i));       % A posteriori covariance matrix storage, forward pass
    x_smoothed = zeros(Nstates,length(t_i));            % State vector storage, backward pass
    P_smoothed = zeros(Nstates,Nstates,length(t_i));    % Covariance matrix storage, backward pass

    %%% Initialization
    xBar = zeros(Nstates,1);
    xBar(1)=y(1);
    xHat=xBar;
    PBar=dynModel.initCov;
    PHat=PBar;
    l=1;
    %%% Kalman filter forward pass
    for k = 1:length(t_i)
        %TU - Time update
        xBar = dynModel.Phi*xHat;
        PBar = dynModel.Phi*PHat*dynModel.Phi' + dynModel.Q;
        %Store
        x_bar_f(:,k)=xBar;
        P_bar_f(:,:,k)=PBar;

        measUpdateDone=0;
        %MU - Measurement Update only when we have a measurement
        while length(t)>=l && t_i(k)>=t(l)  % Interpolated time has passed one 
                                            % of the measurement times, process
                                            % all measurements that has occurred
            if measUpdateDone==1
                % More than one measurement at the current time
                xBar = xHat;
                PBar = PHat;
            end
            %Switch between measurements
            if ~isnan(y_fp(l))
                H = dynmodel.H_fp;
                y = y_fp;
            elseif ~isnan(y_cgm(l))
                H = dynmodel.H_cgm;
                y = y_cgm;
            else
                error('Encountered empty measurement, should have been filtered out')
            end
            %TODO handle that we have both fp and cgm in the same instant
            dz = y(l)-H*xBar;
            R = error2var(y_error(l));
            Pz = (H*PBar*H'+R);
            if parsedArgs.outlierRemoval == 2 
                % Check the innovation
                if(abs(dz)>parsedArgs.outlierSDlimit*sqrt(Pz))
                    output.outliers(l)=true;
                    disp(['Forward pass flagged measurement as outlier: t = ' num2str(t(l)) ' [min], y = ' num2str(y(l)) ' [mmol/L].'])
                end
            end
            if ~output.outliers(l)
                %Measurement update
                K=PBar*dynModel.H'/Pz;
                xHat = xBar + K*dz;
                PHat = (eye(size(PBar))-K*dynModel.H)*PBar;
                measUpdateDone=1;
            end
            l=l+1;
                
        end
        if measUpdateDone==0    % No measurement was available at this time 
            xHat=xBar;
            PHat=PBar;
        end
        % Store
        x_hat_f(:,k)=xHat;
        P_hat_f(:,:,k)=PHat;
    end % for k

    %%% Rauch-Tung-Striebel backward pass
    x_smoothed(:,k)=xHat;
    P_smoothed(:,:,k)=PHat;
    for k = length(t_i)-1:-1:1
        C=P_hat_f(:,:,k)*dynModel.Phi'*inv(P_bar_f(:,:,k+1));
        x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
        P_smoothed(:,:,k)=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';  
    end
    
    %Generate output struct
    output.y_smoothed = x_smoothed(1,:);
    output.y_smoothed_sd = zeros(size(output.y_smoothed));
    for k = 1:length(t_i)
        output.y_smoothed_sd(k) = sqrt(P_smoothed(1,1,k));
    end
    
    if parsedArgs.outlierRemoval == 1
        %Run through all measurements and see if any are outside the error
        %smoothed band, if so they are outliers
        foundNewOutliers = false;
        y_s_mean  = interpolatedValuesAt(t,t_i,output.y_smoothed, parsedArgs.startDateTime);
        y_s_sd  = interpolatedValuesAt(t,t_i,output.y_smoothed_sd, parsedArgs.startDateTime);
        for i = 1:length(y)
            if ~output.outliers(i) && abs(y(i)-y_s_mean(i))/y_s_sd(i)>parsedArgs.outlierSDlimit
               output.outliers(i)=true;
               foundNewOutliers = true;
               disp(['Smoother flagged measurement ' num2str(i) ' as outlier: t = ' num2str(t(i)) ' [min], y = ' num2str(y(i)) ' [mmol/L].']) 
            end
        end
        if ~foundNewOutliers
            doneFindingOutliers = true;
        else
            disp(['Smoother needs a second pass due to outliers detected. Total # outliers in input data:' num2str(sum(output.outliers))]) 
            
        end
    else
        doneFindingOutliers = true;
    end
end %while not doneFindingOutliers
%Smoothing done



output.y_filtered = x_hat_f(1,:);
output.y_filtered_sd = zeros(size(output.y_filtered));
for k = 1:length(t_i)
    output.y_filtered_sd(k) = sqrt(P_hat_f(1,1,k));
end
output.t_i = t_i;

output.x_filtered = x_hat_f;
output.x_smoothed = x_smoothed;

%Add user supplied wanted time
output.y_smoothed_at_tout = interpolatedValuesAt(parsedArgs.tout,t_i,output.y_smoothed, parsedArgs.startDateTime);
output.y_smoothed_sd_at_tout = interpolatedValuesAt(parsedArgs.tout,t_i,output.y_smoothed_sd, parsedArgs.startDateTime);


%Plot if specified
if isdatetime(parsedArgs.startDateTime)
    t_plot = parsedArgs.startDateTime + minutes(t); % Convert to datetime for plotting
    t_i_plot = parsedArgs.startDateTime + minutes(t_i); % Convert to datetime for plotting
else
    t_plot = t;
    t_i_plot = t_i;
end
  

if parsedArgs.plotResult>=1
    figure()
    
    %h1 = errorbar(t_plot,y,y_error,'r.','MarkerSize',10,'LineWidth',1); % Using 2 sigma
    h1 = plot(t_plot,y,'r.','MarkerSize',10); % Using 2 sigma
    
    hold on
    h2 = plot(t_i_plot,output.y_smoothed,'b','LineWidth',2);
    h3 = plot(t_i_plot,output.y_smoothed+sdsInConfInterval*output.y_smoothed_sd,'b--');
    plot(t_i_plot,output.y_smoothed-sdsInConfInterval*output.y_smoothed_sd,'b--')

    xlabel('Time [min]','FontWeight','bold','FontSize',12);
    ylabel('Glucose [mmol/L]','FontWeight','bold','FontSize',12);
    if parsedArgs.plotResult>=2
        h4= plot(t_i_plot,output.y_filtered,'g','LineWidth',2);
        h5 = plot(t_i_plot,output.y_filtered+sdsInConfInterval*output.y_filtered_sd,'g--');
        plot(t_i_plot,output.y_filtered-sdsInConfInterval*output.y_filtered_sd,'g--')
        legend([h1 h2 h3 h4 h5],{'Unfiltered glucose data','Smoothed estimate','\pm2.5 SD of smoothed estimate', 'Filtered estimate', '\pm2.5 SD of filtered estimate' })
    else
        legend([h1 h2 h3],{'Unfiltered glucose data','Smoothed estimate','\pm2.5 SD of smoothed estimate' })
    end
    title('Kalman smoothing','FontWeight','bold','FontSize',14)
    hold off
end

if parsedArgs.plotInternalStates==1
    figure()
    subplot(2,1,1)
    plot(t_i_plot,x_hat_f);
    legend(dynModel.stateNames)
    ylabel('Forward pass internal states')
    subplot(2,1,2)
    plot(t_i_plot,x_smoothed)
    legend(dynModel.stateNames)
    ylabel('Smoothed internal states')
end
%Check if we need to convert back to original unit
if strcmp(parsedArgs.unit,'mg_dL')==1
    output.y_filtered = convertTo_mg_dL(output.y_filtered);
    output.y_filtered_sd = convertTo_mg_dL(output.y_filtered_sd);
    output.y_smoothed = convertTo_mg_dL(output.y_smoothed);
    output.y_smoothed_sd = convertTo_mg_dL(output.y_smoothed_sd);
    output.y_smoothed_at_tout = convertTo_mg_dL(output.y_smoothed_at_tout);
    output.y_smoothed_sd_at_tout = convertTo_mg_dL(output.y_smoothed_sd_at_tout);
end
end%function SmoothGlucoseData

%Helper method to parse the arguments
function parsedArgs = parseInputVarArgs(varargs)
    %Set defaults
    parsedArgs.outlierRemoval = 0;
    parsedArgs.outlierSDlimit = 2.5;
    parsedArgs.plotResult = 0;
    parsedArgs.plotInternalStates = 0;
    parsedArgs.startDateTime = 0;
    parsedArgs.dynamicModel = 1;
    parsedArgs.y_error_fp = [];
    parsedArgs.y_error_cgm = [];
    parsedArgs.tout = [];
    parsedArgs.unit = 'auto';
    
    if mod(size(varargs,2),2)~=0
        error(['key/value arguments are not even:' num2str(size(varargs))])
    end
    for i = 1:2:size(varargs,2)
        if isfield(parsedArgs,varargs{i})
            %disp(['Parsing argument:' varargs{i}])
            parsedArgs.(varargs{i})=varargs{i+1};
        else
            error(['Invalid argument supplied:' num2str(varargs{i})])
        end
    end
end

%Helper method to set dynamic model
function dynModel=setDynamicModel(dynModelNo,delta_t)
    T_isf = 10;
    if(dynModelNo==1) %simple 2.order system where the rate of change of glucose dies out.
        a=-0.025;
        F =[0 1 0 0;
            0 a 0 0;
            1/Tisf 0 -1/Tisf 0;
            0 0 0 0];               % System matrix (continuous)
        dynModel.Q=[0 0 0 0;
                    0 0.02*delta_t 0 0;
                    0 0 0 0;
                    0 0 0 0];     % Process noise covariance matrix.
        dynModel.H_fp=[1 0 0 0];      % Measurement matrix for fingerpricks
        dynModel.H_fp=[0 0 1 1];      % Measurement matrix for CGM
        dynModel.initCov = diag([0.25 1 1 1]);   % Initial covariance
        %%% Discretization
        dynModel.Phi=expm(F*delta_t);        % Discrete state transition matrix
        dynModel.stateNames = {'Gp','dGp','Gisf','b'};
    elseif (dynModelNo==2)
        Td = 15.000; % Time constant describing flow between compartments [min]
        F =[0 0 1;0 -1/Td 0;0 1/Td -1/Td]; % System matrix (continuous)
        dynModel.Q=[0 0 0;0 0.02*delta_t 0;0 0 0]; % Process noise covariance matrix.
        dynModel.H=[1 0 0];                       % Measurement matrix.
        dynModel.initCov = diag([10 1 1]);         % Initial covariance
        dynModel.Phi=expm(F*delta_t);                    % Discrete state transition matrix
        dynModel.stateNames = {'Gp','C','R'};
    else
        error(['Unsupported model:' num2str(model)])
    end
end

%Helper method to convert datetimes to relative time
function t = convertToRelativeTime(datetimes,startdatetime)
    t=zeros(length(datetimes),1);
    for i=1:length(datetimes)
        t(i) = minutes(datetimes(i)-startdatetime);
    end
end

%Helper method to convert glucose unit from mmol/L to mg/dL
function y_mgdL = convertTo_mg_dL(y_mmolL)
    y_mgdL = y_mmolL*18.018;
end

%Helper method to convert glucose unit from mg/dL to mmol/dL
function y_mmolL = convertTo_mmol_L(y_mgdL)
    y_mmolL = y_mgdL/18.018;
end

%Helper method that guesses a error based on the measured value, based on ISO 15197 
function y_error = setIsoError(y)
    y_error = zeros(size(y));              % ISO 15197, set it based on the measured values
    for i = 1:length(y) % Make error bars (assumes fingerprick measurement
                        % errors according to ISO15197)
        if y(i)>5.6
            y_error(i) = 0.2*y(i);
        else
            y_error(i) = 0.83;
        end
    end
end

%Helper function to find closest to wanted output times
function yout = interpolatedValuesAt(tout, t, y, startdatetime)
    if isdatetime(tout)
        tout = convertToRelativeTime(tout, startdatetime);
    end
    yout = zeros(size(tout));
    for i = 1:length(tout)
        if tout(i)>t(end) || tout(i)<t(1)
            yout(i) = NaN;
        else
            for j=1:length(t)
                if t(j)>=tout(i)
                    yout(i) = y(j);
                    break;
                end
            end
        end
    end
end

%Helper function to autodetect glucose unit
function unit = autoDetectUnit(measurements)
    unit = 'mmol_L';
    if mean(measurements)>50
        disp('SmoothGlucoseData autodetected mg/dL as unit')
        unit = 'mg_dL';
    end
end
