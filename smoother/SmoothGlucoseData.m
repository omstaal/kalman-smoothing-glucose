function output=SmoothGlucoseData(t_in,y_in,varargin)
% SMOOTHGLUCOSEDATA Creates a smoothed glucose curve with variance
%                   estimates from input glucose readings
% Usage:
% output=SmoothGlucoseData(t,y,t_i,...)
%   Generates a smoothed estimate from the input data t (array of datetime, 
%   or time in minutes as doubles) and y (glucose values)
%
%   output is a struct with fields:
%       y_smoothed : smoother estimates, mean 
%       y_smoothed_sd : smoother estimates, standard deviation
%       y_filtered : forward pass KF estimates, mean 
%       y_filtered_sd : forward pass KF estimates, standard deviation
%       t_i : interpolated time corresponding to the above vectors
%       The above are all from first to last measurement with 1 sec resolution  
%
%   
%   The supported variable arguments are as follows:
%   'y_error' : [] or an array of same length as y
%   'outlierRemoval' : 0,1 or 2
%   'outlierSDlimit : a number > 0
%   'plotResult' : 0, 1 or 2
%   'plotInternalStates' : 0 or 1
%   'startDateTime' : a datetime
%   'dynamicModel' : 1 or 2
%   't_out' : user-supplied vector of times that an estimate is wanted for, 
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
%   use when removing outliers. 2.5 would be a conservative value while 1
%   would lead to quite aggressive outlier removal
%   Default if not supplied is 2.5 
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

%TODO list for the future
% 1) possibility to specify "no meal/insulin input" periods in the input
% data, where a lower process noise should be used.

%Parse the variable arguments
parsedArgs = parseInputVarArgs(varargin);

%Handle unit
if strcmp(parsedArgs.unit,'auto')==1
   parsedArgs.unit = autoDetectUnit(y_in);
end
if strcmp(parsedArgs.unit,'mgdl')==1 %This code assumes mmol/L, so convert to that and convert back at the end
    y_in = convertTo_mg_dL(y_in);
    parsedArgs.y_error = convertTo_mmol_L(parsedArgs.y_error);
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

%Make input vectors dense (Remove any nan entries in y)
nonNan = ~isnan(y_in);
y = y_in(nonNan);
t = t(nonNan);

delta_t = 1/6;	   %Stepping more often than every second is not recommended, nor more seldom than 1 minute
t_i = t(1):delta_t:t(end);

%Set dynamic model to use
dynModel = setDynamicModel(parsedArgs.dynamicModel, delta_t);
Nstates = length(dynModel.H);

%Set up error to variance computation
sdsInConfInterval = 2.5; %2 for 95% CI, 2.5 for 99% CI
error2var = @(error) (error/sdsInConfInterval)^2; %Assumes normal distribution, and  error is given as a 99% confidence interval

if length(parsedArgs.y_error)==length(y)   % Assume user supplied error estimates
     y_error = parsedArgs.y_error;
elseif isempty(parsedArgs.y_error)         % Empty array supplied, assume error follows ISO15197
    y_error = setIsoError(y);
else
    error('Bad y_array argument supplied')    
end

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
            dz = y(l)-dynModel.H*xBar;
            R = error2var(y_error(l));
            Pz = (dynModel.H*PBar*dynModel.H'+R);
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
if strcmp(parsedArgs.unit,'mgdl')==0 
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
    parsedArgs.y_error = [];
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
    if(dynModelNo==1) %simple 2.order system where the rate of change of glucose dies out.
        a=-0.025;
        F =[0 1;0 a];               % System matrix (continuous)
        dynModel.Q=[0 0;0 0.02*delta_t];     % Process noise covariance matrix.
        dynModel.H=[1 0];                    % Measurement matrix.
        dynModel.initCov = diag([0.25 1]);   % Initial covariance
        %%% Discretization
        dynModel.Phi=expm(F*delta_t);        % Discrete state transition matrix
        dynModel.stateNames = {'Gp','dGp'};
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
    for i=2:length(datetimes)
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

%Helper function to find interpolted values at specific times
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
        unit = 'mg/dL';
    end
end
