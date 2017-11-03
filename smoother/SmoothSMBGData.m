function output=SmoothSMBGData(t_in,y_in,varargin)
% SMOOTHSMBGDATA Creates a smoothed glucose curve from input glucose readings,
% assumed to come from a Self Monitoring Blood Glucose meter
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
% This code has been tested on Matlab R2016b. There are issues with using
% this code on earlier versions.

%TODO list for the future
% 1) possibility to specify "no meal/insulin input" periods in the input
% data, where a lower process noise should be used.

%Parse the variable arguments
parsedArgs = parseInputVarArgs(varargin);

%Handle unit
if strcmp(parsedArgs.unit,'auto')==1
   parsedArgs.unit = autoDetectGlucoseUnit(y_in);
end
if strcmp(parsedArgs.unit,'mg_dL')==1 %This code assumes mmol/L, so convert to that and convert back at the end
    y_in = convertTo_mmol_L(y_in);
    parsedArgs.y_error = convertTo_mmol_L(parsedArgs.y_error);
end



%Handle time
if isdatetime(t_in)
    %convert to relative time
    t_in.TimeZone='';
    startDateTime = t_in(1);
    t_in=convertToRelativeTime(t_in, startDateTime);
else
    %Assume relative time is passed in
    startDateTime = NaN;
end

%Set dynamic model to use
dynModel = setDynamicModel(parsedArgs.dynamicModel);
output.delta_t = dynModel.delta_t;
Nstates = length(dynModel.H);

%Interpolated time vector
t_i = t_in(1):dynModel.delta_t:t_in(end);

%Make input vectors dense (Remove any nan entries in y)
nonNan = ~isnan(y_in);
y = y_in(nonNan);
t = t_in(nonNan);
t_i_first = find(t_i>=t(1), 1);
t_i_last = find(t_i<=t(end), 1, 'last');
if t_i_last<length(t_i) 
    t_i_last = t_i_last+1;
end
t_i_valid = t_i(t_i_first:t_i_last);
%disp(['Smoothing between ' num2str(t_i_valid(1)) ' and ' num2str(t_i_valid(end))])


%Set up error to variance computation
sdsInConfInterval = 2; %2 for 95% CI, 2.5 for 99% CI
error2var = @(error) (error/sdsInConfInterval)^2; %Assumes normal distribution, and  error is given as a 99% confidence interval

if length(parsedArgs.y_error)==length(y)   % Assume user supplied error estimates
     y_error = parsedArgs.y_error;
elseif isempty(parsedArgs.y_error)         % Empty array supplied, assume error follows ISO15197
    y_error = setIsoError(y);
else
    error('Bad y_error argument supplied')    
end

outliers = false(size(y));
doneFindingOutliers=false;

while ~doneFindingOutliers
    %%% Storage
    x_hat_f = zeros(Nstates,length(t_i_valid));         % A priori state vector storage, forward pass
    x_bar_f = zeros(Nstates,length(t_i_valid));         % A posteriori state vector storage, forward pass
    P_hat_f = zeros(Nstates,Nstates,length(t_i_valid));       % A priori covariance matrix storage, forward pass
    P_bar_f = zeros(Nstates,Nstates,length(t_i_valid));       % A posteriori covariance matrix storage, forward pass
    x_smoothed = zeros(Nstates,length(t_i_valid));            % State vector storage, backward pass
    P_smoothed = zeros(Nstates,Nstates,length(t_i_valid));    % Covariance matrix storage, backward pass

    %%% Initialization
    xBar = zeros(Nstates,1);
    xBar(1)=y(1);
    xHat=xBar;
    PBar=dynModel.initCov;
    PHat=PBar;
    l=1;
    %%% Kalman filter forward pass
    for k = 1:length(t_i_valid)
        %TU - Time update
        xBar = dynModel.Phi*xHat;
        PBar = dynModel.Phi*PHat*dynModel.Phi' + dynModel.Q;
        %Store
        x_bar_f(:,k)=xBar;
        P_bar_f(:,:,k)=PBar;

        measUpdateDone=0;
        %MU - Measurement Update only when we have a measurement
        while length(t)>=l && t_i_valid(k)>=t(l)  % Interpolated time has passed one 
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
                    outliers(l)=true;
                    disp(['Forward pass flagged measurement as outlier: t = ' num2str(t(l)) ' [min], y = ' num2str(y(l)) ' [mmol/L].'])
                end
            end
            if ~outliers(l)
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
        %limit to strictly positive for those states that are strictly
        %positive
        xHat(dynModel.strictlyPositiveStates & xHat<0)=0;
        
        % Store
        x_hat_f(:,k)=xHat;
        P_hat_f(:,:,k)=PHat;
    end % for k

    %%% Rauch-Tung-Striebel backward pass
    x_smoothed(:,k)=xHat;
    P_smoothed(:,:,k)=PHat;
    for k = length(t_i_valid)-1:-1:1
        C=(P_hat_f(:,:,k)*dynModel.Phi')/P_bar_f(:,:,k+1);
        x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
        P_smoothed(:,:,k)=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';
        %limit to strictly positive for those states that are strictly
        %positive
        x_smoothed(dynModel.strictlyPositiveStates & x_smoothed<0)=0;
    end
    
    %Generate output struct
    output.y_smoothed = nan(size(t_i));
    output.y_smoothed_sd = nan(size(t_i));
    output.y_smoothed(t_i_first:t_i_last) = x_smoothed(1,:);
    for k = 1:length(t_i_valid)
        output.y_smoothed_sd(t_i_first-1+k) = sqrt(P_smoothed(1,1,k));
    end
    
    if parsedArgs.outlierRemoval == 1
        %Run through all measurements and see if any are outside the error
        %smoothed band, if so they are outliers
        foundNewOutliers = false;
        y_s_mean  = closestValues(t,t_i,output.y_smoothed, startDateTime);
        y_s_sd  = closestValues(t,t_i,output.y_smoothed_sd, startDateTime);
        for i = 1:length(y)
            if ~outliers(i) && abs(y(i)-y_s_mean(i))/y_s_sd(i)>parsedArgs.outlierSDlimit
               outliers(i)=true;
               foundNewOutliers = true;
               disp(['Smoother flagged measurement ' num2str(i) ' as outlier: t = ' num2str(t(i)) ' [min], y = ' num2str(y(i)) ' [mmol/L].']) 
            end
        end
        if ~foundNewOutliers
            doneFindingOutliers = true;
        else
            disp(['Smoother needs a second pass due to outliers detected. Total # outliers in input data:' num2str(sum(outliers))]) 
            
        end
    else
        doneFindingOutliers = true;
    end
end %while not doneFindingOutliers
%Smoothing done

output.outliers = closestValues(t_in,t,outliers,startDateTime)==1;

output.y_filtered = nan(size(t_i));
output.y_filtered_sd = nan(size(t_i));
    
output.y_filtered(t_i_first:t_i_last) = x_hat_f(1,:);
for k = 1:length(t_i_valid)
    output.y_filtered_sd(t_i_first-1+k) = sqrt(P_hat_f(1,1,k));
end
output.t_i = t_i;
if isdatetime(startDateTime)
    output.t_i_relative = t_i;
    output.t_i = convertToAbsoluteTime(t_i, startDateTime);
end

%Add internal states
output.x_filtered = nan(size(x_hat_f,1),size(t_i,2));
output.x_smoothed = nan(size(x_smoothed,1),size(t_i,2));
output.x_filtered(:,t_i_first:t_i_last) = x_hat_f;
output.x_smoothed(:,t_i_first:t_i_last) = x_smoothed;

%Add user supplied wanted times
if length(parsedArgs.tout)==0
    parsedArgs.tout = t_in;
end
output.y_smoothed_at_tout = interpolatedValues(parsedArgs.tout,t_i,output.y_smoothed, startDateTime);
output.y_smoothed_sd_at_tout = interpolatedValues(parsedArgs.tout,t_i,output.y_smoothed_sd, startDateTime);

%Add dynModel
output.dynModel = dynModel;

%Plot if specified
if isdatetime(startDateTime)
    t_plot = startDateTime + minutes(t); % Convert to datetime for plotting
    t_i_plot = startDateTime + minutes(t_i); % Convert to datetime for plotting
else
    t_plot = t;
    t_i_plot = t_i;
end
  

if parsedArgs.plotResult>=1
    figure()
    
    %h1 = errorbar(t_plot,y,y_error,'r.','MarkerSize',10,'LineWidth',1); 
    h1 = plot(t_plot,y,'r.','MarkerSize',10); 
    
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
    figure();
    N=length(dynModel.stateNames);
    for sp=1:N
        subplot(N,1,sp)
        plot(t_i_plot,x_hat_f(sp,:));
        hold on
        plot(t_i_plot,x_smoothed(sp,:))
        ylabel(dynModel.stateNames(sp))
        legend('filtered','smoothed')
        if sp==1
            title(['Internal states model ' num2str(dynModel.id)])
            plot(t_plot,y,'r.','MarkerSize',15)
            legend('filtered','smoothed','measurements')
            if sum(outliers)>0
                plot(t_plot(outliers),y(outliers),'kx','MarkerSize',15)
                legend('filtered','smoothed','measurements','outliers')
            end
        else
            if parsedArgs.dynamicModel==2
                plot([t_plot(1) t_plot(end)],[0 0],'k-')
                legend('filtered','smoothed','zero line')
            end
        end
        hold off
            
    end
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
end%function 

%Helper method to parse the arguments
function parsedArgs = parseInputVarArgs(varargs)
    %Set defaults
    parsedArgs.outlierRemoval = 0;
    parsedArgs.outlierSDlimit = 2;
    parsedArgs.plotResult = 0;
    parsedArgs.plotInternalStates = 0;
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









