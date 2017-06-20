function [y_smoothed,y_smoothed_sd, varargout]=SmoothGlucoseData(t,y,t_i,outlierRemoval,plotResult)
% SMOOTHGLUCOSEDATA Creates a smoothed glucose curve with variance estimates
% Usage:
% [y_smoothed,y_smoothed_sd]=SmoothGlucoseData(t,y,t_i,plot)
%   Generates a smoothed estimate from the input data t (time in minutes) 
%   and y (glucose values) with at the times given by t_i (minutes, should 
%   have a uniform sampling intervals less than those in t). 
%   mean y_smoothed and standard deviation y_smoothed_sd
%
%   The outlierRemoval parameter controls outlier removal:
%   Measurements are considered to be outliers if the innovation is outside
%   X std devs of the innovation variance in the forward pass Kalmanfilter
%   If outlierRemoval==1, a limit of 5 std dev is used (conservative) 
%   If outlierRemoval==2, a limit of 3 std devs (medium)
%   If outlierRemoval==3, a limit of 2 std devs (aggressive)
%   If outlierRemoval is any other value, outlier removal is not done. 
%   
% If plotResult==1, a plot will be produced showing the result of the
%   smoothing
% If the four-output variant is used:
% [y_smoothed,y_smoothed_sd,y_filtered, y_filtered_sd]=SmoothGlucoseData(t,y,t_i,plot)
% ,the data from the forward pass is
% also returned (and plotted if specified)

%%TODO time handling should be better. Should handle datetime arrays

%%% Tunable parameters
% The system used by the Kalman smoother is very simple and only describes a state that 
% has a rate of change that decays 
% $$ \dot{x_1} = x_2 $$
% $$ \dot{x_2} = -ax_2 $$
outlierStds=NaN;
if outlierRemoval==1
    outlierStds=5;
elseif outlierRemoval==2
    outlierStds=3;
elseif outlierRemoval==3
    outlierStds=2;    
end
delta_t = (max(t_i) - min(t_i))/(length(t_i)-1);%Time step of interpolated signal (minutes)
if delta_t>0.5
    error('delta_t larger than half a minute, not recommended')
end
a=-0.05;
F =[0 1;0 a];%System matrix (continuous)- simple 2.order system where the rate of change of glucose dies out
Q=[0 0;0 0.05*delta_t];%Process noise covariance matrix
H=[1 0];%Measurement matrix
%According to ISO15197
Rlow = 0.11;%Measurement variance - (mmol/L)^2
isoLimit = 5.6;%mmol/L
RhighFact = 0.0064;

%%% Discretization
Phi=expm(F*delta_t);%Discrete state transition matrix

%%% Storage
x_hat_f = zeros(2,length(t_i));%A priori state vector storage, forward pass
x_bar_f = zeros(2,length(t_i));%A posteriori state vector storage, forward pass
P_hat_f = zeros(2,2,length(t_i));%A priori covariance matrix storage, forward pass
P_bar_f = zeros(2,2,length(t_i));%A posteriori covariance matrix storage, forward pass
x_smoothed = zeros(2,length(t_i));%State vector storage, backward pass
P_smoothed = zeros(2,2,length(t_i));%Covariance matrix storage, backward pass


%%% Initialization
xBar=[y(1);0];
xHat=xBar;
PBar =[0.25 0;0 1];
PHat=PBar;

l=1;
%Remove any nan entries
nonNan = ~isnan(y) & ~isnan(t);
y = y(nonNan);
t = t(nonNan);

%%% Kalman filter forward pass
for k = 1:length(t_i)
    %TU - Time update
    xBar = Phi*xHat;
    PBar = Phi*PHat*Phi' + Q;
    %Store
    x_bar_f(:,k)=xBar;
    P_bar_f(:,:,k)=PBar;
    
    measUpdateDone=0;
    %MU - Measurement Update only when we have a measurement
    while length(t)>=l && t_i(k)>=t(l)%Interpolated time has passed one of the measurement times
        if measUpdateDone==1
            %More than one measurement at the current time
            xBar = xHat;
            PBar = xBar;
        end
        dz = y(l)-H*xBar;
        if(y(l)>isoLimit)%Set measurement variance for this measurement
            R=RhighFact*y(l)^2;
        else
            R=Rlow;
        end
        Pz = (H*PBar*H'+R);
        if isnan(outlierStds)
            isOutlier = 0;
        else
            %Check the innovation
            if(abs(dz)>outlierStds*sqrt(Pz))
                isOutlier=1;
                disp(['Flagged measurement as outlier: (t,y)=' num2str([t(l) y(l)]) ])
            else
                isOutlier=0;
            end
        end
        if isOutlier==0
            K=PBar*H'/Pz;
            xHat = xBar + K*dz;
            PHat = (eye(size(PBar))-K*H)*PBar;
            measUpdateDone=1;
        end
        l=l+1;
    end
    if measUpdateDone==0%No measurement was available at this time 
        xHat=xBar;
        PHat=PBar;
    end
    %Store
    x_hat_f(:,k)=xHat;
    P_hat_f(:,:,k)=PHat;
end %for k

%%% Rauch-Tung-Striebel backward pass
x_smoothed(:,k)=xHat;
P_smoothed(:,:,k)=PHat;
for k = length(t_i)-1:-1:1
    C=P_hat_f(:,:,k)*Phi'*inv(P_bar_f(:,:,k+1));
    x_smoothed(:,k)=x_hat_f(:,k)+C*(x_smoothed(:,k+1)-x_bar_f(:,k+1));
    P_smoothed(:,:,k)=P_hat_f(:,:,k)+C*(P_smoothed(:,:,k+1)-P_bar_f(:,:,k+1))*C';  
end
y_smoothed = x_smoothed(1,:);
y_smoothed_sd = zeros(size(y_smoothed));
for k = 1:length(t_i)
    y_smoothed_sd(k) = sqrt(P_smoothed(1,1,k));
end

if nargout==4
    disp('Adding filtered (intermediate result before smoothing step)')
    y_filtered = x_hat_f(1,:);
    y_filtered_sd = zeros(size(y_filtered));
    for k = 1:length(t_i)
        y_filtered_sd(k) = sqrt(P_hat_f(1,1,k));
    end
    varargout{1} = y_filtered;
    varargout{2} = y_filtered_sd;
end

if plotResult==1
    figure()
    y_error = zeros(size(y));
    for i = 1:length(y)%Make error bars (assumes fingerprick measurement errors according to ISO15197)
        if y(i)>5.6
            y_error(i) = 0.2*y(i);
        else
            y_error(i) = 0.83;
        end
    end
    h1 = errorbar(t,y,y_error,'r.','MarkerSize',10,'LineWidth',1);%Using 2 sigma
    hold on
    h2 = plot(t_i,y_smoothed,'b','LineWidth',2);
    h3 = plot(t_i,y_smoothed+2*y_smoothed_sd,'b--');
    plot(t_i,y_smoothed-2*y_smoothed_sd,'b--')

    xlabel('Time [min]')
    ylabel('Glucose [mmol/L]')
    if nargout == 4
        h4= plot(t_i,y_filtered,'g','LineWidth',2);
        h5 = plot(t_i,y_filtered+2*y_filtered_sd,'g--');
        plot(t_i,y_filtered-2*y_filtered_sd,'g--')
        legend([h1 h2 h3 h4 h5],{'Fingerprick measurements w/error','Smoothed estimate','\pm2 SD of smoothed estimate', 'Filtered estimate', '\pm2 SD of filtered estimate' })
    else
        legend([h1 h2 h3],{'Fingerprick measurements w/error','Smoothed estimate','\pm2 SD of smoothed estimate' })
    end
    hold off
end

end%function

%%% R value rationale: 
% Blood glucose meters are supposed to lie within 20% of the real
% value 
% The ISO15197:2013 standard specifies error less than +-0.83 for BG <= 5.6,
% and less than 20% of the true value for BG>5.6 mmol/L.
% Low level: Assuming normal distribution, this means that 2 standard deviations (95% CI)
% at a normal level of 5 mmol/L shall be less than 0.83 mmol/L => 2SD = 0.415mmol/L 
% => SD^2=VAR=0.172 (mmol/L)^2
% High level: 2SD=0.2*y, SD=0.1y, VAR = 0.01y^2



