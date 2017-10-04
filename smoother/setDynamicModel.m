%Helper method to set dynamic model
function dynModel=setDynamicModel(dynModelNo)
    dynModel.delta_t=1/6; %1 min/6 = 10 sec
    if(dynModelNo==1) %simple 2.order system where the rate of change of glucose dies out.
        dynModel.id=1;
        a=-0.05;
        qm1 = 0.005*dynModel.delta_t;
        dynModel.F =[0 1;0 a];               % System matrix (continuous)
        dynModel.Q=[0 0;0 qm1];     % Process noise covariance matrix.
        dynModel.H=[1 0];                    % Measurement matrix.
        dynModel.initCov = diag([0.25 1]);   % Initial covariance
        %%% Discretization
        dynModel.Phi=expm(dynModel.F*dynModel.delta_t);        % Discrete state transition matrix
        dynModel.stateNames = {'Gp','dGp'};
        dynModel.strictlyPositiveStates = [true;false];
    elseif (dynModelNo==2) %Lumped Tinsulin/meal state that is central/remote
        dynModel.id = 2;
        Td = 10.000; % Time constant [min] describing flow between central and remote compartments
        qm2 = 0.02*dynModel.delta_t;
        dynModel.F =[0 0 1;0 -1/Td 0;0 1/Td -1/Td]; % System matrix (continuous)
        dynModel.Q=[0 0 0;0 qm2 0;0 0 0]; % Process noise covariance matrix.
        dynModel.H=[1 0 0];                       % Measurement matrix.
        dynModel.initCov = diag([10 1 1]);         % Initial covariance
        dynModel.Phi=expm(dynModel.F*dynModel.delta_t);                    % Discrete state transition matrix
        dynModel.stateNames = {'Gp','C','R'};
        dynModel.strictlyPositiveStates = [true;false;false];
    elseif (dynModelNo==3) %Insulin and meal, central
        dynModel.id = 3;
        Ti = 20.000; % Time constant describing insulin flow between compartments [min]
        Tm = 10.000; % Time constant describing meal flow between compartments [min]
        qm3i = 0.01*dynModel.delta_t;
        qm3m = 0.01*dynModel.delta_t;
        dynModel.F =[0 -1 1;0 -1/Ti 0;0 0 -1/Tm]; % System matrix (continuous)
        dynModel.Q=[0 0 0;0 qm3i 0;0 0 qm3m]; % Process noise covariance matrix.
        dynModel.H=[1 0 0];                        % Measurement matrix.
        dynModel.initCov = diag([10 1 1]);         % Initial covariance
        dynModel.Phi=expm(dynModel.F*dynModel.delta_t);              % Discrete state transition matrix
        dynModel.stateNames = {'Gp','I','M'};
        dynModel.strictlyPositiveStates = true(3,1);
    else
        error(['Unsupported model:' num2str(model)])
    end
    
    if ~rank(obsv(dynModel.F,dynModel.H))==size(dynModel.F,1)
        error('System is not observable')
    end
end

