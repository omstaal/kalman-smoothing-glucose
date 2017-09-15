function y_error = setIsoError(y)
    %SETISOERROR Helper method that guesses a error based on the measured value, based on ISO 15197 
    %Assumes y is given in mmol/L
    y_error = zeros(size(y));              % ISO 15197, set it based on the measured values
    for i = 1:length(y) % Make error bars (assumes fingerprick measurement
                        % errors according to ISO15197)
        if y(i)>5.55
            y_error(i) = 0.15*y(i);
        else
            y_error(i) = 0.83;
        end
    end
end
