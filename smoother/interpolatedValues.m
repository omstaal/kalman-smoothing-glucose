function yout = interpolatedValues(tout, t, y, startdatetime)
    %INTERPOLATEDVALUES Helper function to find values at interpolated output times
    if isdatetime(tout)
        tout.TimeZone='';
        tout = convertToRelativeTime(tout, startdatetime);
    end
    if sum(diff(tout)<0)>0
        error('tout was not monotonously increasing')
    end
    yout = zeros(size(tout));
    lastj = 1;
    for i = 1:length(tout)
        if tout(i)>t(end) || tout(i)<t(1)
            %Do not extrapolate
            yout(i) = NaN;
        else
            for j=lastj:length(t)
                if t(j)==tout(i)
                    yout(i) = y(j);
                    lastj = j;
                elseif t(j)>tout(i)
                    y1 = y(j-1);
                    y2 = y(j);
                    t1 = t(j-1);
                    t2 = t(j);
                    
                    yout(i) = y1 + (tout(i)-t1)*(y2-y1)/(t2-t1);
                    %yout(i) = y1*(t2-tout(i)) + y2*(tout(i) - t1); %Should
                    %be the same
                    
                    lastj = j;
                    break;
                end
            end
        end
    end
end

