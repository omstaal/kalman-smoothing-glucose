function yout = closestValues(tout, t, y, startdatetime)
    %CLOSESTVALUES Helper function to find closest to wanted output times
    if isdatetime(tout)
        tout = convertToRelativeTime(tout, startdatetime);
    end
    if sum(diff(tout)<0)>0
        error('tout was not monotonously increasing')
    end
    yout = zeros(size(tout));
    lastj = 1;
    for i = 1:length(tout)
        if tout(i)>t(end) || tout(i)<t(1)
            yout(i) = NaN;
        else
            for j=lastj:length(t)
                if t(j)>=tout(i)
                    yout(i) = y(j);
                    lastj = j;
                    break;
                end
            end
        end
    end
end

