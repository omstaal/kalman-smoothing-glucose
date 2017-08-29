function yout = closestValues(tout, t, y, startdatetime)
    %CLOSESTVALUES Helper function to find closest to wanted output times
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

