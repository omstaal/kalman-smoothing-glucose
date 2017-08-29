%Helper method to convert datetimes to relative time
function t = convertToRelativeTime(datetimes,startdatetime)
    t=zeros(length(datetimes),1);
    for i=1:length(datetimes)
        t(i) = minutes(datetimes(i)-startdatetime);
    end
end
