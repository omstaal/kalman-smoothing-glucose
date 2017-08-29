%Helper method to convert relative time to absolute time
function datetimes = convertToAbsoluteTime(reltimes,startdatetime)
    datetimes=datetime(zeros(length(reltimes),3));
    for i=1:length(datetimes)
        datetimes(i) = startdatetime + minutes(reltimes(i));
    end
end
