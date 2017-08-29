function unit = autoDetectGlucoseUnit( measurements )
%AUTODETECTGLUCOSEUNIT Helper function to autodetect glucose unit
    unit = 'mmol_L';
    if mean(measurements)>50
        disp('SmoothGlucoseData autodetected mg/dL as unit')
        unit = 'mg_dL';
    end
end

