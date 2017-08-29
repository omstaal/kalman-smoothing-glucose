%Helper method to convert glucose unit from mmol/L to mg/dL
function y_mgdL = convertTo_mg_dL(y_mmolL)
    y_mgdL = y_mmolL*18.018;
end

