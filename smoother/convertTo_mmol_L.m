%Helper method to convert glucose unit from mg/dL to mmol/dL
function y_mmolL = convertTo_mmol_L(y_mgdL)
    y_mmolL = y_mgdL/18.018;
end
