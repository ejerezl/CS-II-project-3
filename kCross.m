function value = kCross(x, y, crossLength, crossWidth, crossValue, restValue)
value = restValue;
if ((abs(0.5-x)<crossLength/2 && abs(0.5-y)<crossWidth/2)...
        || (abs(0.5-x)<crossWidth/2 && abs(0.5-y)<crossLength/2))
    value = crossValue;
end
end