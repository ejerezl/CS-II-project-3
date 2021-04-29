function value = kSquare(x, y, squareLength, innerLength, squareValue, restValue)
value = restValue;
if ((abs(0.5-x)<squareLength/2 && abs(0.5-y)>innerLength/2 && abs(0.5-y)<squareLength/2)...
        || (abs(0.5-x)<squareLength/2 && abs(0.5-y)<squareLength/2 && abs(0.5-x)>innerLength/2))
    value = squareValue;
end
end