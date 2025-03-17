function y = tableInterpolation(x, xList, yList)
    y = yList(1) + (x - xList(1)) * (yList(2) - yList(1)) / ...
                                    (xList(2) - xList(1));
end