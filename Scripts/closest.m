function [cl,closestIndex] = closest(arr,val) 
%{
 Find corresponding index:
    https://www.mathworks.com/matlabcentral/answers/
    152301-find-closest-value-in-array#answer_1559017
%}
    [~,closestIndex] = min(arr-val.', [], ComparisonMethod = "abs");
    cl = arr(closestIndex);
end