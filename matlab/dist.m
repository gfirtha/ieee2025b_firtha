function [d] = dist(x)
%DIST Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    x
end

arguments (Output)
    d
end

    d = sqrt( sum(x.^2,2) );
end