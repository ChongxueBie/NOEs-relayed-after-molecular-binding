% function [f] = KD_fitting_lsq(x,con, rawData, na, B1)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% a = x(2);
% b = x(3);
% c = x(1);
% 
% signal = (B1^2 / (c^2 + B1^2))*a*na*con./(con+b);
% f = signal - rawData;
% 
% 
% end

function [f] = KD_fitting_lsq(x,B1, rawData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = x(1);
b = x(2);

signal = (B1.^2 ./ (a + B1.^2))*b;
f = signal - rawData;


end
