% function [signal] = KD_fitting(xfit, con, na, B1)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% a = xfit(2);
% b = xfit(3);
% c = xfit(1);
% 
% signal = (B1^2 / (c^2 + B1^2))*a*na*con./(con+b);
% 
% 
% end

function [signal] = KD_fitting(xfit, B1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a = xfit(1);
b = xfit(2);

signal = (B1.^2 ./ (a + B1.^2))*b;

end