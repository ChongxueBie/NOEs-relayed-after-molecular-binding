% 
% function [f] = KD_fitting_w_lsq(x,B1, rawData)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% a = x(1);
% b = x(2);
% c = x(3);
% 
% signal = (0.5*(a*B1.^2) ./ (1+a*B1.^2) + 0.5*B1.^2 ./ (b^2+B1.^2)) * c;
% 
% f = signal - rawData;
% 
% 
% end

% 
% function [f] = KD_fitting_w_lsq(x,B1, rawData, kbind, Rou)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% Ffree = x(1);
% T1aT2a = x(2);
% Fbound = 1-Ffree;
% C = x(3);
% 
% signal = (Ffree*((T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2)) + Fbound*(B1.^2 ./ (kbind*Rou+B1.^2))) * C;
% 
% f = sqrt((signal - rawData).^2);
% 
% 
% end

% function [f] = KD_fitting_w_lsq(x,B1, rawData, Ffree, Rou,kbind)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% % kbind = x(1);
% T1aT2a = x(1);
% Fbound = 1-Ffree;
% C = x(2);
% 
% signal = (Ffree*((T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2)) + Fbound*(B1.^2 ./ (kbind*Rou+B1.^2))) * C;
% 
% f = sqrt((signal - rawData).^2);
% 
% 
% end

% function [f] = KD_fitting_w_lsq(x,B1, rawData, kbind)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% T1aT2a = x(1);
% C = x(2);
% 
% signal = ((T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2)) * C;
% 
% f = sqrt((signal - rawData).^2);
% 
% 
% end

function [f] = KD_fitting_w_lsq(x, B1, rawData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% kbind = x(1);
pq = x(1);
C = x(2);

signal = B1.^2 ./ (pq+B1.^2) * C;

f = sqrt((signal - rawData).^2);


end