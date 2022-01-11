% function [signal] = KD_fitting_w(xfit, B1)
% %UNTITLED Summary of this function goes here
% %   Detailed explanation goes here
% a = xfit(1);
% b = xfit(2);
% c = xfit(3);
% 
% signal = (0.5*(a*B1.^2) ./ (1+a*B1.^2) + 0.5*B1.^2 ./ (b^2+B1.^2)) * c;
% 
% end

% function [signal] = KD_fitting_w(xfit, B1, kbind, Rou)
% Ffree = xfit(1);
% T1aT2a = xfit(2);
% Fbound = 1-Ffree;
% C = xfit(3);
% 
% signal = (Ffree*(T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2) + Fbound*B1.^2 ./ (kbind*Rou+B1.^2)) * C;
% 
% end

% function [signal] = KD_fitting_w(xfit, B1, Ffree, Rou,kbind)
% % kbind = xfit(1);
% T1aT2a = xfit(1);
% Fbound = 1-Ffree;
% C = xfit(2);
% 
% signal = (Ffree*(T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2) + Fbound*B1.^2 ./ (kbind*Rou+B1.^2)) * C;
% 
% end


% function [signal] = KD_fitting_w(xfit, B1, kbind)
% 
% T1aT2a = xfit(1);
% C = xfit(2);
% 
% signal = ((T1aT2a*B1.^2) ./ (1+T1aT2a*B1.^2)) * C;
% 
% end

function [signal] = KD_fitting_w(xfit, B1)
% kbind = xfit(1);
pq = xfit(1);
C = xfit(2);

signal = B1.^2 ./ (pq+B1.^2) * C;

end
