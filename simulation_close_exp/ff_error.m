function [f]  = ff_error(x,Con,na,rawData)
a=x(1);
b=x(2);

signal = a*na * Con./(Con+b);
% f = signal - rawData;
f = sqrt((signal - rawData).^2);

end