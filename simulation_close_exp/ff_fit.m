function [signal]  = ff_fit(x,Con,na)
a=x(1);
b=x(2);

signal = a*na .* Con./(Con+b);

end