function [d, rate] = ric_rate(z,bet)
b = [bet(1)-bet(2),bet(1)+bet(2)]-z;
if(imag(z)<0)
    b = b - 1j*bet(3);
end

ba = abs(b);
d  = sum(ba./b) / sum(ba);
rate = abs(b(2)-b(1)) / sum(ba);
end