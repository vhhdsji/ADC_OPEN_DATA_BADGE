function [offset,amplitude,phase]=fit_sine(x,f)
x=double(x(:));
n=(0:numel(x)-1).';
A=[cos(2*pi*f*n) sin(2*pi*f*n) ones(size(n))];
c=A\x;
amplitude=hypot(c(1),c(2));
phase=atan2(-c(2),c(1));
offset=c(3);
end
