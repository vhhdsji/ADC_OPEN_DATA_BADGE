function y=local_rms(x)
x=double(x);
y=sqrt(mean(abs(x(:)).^2));
end
