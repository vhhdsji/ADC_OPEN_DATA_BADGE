function f=alias_frequency(f)
f=mod(f,1);
if f>0.5
    f=1-f;
end
end
