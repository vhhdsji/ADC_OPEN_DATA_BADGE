function [frequency_MHz,fft_db,fft_magnitude,sndr_db,sfdr_db]=fft_metrics(x,Fs,tone_bin)
x=double(x(:));
N=numel(x);
if N<4 || mod(N,2)~=0
    error('fft_metrics:InvalidLength','Input length must be an even integer of at least four samples.');
end
if ~(isscalar(Fs)&&isfinite(Fs)&&Fs>0)
    error('fft_metrics:InvalidFs','Fs must be a positive finite scalar.');
end
if ~(isscalar(tone_bin)&&isfinite(tone_bin)&&tone_bin==round(tone_bin)&&tone_bin>0&&tone_bin<N/2)
    error('fft_metrics:InvalidToneBin','tone_bin must be an integer in the range 1 to N/2-1.');
end
x=x-mean(x);
Y=fft(x);
fft_magnitude=abs(Y(1:N/2+1))/N;
fft_magnitude(2:end-1)=2*fft_magnitude(2:end-1);
power_spectrum=abs(Y(1:N/2+1)/N).^2;
power_spectrum(2:end-1)=2*power_spectrum(2:end-1);
fft_db=20*log10(max(fft_magnitude,10^(-140/20)));
frequency_MHz=(0:N/2).'/N*Fs/1e6;
k=tone_bin+1;
mask=true(size(power_spectrum));
mask([1 k])=false;
noise_distortion_power=sum(power_spectrum(mask));
spur_power=max(power_spectrum(mask));
sndr_db=10*log10(power_spectrum(k)/noise_distortion_power);
sfdr_db=10*log10(power_spectrum(k)/spur_power);
end
