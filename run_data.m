clear
clc
close all

run_file=mfilename('fullpath');
if isempty(run_file)
    run_file=which('run_data.m');
end
if isempty(run_file)
    error('Cannot determine the location of run_data.m.');
end

project_dir=fileparts(run_file);
data_dir=fullfile(project_dir,'data');
src_dir=fullfile(project_dir,'src');

assert(isfolder(src_dir),'Cannot find source folder: %s',src_dir);
addpath(src_dir,'-begin');
rehash;

required_functions={'wcalsin','fit_sine','fft_metrics','alias_frequency','plot_fft','findbin','findfreq','local_rms'};
for i=1:numel(required_functions)
    assert(exist(required_functions{i},'file')==2,'Cannot find %s.m',required_functions{i});
end

plot_style=struct(...
    'blue',[0 0.4470 0.7410],...
    'grid',[0.68 0.68 0.68],...
    'font','Arial',...
    'tick',10,...
    'label',11,...
    'text',10,...
    'width',1.1,...
    'trace',0.425);

lf_file=load(fullfile(data_dir,'lf_measured_raw.mat'),'DATA');
hf_file=load(fullfile(data_dir,'hf_measured_raw.mat'),'DATA');
LF=lf_file.DATA;
HF=hf_file.DATA;

lf_bits=double(LF.bits);
hf_bits=double(HF.bits);
channel_count=4;
sample_count=size(lf_bits,1);
bit_count=size(lf_bits,2);

assert(sample_count==8192 && size(hf_bits,1)==8192,'This release expects 8192 samples per capture.');
assert(bit_count==20 && size(hf_bits,2)==20,'This release expects 20 internal bit columns.');
assert(LF.config.Fs_Hz==4e9 && HF.config.Fs_Hz==4e9,'This release expects Fs = 4 GS/s.');
assert(all(lf_bits(:)==0 | lf_bits(:)==1),'LF DATA.bits must contain only 0/1.');
assert(all(hf_bits(:)==0 | hf_bits(:)==1),'HF DATA.bits must contain only 0/1.');

lf_bin_error=abs(LF.config.normalized_frequency*sample_count-LF.config.tone_bin);
hf_bin_error=abs(HF.config.normalized_frequency*sample_count-HF.config.tone_bin);
if lf_bin_error>1e-8
    warning('run_data:LFNoncoherent','LF tone is not coherent with the supplied tone_bin.');
end
if hf_bin_error>1e-8
    warning('run_data:HFNoncoherent','HF tone is not coherent with the supplied tone_bin.');
end

lf_channel_frequency=alias_frequency(LF.config.normalized_frequency*channel_count);

lf_channel_bits=cell(1,channel_count);
hf_channel_bits=cell(1,channel_count);
weight=zeros(channel_count,bit_count);
frequency_used_lf=zeros(1,channel_count);

for channel=1:channel_count
    lf_index=channel:channel_count:sample_count;
    hf_index=channel:channel_count:sample_count;
    lf_channel_bits{channel}=lf_bits(lf_index,:);
    hf_channel_bits{channel}=hf_bits(hf_index,:);
    [weight(channel,:),~,~,~,~,frequency_used_lf(channel)]=wcalsin(...
        lf_channel_bits{channel},lf_channel_frequency,0.5,1e-12,100,1,0,0,1,1,[]);
end

lf_weighted=cell(1,channel_count);
hf_weighted=cell(1,channel_count);
for channel=1:channel_count
    lf_weighted{channel}=lf_channel_bits{channel}*weight(channel,:).';
    hf_weighted{channel}=hf_channel_bits{channel}*weight(channel,:).';
end

channel_fixed_offset=zeros(1,channel_count);
channel_fixed_amplitude=zeros(1,channel_count);
channel_gain=zeros(1,channel_count);
for channel=1:channel_count
    [channel_fixed_offset(channel),channel_fixed_amplitude(channel)]=fit_sine(...
        lf_weighted{channel},lf_channel_frequency);
    if ~(isfinite(channel_fixed_amplitude(channel))&&channel_fixed_amplitude(channel)>0)
        error('Invalid LF fitted amplitude for channel %d.',channel);
    end
    channel_gain(channel)=1/channel_fixed_amplitude(channel);
end

lf_channel_postcal=cell(1,channel_count);
hf_channel_postcal=cell(1,channel_count);
for channel=1:channel_count
    lf_channel_postcal{channel}=(lf_weighted{channel}-channel_fixed_offset(channel))*channel_gain(channel);
    hf_channel_postcal{channel}=(hf_weighted{channel}-channel_fixed_offset(channel))*channel_gain(channel);
end

lf_postcal=zeros(sample_count,1);
hf_postcal=zeros(sample_count,1);
for channel=1:channel_count
    lf_postcal(channel:channel_count:end)=lf_channel_postcal{channel};
    hf_postcal(channel:channel_count:end)=hf_channel_postcal{channel};
end

lf_tone_bin=resolve_tone_bin(LF.config,lf_postcal);
hf_tone_bin=resolve_tone_bin(HF.config,hf_postcal);

[lf_frequency_MHz,lf_fft_db,~,lf_sndr_db,lf_sfdr_db]=fft_metrics(...
    lf_postcal,4e9,lf_tone_bin);
[hf_frequency_MHz,hf_fft_db,~,hf_sndr_db,hf_sfdr_db]=fft_metrics(...
    hf_postcal,4e9,hf_tone_bin);

result_figure=figure(...
    'Color','w',...
    'Visible','on',...
    'Units','pixels',...
    'Position',[100 80 1060 430],...
    'Name','ADC Calibration',...
    'NumberTitle','off');

ax_lf=axes(result_figure,'Position',[0.075 0.145 0.405 0.805]);
ax_hf=axes(result_figure,'Position',[0.565 0.145 0.405 0.805]);
lf_plot_config=LF.config;
hf_plot_config=HF.config;
lf_plot_config.tone_bin=lf_tone_bin;
hf_plot_config.tone_bin=hf_tone_bin;
plot_fft(ax_lf,lf_frequency_MHz,lf_fft_db,lf_plot_config,lf_sndr_db,lf_sfdr_db,plot_style);
plot_fft(ax_hf,hf_frequency_MHz,hf_fft_db,hf_plot_config,hf_sndr_db,hf_sfdr_db,plot_style);

drawnow;

fprintf('\n==============================================\n');
fprintf('LF CAL : SNDR %.4f dB, SFDR %.4f dB\n',lf_sndr_db,lf_sfdr_db);
fprintf('HF VAL : SNDR %.4f dB, SFDR %.4f dB\n',hf_sndr_db,hf_sfdr_db);
fprintf('==============================================\n');

function tone_bin=resolve_tone_bin(config,x)
if isfield(config,'tone_bin') && isscalar(config.tone_bin) && isfinite(config.tone_bin)
    tone_bin=round(double(config.tone_bin));
elseif isfield(config,'Fin_Hz') && isscalar(config.Fin_Hz) && isfinite(config.Fin_Hz)
    tone_bin=findbin(4e9,double(config.Fin_Hz),8192);
else
    fin=findfreq(x,4e9);
    tone_bin=findbin(4e9,fin,8192);
end
end
