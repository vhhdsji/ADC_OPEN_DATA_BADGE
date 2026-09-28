function plot_figure5_report
% Figure 5 report: display the four measurement plots and power breakdown.
% Measured SNDR/SFDR versus input frequency, multi-chip SNDR under supply
% voltage and temperature variations, DR at a 170MHz input, and power breakdown.
% All curve values, ticks, reference values and annotation positions are in MAT files.
% Numeric constants below control drawing styles, canvas size and unit conversion only.
% Fixed ratios are written as decimal constants. Variable divisions calculate data ratios.

%% 0. Load one MAT file per subplot
root=fileparts(mfilename('fullpath'));
F=load(fullfile(root,'data','frequency.mat'));
V=load(fullfile(root,'data','supply_voltage.mat'));
T=load(fullfile(root,'data','temperature.mat'));
D=load(fullfile(root,'data','dynamic_range_170MHz.mat'));
P=load(fullfile(root,'data','power_breakdown.mat'));
out=fullfile(root,'plots');
if ~exist(out,'dir'), mkdir(out); end
f=figure('Color','w','Units','inches','Position',[1 1 4.52 3.34],'Renderer','painters');
set(f,'DefaultTextFontName','Arial','DefaultTextFontSize',6.5,'DefaultTextFontWeight','bold');
set(f,'DefaultAxesLabelFontSizeMultiplier',1.12);
set(f,'PaperUnits','inches','PaperPosition',[0 0 4.52 3.34],'PaperSize',[4.52 3.34]);

%% 1. Measured SNDR/SFDR versus input frequency
ax=panelaxes(f,F.AxesRect_px);
drawpair(ax,F.Frequency_MHz,F.SNDR_dB,F.SFDR_dB);
configure(ax,F);
xline(ax,F.ReferenceFrequency_MHz,':k','HandleVisibility','off');
text(ax,F.ReferenceLabelPosition(1),F.ReferenceLabelPosition(2),F.ReferenceLabel,'FontSize',5.8);
plot(ax,F.ArrowX,F.ArrowY,'k-','LineWidth',0.7);
patch(ax,F.ArrowHeadX,F.ArrowHeadY,'k','EdgeColor','none');
text(ax,F.ValueLabelPosition(1),F.ValueLabelPosition(2),sprintf('%.1f',F.SNDR_dB(end)));
pairlegend(f,F.LegendRect_px);

%% 2. Multi-chip SNDR under supply voltage variations
comparison(f,V,V.SupplyVoltage_mV,V.XTick_mV);

%% 3. Multi-chip SNDR under temperature variations
comparison(f,T,T.Temperature_C,T.XTick_C);

%% 4. DR at a 170MHz input
ax=panelaxes(f,D.AxesRect_px);
drawpair(ax,D.SignalPower_dBFS,D.SNDR_dB_recovered,D.SFDR_dB_recovered);
configure(ax,D);
plot(ax,[-D.DR_dB 0],[0 0],'k-');
text(ax,D.DRLabelPosition(1),D.DRLabelPosition(2),sprintf('DR = %.1f dB',D.DR_dB));
pairlegend(f,D.LegendRect_px);

%% 4a. Enlarged inset near full-scale input
ax=panelaxes(f,D.InsetRect_px);
drawpair(ax,D.Inset_SignalPower_dBFS,D.Inset_SNDR_dB,D.Inset_SFDR_dB);
set(ax,'XLim',D.InsetXLim,'YLim',D.InsetYLim,'XTick',D.InsetXTick,'YTick',D.InsetYTick,'FontSize',4.5);

%% 5. Power breakdown
ax=axes(f,'Position',P.AxesPosition_percent*0.01);
hold(ax,'on');
colors=[0.3568627451 0.6078431373 0.8352941176; ...
    0.9294117647 0.4901960784 0.1921568627; ...
    0.6470588235 0.6470588235 0.6470588235; ...
    1.0 0.7529411765 0.0; ...
    0.4392156863 0.6784313725 0.2784313725];
angle=P.StartAngle_deg;
values=P.Power_mW_original;
for k=1:numel(values)
    span=360*values(k)/sum(values);
    theta=linspace(angle,angle-span,180)*0.0174532925199433;
    patch(ax,[cos(theta) .55*cos(fliplr(theta))],[sin(theta) .55*sin(fliplr(theta))],colors(k,:),'EdgeColor','none');
    mid=(angle-span*0.5)*0.0174532925199433;
    text(ax,.79*cos(mid),.79*sin(mid),sprintf('%.0f%%',P.Percent_display(k)),...
        'Color','w','HorizontalAlignment','center');
    angle=angle-span;
end
text(ax,0,0,sprintf('Total =\n%.1fmW',P.TotalPower_mW_display),'HorizontalAlignment','center');
positions=P.LabelPosition_percent*0.01;
for k=1:numel(values)
    block=P.BlockNames{k};
    if strcmp(block,'Reference Buffer'), block=sprintf('Reference\nBuffer'); end
    label=sprintf('%s\n%.1fmW',block,P.Power_mW_display(k));
    alignment=strtrim(P.LabelAlignment(k,:));
    text(ax,positions(k,1),positions(k,2),label,'HorizontalAlignment',alignment);
end
axis(ax,'equal'); axis(ax,'off');
xlim(ax,P.XLim_percent*0.01); ylim(ax,P.YLim_percent*0.01);

%% 6. Export editable figure and publication formats
savefig(f,fullfile(out,'figure5_report.fig'));
print(f,fullfile(out,'figure5_report.png'),'-dpng','-r400');
print(f,fullfile(out,'figure5_report.svg'),'-dsvg','-painters');
print(f,fullfile(out,'figure5_report.pdf'),'-dpdf','-painters');
end

%% Shared drawing styles
function ax=panelaxes(f,rect)
ax=axes(f,'Position',pixelposition(rect));
hold(ax,'on');
set(ax,'FontName','Arial','FontSize',5.8,'FontWeight','bold','TickDir','in','LineWidth',0.5,'GridAlpha',0.15);
grid(ax,'on'); box(ax,'on');
end

function position=pixelposition(rect)
position=[rect(1)*0.0022123893805310, ...
    1.0-(rect(2)+rect(4))*0.0029940119760479, ...
    rect(3)*0.0022123893805310,rect(4)*0.0029940119760479];
end

function configure(ax,S)
set(ax,'XLim',S.XLim,'YLim',S.YLim,'XTick',S.XTick,'YTick',S.YTick);
xlabel(ax,S.XLabel); ylabel(ax,S.YLabel);
end

function drawpair(ax,x,sndr,sfdr)
plot(ax,x,sfdr,'-^','Color',[0 .447 .741],'MarkerFaceColor','none','MarkerSize',1.65,'LineWidth',0.55);
plot(ax,x,sndr,'-o','Color',[.85 .325 .098],'MarkerFaceColor','none','MarkerSize',1.65,'LineWidth',0.55);
end

%% Shared multi-chip comparison
function comparison(f,S,x,ticks)
ax=panelaxes(f,S.AxesRect_px);
x=x/S.XUnitDivisor; ticks=ticks/S.XUnitDivisor;
colors=[0 .447 .741;.85 .325 .098;.929 .694 .125];
grey=[0.6470588235 0.7490196078 0.4980392157];
for chip=1:size(S.Y_one_time_percent,2)
    plot(ax,x,S.Y_recal_percent(:,chip),'-o','Color',grey(chip)*ones(1,3),...
        'MarkerSize',1.65,'LineWidth',0.55,'MarkerFaceColor','none');
    plot(ax,x,S.Y_one_time_percent(:,chip),'-o','Color',colors(chip,:),...
        'MarkerSize',1.65,'LineWidth',0.55,'MarkerFaceColor','none');
end
margin=.05*(ticks(end)-ticks(1));
set(ax,'XLim',[ticks(1)-margin ticks(end)+margin],'XTick',ticks,...
    'YLim',S.YLim,'YTick',S.YTick,'YTickLabel',S.YTickLabel);
xlabel(ax,S.XLabel); ylabel(ax,S.YLabel);
fixedlegend(f,S.LegendRect_px,S.ChipLabels,colors,{'o','o','o'},'horizontal');
text(ax,S.RecalLabelPosition(1)/S.XUnitDivisor,S.RecalLabelPosition(2),'w/ recal',...
    'Color',[.65 .65 .65],'HorizontalAlignment','center','FontSize',5.8);
text(ax,S.CalibrationLabelPosition(1)/S.XUnitDivisor,S.CalibrationLabelPosition(2),...
    S.CalibrationLabel,'HorizontalAlignment','center','FontSize',5.8);
end

%% Shared fixed legend geometry
function pairlegend(f,rect)
fixedlegend(f,rect,{'SFDR','SNDR'},[0 .447 .741;.85 .325 .098],{'^','o'},'vertical');
end

function fixedlegend(f,rect,labels,colors,markers,orientation)
ax=axes('Parent',f,'Position',pixelposition(rect),'XLim',[0 rect(3)],...
    'YLim',[0 rect(4)],'YDir','reverse','Visible','off','Clipping','off');
hold(ax,'on');
rectangle(ax,'Position',[0 0 rect(3) rect(4)],'FaceColor','white','EdgeColor',[.16 .16 .16],'LineWidth',0.45);
count=numel(labels);
for k=1:count
    if strcmp(orientation,'horizontal')
        left=(k-1)*rect(3)/count; y=rect(4)*0.5;
        lineX=left+[3 16]; markerX=left+9; textX=left+19;
    else
        y=(k-.5)*rect(4)/count;
        lineX=[3 15]; markerX=8; textX=18;
    end
    plot(ax,lineX,[y y],'-','Color',colors(k,:),'LineWidth',0.55);
    plot(ax,markerX,y,'LineStyle','none','Marker',markers{k},'Color',colors(k,:),...
        'MarkerSize',1.8,'MarkerFaceColor','none','LineWidth',0.55);
    text(ax,textX,y,labels{k},'FontName','Arial','FontSize',5.8,'FontWeight','bold',...
        'Color',[.12 .12 .12],'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','none');
end
end
