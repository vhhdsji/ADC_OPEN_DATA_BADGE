function plot_fft(ax,frequency_MHz,fft_db,config,sndr_db,sfdr_db,s)
plot(ax,frequency_MHz,fft_db,'Color',s.blue,'LineWidth',s.trace)
hold(ax,'on')
set(ax,'Color','w','FontName',s.font,'FontSize',s.tick,...
    'FontWeight','bold','LineWidth',s.width,'TickDir','in',...
    'TickLength',[0.012 0.012],'XGrid','on','YGrid','on',...
    'GridColor',s.grid,'GridAlpha',0.35,'Box','on','Layer','top')
pbaspect(ax,[1.14 1 1])
xlim(ax,[0 2000])
ylim(ax,[-120 0])
xticks(ax,0:500:2000)
yticks(ax,-120:20:0)
xlabel(ax,'Frequency (MHz)','FontSize',s.label,'FontWeight','bold')
ylabel(ax,'Amplitude (dBFS)','FontSize',s.label,'FontWeight','bold')
N=8192;
for h=2:3
    b=mod(h*config.tone_bin,N);
    if b>N/2
        b=N-b;
    end
    xv=b/N*4e9/1e6;
    yv=fft_db(b+1);
    plot(ax,xv,yv,'s','MarkerSize',4.5,'LineWidth',0.8,...
        'MarkerEdgeColor','r','MarkerFaceColor','w','HandleVisibility','off')
    text(ax,xv,min(yv+6,-5),sprintf('%d',h),'HorizontalAlignment','center',...
        'FontName',s.font,'FontSize',s.text,'FontWeight','bold')
end
t=sprintf(['8192 pts FFT\nFs = 4 GS/s\nFin = %.1f MHz\n'...
    'SNDR = %.1f dB\nSFDR = %.1f dB'],...
    config.Fin_Hz/1e6,sndr_db,sfdr_db);
if config.Fin_Hz<1e9
    p=[0.455 0.575 0.445 0.385];
else
    p=[0.425 0.575 0.445 0.385];
end
q=get(ax,'Position');
p=[q(1)+p(1)*q(3) q(2)+p(2)*q(4) p(3)*q(3) p(4)*q(4)];
annotation(ancestor(ax,'figure'),'textbox',p,'String',t,...
    'FontName',s.font,'FontSize',s.text,'FontWeight','bold',...
    'LineWidth',1,'BackgroundColor','w','FitBoxToText','off',...
    'Margin',3,'VerticalAlignment','top')
end
