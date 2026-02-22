clear
close all
addpath(genpath('Functions\'))

Lssd = 2;                 % Half-length of integration domain in x-direction
dxSSD = 1e-3;

taper = 0.25;
freq = logspace(0,log10(10e3),256)'; % Frequency vector (logarithmic scale)
c = 343;                    % Speed of sound in air (m/s)
f0  = 1000;

xRec = [0, 1.5;
        2,1.5];
xRef = [0,1.5];       % Reference position
xSource = [2,-1];        % Source position
x0 = (-Lssd:dxSSD:Lssd)';        % x-axis spatial vector
y0 = 0;                   % Fixed y-position of the integration line
win = tukeywin(length(x0),taper); % Tukey window applied to x-axis

%% Calculate radiated field
F = sign(y0-xSource(2));
rS = sqrt( (x0-xSource(1)).^2 + (y0-xSource(2)).^2 );

dref = 1./(abs(xRef(2)./(xRef(2)-xSource(2))))./rS;
SPA_corr = sqrt(2*pi./abs(dref)).*exp(-1i*pi/4*sign(dref));
khn = (y0-xSource(2))./rS;
Dx0_prop = 2/(4*pi)*win.*SPA_corr.*khn./rS;

dxMesh = 1e-2;
x = (-2.5:dxMesh:2.5)';
y = (-1:dxMesh:2.5)';
[X,Y] = meshgrid(x,y);

k0 = 2*pi*f0/c;
Dx0 = sqrt(1./k0).*(1i*k0.*Dx0_prop).*exp(-F*1i*k0.*rS);
P_field = zeros(size(X));
for xi = 1  : length(x0)
    progbar(1,length(x0),xi);
    R0 = sqrt((X-x0(xi)).^2 + Y.^2);
    P_field = P_field + 1/4/pi*exp(-1i*k0*R0)./R0.*Dx0(xi)*(dxSSD);
end

%%
R = sqrt((X-xSource(1)).^2 + (Y - xSource(2)).^2 );
P_target = 1./4/pi.*exp(-1i*k0*R)./R;



x0Stat = xSource(1)-(xRec(:,1)-xSource(1))./(xRec(:,2)-xSource(2)).*xSource(2);
xRefStat = (xRec(:,1)-xSource(1))./(xRec(:,2)-xSource(2))*(xRef(2)-xSource(2)) + xSource(1);

%%
trans = 0.8;
scale = 0.2;
LW = 0.75;
LW2 = 2.5;
q= 7;
pos = [ 0.15 0.175 0.8 0.8 ];
ftsize = 13;

ri = 1;
rRec = sqrt((x0-xRec(ri,1)).^2 + (y0-xRec(ri,2)).^2);
A0 = abs(Dx0_prop*1./rRec/4/pi);
[~,ix_stat] = min(abs(x0Stat(ri)-x0));
A0norm = A0/abs(A0(ix_stat));
Leff = [-sum(A0norm(1:ix_stat)) sum(A0norm(ix_stat+1:end))]*dxSSD;

f = figure('Units','points','Position',[150,150,250,200]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(P_field),'FaceAlpha',0.75);
hold on
shading interp
axis equal tight
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')

cax = clim;
hold on
q= 100;
ixs = find(x0<=x(end)&x0>=x(1));
ixs = ixs(1:q:end);
n0 = repmat([0,1],length(x0(ixs)));
draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.03);
off = 0.02;

scale = 0.5;
fill([x0; (x0)], [zeros(size(A0/max(abs(A0)))); A0/max(abs(A0))*scale], [1 1 1]*0.95, ...
    'FaceAlpha', 0.9, 'EdgeColor', 'none');
plot(x0, A0/max(abs(A0))*scale, 'Color',[100,100,100]/255, 'LineWidth', LW2)

clim(cax) % restore caxis


scale = 0.1;
fill([Leff'+x0Stat(ri); flipud(Leff'+x0Stat(ri))], [0;0;scale;scale], [245,154,159]/255, ...
    'FaceAlpha', trans, 'EdgeColor',[209,151,137]/255, 'LineWidth',LW2);
plot(Leff'+x0Stat(ri), [0;0], 'Color',[162,50,52]/255, 'LineWidth', LW2)
clim(cax) % restore caxis


caxis([-1,1]*1/4/pi/(xRef(2)-xSource(2)))



scatter(xSource(1),xSource(2),20,"black",'filled')
scatter(xRec(ri,1),xRec(ri,2),20,"black",'filled')
scatter(x0Stat(ri,1),0,20,"black",'filled')
scatter(xRefStat(ri,1),xRef( 2),20,"black",'filled')
line([xSource(1),xRec(ri,1)],[xSource(2),xRec(ri,2)],'color','black','LineWidth',1.5)


line(gca().XLim,[0,0],'LineWidth',0.1,'Color',[0,0,0]);
line(gca().XLim,[1,1]*xRef(2),'LineWidth',1,'Color',[0,0,0], 'LineStyle','--');

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

set(gcf,'PaperPositionMode','auto');

%%
rRec = sqrt((x0-xRec(:,1)').^2 + (y0-xRec(:,2)').^2);
k = 2*pi*freq'/c;
P_field_rec = zeros(length(k), size(xRec,1));
for ri = 1 : size(xRec,1)
    P_field_rec(:,ri) = sum( 1/4/pi*sqrt(1./k).*(1i*k.*Dx0_prop).*exp(-F*1i*k.*rS).*exp(-1i*rRec(:,ri)*k)./rRec(:,ri)*dxSSD, 1).';
end

pos = [ 0.13 0.27 0.84 0.7];
ftsize = 16;

ri = 1;

Aref = 1/4/pi/norm(xRec(ri,:)-xSource);

f = figure('Units','points','Position',[150,150,400,180]);
p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq, 20*log10(abs(P_field_rec(:,ri)/Aref)), 'k', 'LineWidth', 2);
ylim([-20,5]);


rRec = sqrt((x0-xRec(ri,1)).^2 + (y0-xRec(ri,2)).^2);
A0 = abs(Dx0_prop*1./rRec/4/pi);
[~,ix_stat] = min(abs(x0Stat(ri)-x0));
A0norm = A0/abs(A0(ix_stat));
Leff = [-sum(A0norm(1:ix_stat)) sum(A0norm(ix_stat+1:end))]*dxSSD;
rhoP = norm(xSource-x0(ix_stat));
rhoG = norm(xRec(ri,:)-x0(ix_stat));

d_rec = rhoP*rhoG/(rhoP+F*rhoG);
kc = pi/2*d_rec./(Leff.^2*khn(ix_stat).^2);
fc = abs( kc*c/2/pi );

hold on;
yl = gca().YLim;
sl1 = plot([fc(1),fc(1)],yl,'LineWidth',1.5,'LineStyle',':','Color','black');
sl2 = plot([fc(2),fc(2)],yl,'LineWidth',1.5,'LineStyle','--','Color','black');
%ylim([-10,2])
%xlim([20,10e3])


grid on
xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel(sprintf('$|\\hat{P}_{\\mathrm{synth}}(\\mathbf{x}_{\\mathbf{r}_%d},f)|$ [dB]', ri), ...
       'FontSize', ftsize, 'Interpreter', 'latex')
legend([sl1,sl2],{'$f_{c,1}$','$f_{c,2}$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');

%print( '-r300', sprintf('linear_wfs_transfer_taper_%i',taper0*100) ,'-dpng')
