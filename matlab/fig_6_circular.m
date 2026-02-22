clear
close all
addpath(genpath('Functions'))
rSSD = 2;         % SSD radius
nSSD = 360*4;       % secondary source number
xSource = [1, 0 ];    % source position
xRec = [0,0;
    0, 1.5];
c = 343;
dx = 8*0.0025;
taper = 0.25;
freq_0  = 1000;

phi = (0:nSSD-1)'/nSSD'*2*pi + pi;
dphi = mean(diff(phi));
x0 = rSSD*[cos(phi) sin(phi)];
n0 = [ -cos(phi) -sin(phi) ];
Rs = norm(xSource);
rhoP = sqrt(sum( (x0-xSource).^2,2 ));
rhoG = sqrt(sum(x0.^2,2));
dl = mean(rhoG*2*pi/size(x0,1));
%
x = (-1.25*rSSD:dx:1.25*rSSD);
y = (-1.25*rSSD:dx:1.25*rSSD);
[X,Y] = meshgrid(x,y);
Xfield = [X(:),Y(:)];
field = zeros(size(Xfield,1),1);
R = sqrt(  (Xfield(:,1)-x0(:,1)').^2 + (Xfield(:,2)-x0(:,2)').^2  );

kh = (x0-xSource)./rhoP;
khn = sum(n0 .* kh,2);
F = -(-1+2*all(khn<0)); % Focused flag
if F == 1
    win = double(khn>=0);
elseif F == -1
    ks = (xSource)/norm(xSource);
    %    win = (kh*ks').^2.*((kh*ks')>0);
    win = (kh*ks').^0.*((kh*ks')>0);
end
[~,max_ix] = max(abs(khn.*win));
win0 = tukeywin(length(find(win==1)),0.5);
win0 = [win0;zeros(length(find(win~=1)),1)];
win0 = circshift(win0,max_ix-round(length(find(win==1))/2)-1);
win = win.*win0;

rRef = rhoG;
dref = sqrt(rRef.*rhoP./(F*rRef+rhoP));
k0 = 2*pi*freq_0/c;
Px0 = 1/(4*pi)*exp(-1i*F*k0*rhoP)./rhoP;
Dx = sqrt(8*pi*1i*k0).*(win.*dref.*khn.*Px0*dl);
Psynth = reshape( sum(Dx.'.*exp(-1i*k0*R)./R/4/pi,2) ,length(y),length(x));
R0 = sqrt( (X-xSource(1)).^2 + (Y-xSource(2)).^2  );
Pref  = 1/(4*pi)*exp(-1i*k0*R0)./R0;

%%
[x0Stat,xRefStat] = get_stat_pos(xRec, xSource, rSSD,F);

for ri = 1 :  size(xRec,1)
    [~,stat_ix(ri)] = min(sum((x0-x0Stat(ri,:)).^2,2));
    win_A = sqrt(8*pi).*abs(win.*dref.*khn.*Px0*dl)./sqrt( sum((x0-xRec(ri,:)).^2,2) )/4/pi;
    win_A = win_A./win_A(stat_ix(ri));
    win_A = circshift(win_A, -stat_ix(ri)+1);
    Leff(ri,:) = ( [sum( win_A(2:end/2) ), sum( win_A(end/2+1:end) )] + win_A(1)/2 )  *dl;

    rG = norm(xRec(ri,:) - x0Stat(ri,:));
    rP = norm(xSource - x0Stat(ri,:));
    d = rP*rG./(F*rP + rG);
    kc(ri,:) = pi/2*d./Leff(ri,:).^2./khn(stat_ix(ri)).^2;
    fc(ri,:) = kc(ri,:)*c/2/pi;

end
A0 = win(stat_ix).*sqrt( dist(xRefStat-x0Stat)./dist(xRefStat-xSource).*...
    dist(xRec-xSource)./dist(xRec-x0Stat) );

freq_vec = logspace(0,log10(10e3),256)'; % Frequency vector (logarithmic scale)

k = 2*pi*freq_vec/c;
Px0f = 1/(4*pi)*exp(-1i*F*k*rhoP')./rhoP';
Dx_f = sqrt(8*pi*1i*k).*(win.*dref.*khn.*Px0f.'*dl).';
P_field_rec = zeros(length(k), size(xRec,1));
for ri = 1 : size(xRec,1)

    rRec = sqrt((x0(:,1)-xRec(ri,1)').^2 + (x0(:,2)-xRec(ri,2)').^2);
    P_field_rec(:,ri) = sum(1/4/pi*Dx_f.*exp(-1i*k*rRec')./rRec', 2).';
    Pref(ri) = A0(ri).*1/4/pi/norm(xRec(ri,:)-xSource);
end

%%
ri = 2;
LW = 1;
LW2 = 2;
pos = [ 0.055 0.14 0.375 0.85;
   0.5325 0.65 0.45 0.325;
   0.5325 0.15 0.45 0.325 ];
ftsize = 12;

fig = figure('Units','points','Position',[150,150,550,250]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(Psynth),'FaceAlpha',0.75);
xl = gca().XLim;
yl = gca().YLim;
shading interp
axis equal tight
hold on
xRef = x0 + F*kh*rSSD;
xRef((win==0)) = nan;
scale_window= 0.5;

win_A = sqrt(8*pi).*abs(win.*dref.*khn.*Px0*dl)./sqrt( sum((x0-xRec(ri,:)).^2,2) )/4/pi;
win_x = x0(:,1) + win_A/max(abs(win_A)) .* n0(:,1) * scale_window;
win_y = x0(:,2) + win_A/max(abs(win_A)) .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [1,1,1]*0.95, ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [100,100,100]/255, 'LineWidth', 1.5);
scale_window = 0.1;
win_eff = zeros(size(x0,1),1);
win_eff( stat_ix(ri)-round(Leff(ri,2)/dl) : stat_ix(ri)+1+round(Leff(ri,1)/dl) ) = 1;
win_effx = x0(:,1) + win_eff .* n0(:,1) * scale_window;
win_effy = x0(:,2) + win_eff .* n0(:,2) * scale_window;

fill([x0(:,1); flipud(win_effx)], [x0(:,2); flipud(win_effy)], [209,151,137]/255, ...
   'FaceAlpha', 0.9, 'EdgeColor',[162,50,52]/255,'EdgeAlpha',0.5);
q = 15;
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.04 );

plot(x0(win_eff~=0,1), x0(win_eff~=0,2), 'Color',[162,50,52]/255, 'LineWidth', LW2)
plot(xRef(:,1),xRef(:,2),'Color','black','LineStyle','--','LineWidth',1.5);
scatter(xSource(1),xSource(2),20,"black",'filled')
scatter(xRec(:,1),xRec(:,2),20,"black",'filled')
scatter(xRefStat(ri,1),xRefStat(ri,2),20,"black",'filled')
scatter(x0Stat(ri,1),x0Stat(ri,2),20,"black",'filled')
line([xSource(1),xRec(ri,1)],[xSource(2),xRec(ri,2)],'color','black','LineWidth',1.5)
line([xRefStat(ri,1),xRec(ri,1)],[xRefStat(ri,2),xRec(ri,2)],'color','black','LineWidth',1.5,'LineStyle',':')
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
clim([-1,1]*1/4/pi/norm(xSource));
%clim([-70,0])
xlim(xl);
ylim(yl);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');

ri = 1;
p1 = axes('Units','normalized','Position',pos(2,:));
semilogx(freq_vec, 20*log10(abs(P_field_rec(:,ri)/Pref(ri))), 'k', 'LineWidth', 2);
ylim([-20,10])
hold on;
yl = gca().YLim;
sl1 = plot([fc(ri,1),fc(ri,1)],yl,'LineWidth',1.5,'LineStyle',':','Color','black');
sl2 = plot([fc(ri,2),fc(ri,2)],yl,'LineWidth',1.5,'LineStyle','--','Color','black');
grid on
xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel(sprintf('$|\\hat{P}_{\\mathrm{synth}}(\\mathbf{x}_{\\mathbf{r}_%d})|$ [dB]', ri), ...
       'FontSize', ftsize, 'Interpreter', 'latex')
legend([sl1,sl2],{'$f_{c,1}$','$f_{c,2}$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);

ri = 2;
p2 = axes('Units','normalized','Position',pos(3,:));
semilogx(freq_vec, 20*log10(abs(P_field_rec(:,ri)/Pref(ri))), 'k', 'LineWidth', 2);
ylim([-20,10])
hold on;
yl = gca().YLim;
sl1 = plot([fc(ri,1),fc(ri,1)],yl,'LineWidth',1.5,'LineStyle',':','Color','black');
sl2 = plot([fc(ri,2),fc(ri,2)],yl,'LineWidth',1.5,'LineStyle','--','Color','black');
grid on
xlabel('$f$ [Hz]', 'FontSize',ftsize,'Interpreter','latex')
ylabel(sprintf('$|\\hat{P}_{\\mathrm{synth}}(\\mathbf{x}_{\\mathbf{r}_%d})|$ [dB]', ri), ...
       'FontSize', ftsize, 'Interpreter', 'latex')
legend([sl1,sl2],{'$f_{c,1}$','$f_{c,2}$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
%print( '-r300', 'focused_vused_fres_unfocq','-dpng')
