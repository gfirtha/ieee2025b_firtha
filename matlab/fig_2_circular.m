clear
close all

rSSD = 2;         % SSD radius
nSSD = 360*3;       % secondary source number
xSource = [3, 0,0 ];    % source position
xRec = [0,0,0];
c = 343;
dx = 8*0.0025;
taper = 0.25;

phi = (0:nSSD-1)'/nSSD'*2*pi + pi;
dphi = mean(diff(phi));

Lz = 5;
dz = 0.01;
z0 = (-Lz:dz:Lz)';

%% 3D KA
[Phi,Z0] = meshgrid(phi,z0);
X0 = [rSSD*cos(Phi(:)) rSSD*sin(Phi(:)) Z0(:)];
N0 = [ -cos(Phi(:)) -sin(Phi(:)) zeros(numel(Phi),1)];

rG = sqrt(sum((xRec-X0).^2,2));
rP = sqrt(sum((X0-xSource).^2,2));

freq_vec = logspace(0,log10(10e3),256)'; 

k = 2*pi*freq_vec / c;
G0 = 1/4/pi*exp(-1i*k*rG.')./rG.';
P0  = 1/4/pi*exp(-1i*k*rP.')./rP.';

P_ref = 1/4/pi/norm(xRec- xSource);

kP = (X0-xSource)./rP;
kN = sum(kP.*N0,2);
dP_near = (1./rP.*kN).'.*P0;
dP_far = (1i.*k*kN.').*P0;
win = kN>=0;

dA = dphi*rSSD*dz;
P_3D_near = sum( 2*win.'.*dP_near.*G0 , 2)*dA;
P_3D_far = sum( 2*win.'.*dP_far.*G0 , 2)*dA;
P_3D_total = P_3D_near + P_3D_far;
%%
pos = [ 0.145 0.27 0.81 0.7];
ftsize = 16;

f = figure('Units','points','Position',[150,150,370,180]);
p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq_vec, 20*log10(abs(P_3D_near/P_ref)),'LineWidth',2,'Color',[1 1 1]*0.75);
hold on
semilogx(freq_vec, 20*log10(abs(P_3D_far/P_ref)),'LineWidth',2,'Color',[1 1 1]*0.5);
semilogx(freq_vec, 20*log10(abs(P_3D_total/P_ref)),'LineWidth',2,'Color',[1 1 1]*0,'LineStyle','--');
grid on
xlabel('$f$ [Hz]','Interpreter','latex');
ylabel('$|\hat{P}_{\mathrm{3D}}(\mathbf{0},f)|$ [dB]','Interpreter','latex');
ylim([-20,20])

%legend([sl1,sl2],{'$f_{c,1}$','$f_{c,2}$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'fig_2_b_2','-dpng')


%% 2.5D KA
[Phi,Z0] = meshgrid(phi,z0);
x0 = rSSD*[cos(phi) sin(phi) zeros(numel(phi),1)];
n0 = [ -cos(phi) -sin(phi) zeros(numel(phi),1) ];

rG = sqrt(sum((xRec-x0).^2,2));
rP = sqrt(sum((x0-xSource).^2,2));
freq_vec = logspace(0,log10(10e3),256)'; 

k = 2*pi*freq_vec / c;
G0 = 1/4/pi*exp(-1i*k*rG.')./rG.';
P0  = 1/4/pi*exp(-1i*k*rP.')./rP.';

P_ref = 1/4/pi/norm(xRec- xSource);

kP = (x0-xSource)./rP;
kN = sum(kP.*n0,2);
dP_near = (1./rP.*kN).'.*P0;
dP_far = (1i.*k*kN.').*P0;
win = kN>=0;

dref = rP.*rG./(rP+rG);
ds = dphi*rSSD;
P_25D_near = sqrt(8*pi./1i./k).*sum( sqrt(dref.').*win.'.*dP_near.*G0 , 2)*ds;
P_25D_far = sqrt(8*pi./1i./k).*sum( sqrt(dref.').*win.'.*dP_far.*G0 , 2)*ds;
P_25D_total = P_25D_near + P_25D_far;
%%
pos = [ 0.145 0.27 0.81 0.7];
ftsize = 16;

f = figure('Units','points','Position',[150,150,370,180]);
p1 = axes('Units','normalized','Position',pos(1,:));
semilogx(freq_vec, 20*log10(abs(P_25D_near/P_ref)),'LineWidth',2,'Color',[1 1 1]*0.75);
hold on
semilogx(freq_vec, 20*log10(abs(P_25D_far/P_ref)),'LineWidth',2,'Color',[1 1 1]*0.5);
semilogx(freq_vec, 20*log10(abs(P_25D_total/P_ref)),'LineWidth',2,'Color',[1 1 1]*0,'LineStyle','--');
grid on
xlabel('$f$ [Hz]','Interpreter','latex');
ylabel('$|\hat{P}_{\mathrm{2.5D}}(\mathbf{0},f)|$ [dB]','Interpreter','latex');
ylim([-20,20])

%legend([sl1,sl2],{'$f_{c,1}$','$f_{c,2}$'},'FontSize',ftsize,'Interpreter','latex','Location','southeast')

set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'fig_2_c_2','-dpng')

%%
f0 = 1e3;
k0 = 2*pi*f0 / c;

x0 = rSSD*[  cos(phi)  sin(phi) zeros(numel(phi),1) ];
n0 =      [ -cos(phi) -sin(phi) zeros(numel(phi),1) ];

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

kP = (x0-xSource)./rP;
kN = sum(kP.*n0,2);
win = kN>=0;


dref = rhoP.*rhoG./(rhoP+rhoG);
ds = dphi*rSSD;
P0  = 1/4/pi*exp(-1i*k0*rP.')./rP.';
dP_far = (1i.*k0*kN.').*P0;
D_25D = sqrt(8*pi./1i./k0).* sqrt(dref.').*win.'.*dP_far;

R0 = sqrt( (Xfield(:,1)-x0(:,1)').^2 +  (Xfield(:,2)-x0(:,2)').^2  );
G0 = 1/4/pi*exp(-1i*k0*R0)./R0;

Pfield = reshape(  sum( D_25D.*G0.*ds , 2) ,length(y), length(x) );
%%
ri = 2;
LW = 1;
LW2 = 2;
pos = [ 0.067 0.16 0.9 0.82];
ftsize = 14;

fig = figure('Units','points','Position',[150,150,310,250]);
p1 = axes('Units','normalized','Position',pos(1,:));
pcolor(x,y,real(Pfield./P_ref),'FaceAlpha',1);
xl = gca().XLim;
yl = gca().YLim;
shading interp
axis equal tight
hold on

clim([-1,1]);
%clim([-70,0])
xlim(xl);
ylim(yl);
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');

scale_window = 1;
win_A = abs(D_25D)';
win_x = x0(:,1) + win_A/max(abs(win_A)) .* n0(:,1) * scale_window;
win_y = x0(:,2) + win_A/max(abs(win_A)) .* n0(:,2) * scale_window;
fill([x0(:,1); flipud(win_x)], [x0(:,2); flipud(win_y)], [1,1,1]*0.95, ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot(win_x, win_y, 'Color', [100,100,100]/255, 'LineWidth', 1.5);


q = 12;
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.04 );

scatter(xRec(:,1),xRec(:,2),20,"black",'filled')
text(xRec(:,1)-0.5,xRec(:,2)+0.25,'$\mathbf{x} = \mathbf{0}$','FontSize',ftsize,'Interpreter','latex')

cb = colorbar;
cb.Label.String = '[]'; % Or use '$|\hat{P}|/P_{\mathrm{ref}}$'
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = ftsize;

%clim([-70,0])
xlim(xl);
ylim(yl);
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
print( '-r300', 'fig_2_a_2','-dpng')
