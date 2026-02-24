clear
close all
addpath(genpath('Functions\'))

Lssd = 2;                 % Half-length of integration domain in x-direction
dxSSD = 1e-3;

taper = 0.25;
freq = logspace(0,log10(10e3),256)'; % Frequency vector (logarithmic scale)
c = 343;                    % Speed of sound in air (m/s)
f0  = 1000;

xRef = [0,1.5];       % Reference position
xSource = [0,-1];        % Source position
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
xRec = [X(:),Y(:)];

k0 = 2*pi*f0/c;
Dx0 = sqrt(1./k0).*(1i*k0.*Dx0_prop).*exp(-F*1i*k0.*rS);

%%

x0Stat = xSource(1)-(xRec(:,1)-xSource(1))./(xRec(:,2)-xSource(2)).*xSource(2);
x0Stat(x0Stat<-Lssd|x0Stat>Lssd) = nan;
xRefStat = (xRec(:,1)-xSource(1))./(xRec(:,2)-xSource(2))*(xRef(2)-xSource(2)) + xSource(1);
%%
for ri = 1 : size(xRec,1)
    ri
    rRec = sqrt((x0-xRec(ri,1)).^2 + (y0-xRec(ri,2)).^2);
    A0 = abs(Dx0_prop*1./rRec/4/pi);
    [~,ix_stat] = min(abs(x0Stat(ri)-x0));
    A0norm = A0/abs(A0(ix_stat));
    Leff = (A0norm(ix_stat)/2 +  [sum(A0norm(1:ix_stat-1)) sum(A0norm(ix_stat+1:end))])*dxSSD;
    rhoP = norm(xSource-[x0(ix_stat),y0]);
    rhoG = norm(xRec(ri,:)-[x0(ix_stat),y0]);

    d_rec = rhoP*rhoG/(rhoP+F*rhoG);
    kc = pi/2*d_rec./(min(abs(Leff)).^2*khn(ix_stat).^2);
    fc(ri,:) = abs( (kc)*c/2/pi );
    if isnan(x0Stat(ri))
        fc(ri,:) = nan;
    end
    if xRec(ri,2)<y0
        fc(ri,:) = nan;
    end
end
%%

R = sqrt((X-xSource(1)).^2 + (Y - xSource(2)).^2 );
P_target = 1./4/pi.*exp(-1i*k0*R)./R;
P_target = P_target / (1/4/pi/norm( xRef-xSource ));

cbar = true;

trans = 0.8;
scale = 0.2;
LW = 0.75;
LW2 = 2.5;
q = 7;
ftsize = 12;

if cbar
    pos = [ 0.135 0.125 0.81 0.9 ];
    f = figure('Units','points','Position',[150,150,330,200]);
else
    pos = [ 0.165 0.125 0.83 0.9 ];
    f = figure('Units','points','Position',[150,150,270,200]);
end
ax1 = axes('Units','normalized','Position',pos(1,:));
pc1 = pcolor(x,y,real(P_target));
pc1.FaceAlpha = 0.5;
hold on
clim([-1,1]*1)

shading interp
axis equal tight
xlabel('$x_{\mathrm{r}}$ [m]', 'FontSize',ftsize)
ylabel('$y_{\mathrm{r}}$ [m]', 'FontSize',ftsize)
set(gca,'FontName','times new roman');
allAxesInFigure = findall(f,'type','axes');
%set(allAxesInFigure, 'TickLabelInterpreter', 'latex');
set(allAxesInFigure, 'FontSize', ftsize);

ax2 = axes('Units','normalized','Position',pos);

pc2 = pcolor( x,y, reshape(fc, length(y), length(x)));
axis equal tight
shading interp
clim([50,1e3])
set(gca, 'ColorScale', 'log');

xl = gca().XLim;
yl = gca().YLim;

% --- Add Contour Lines ---

% --- Add Contour Lines ---
hold(ax2, 'on');
levels = [50, 100, 200, 500];
[C, h] = contour(ax2, x, y, reshape(fc, length(y), length(x)), levels);

% Force these properties explicitly
set(h, 'LineColor', [1 1 1], 'LineWidth', 1.0, 'ShowText', 'on');

% Fix clabel color
cl = clabel(C, h, 'FontSize', 8, 'Interpreter', 'latex', 'Color', [1 1 1]);

% -------------------------
linkaxes([ax1, ax2]);

q= 100;
ixs = find(x0<=x(end)&x0>=x(1));
ixs = ixs(1:q:end);
n0 = repmat([0,1],length(x0(ixs)));
draw_ssd( f, [x0(ixs),0*x0(ixs)] , n0, 0.03);
% Instead of ax2.Visible = 'off' (which can hide the contours during print)
% we hide the axes' internal markings
set(ax2, 'Color', 'none', ...
    'XColor', 'none', ...
    'YColor', 'none', ...
    'XTick', [], ...
    'YTick', []);

if cbar
    % Add the colorbar to ax2
    cb = colorbar(ax2);
    cb.Label.String = '[Hz]';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = ftsize;

    % MANDATORY STEP: Re-align the axes.
    % Adding the colorbar resized ax2. We must now set ax1 to match ax2.
    newPos = get(ax2, 'Position'); % Get the shrunk position of ax2
    set(ax1, 'Position', newPos);  % Apply it to ax1
else
end

% Set the fonts/interpreters for both axes consistently
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure, 'TickLabelInterpreter', 'latex'); % Use latex for the numbers
set(allAxesInFigure, 'FontName', 'Times New Roman', 'FontSize', ftsize);
% Ensure the figure background is preserved
set(f, 'InvertHardcopy', 'off', 'Color', 'w');
set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
% Print with OpenGL to preserve overlapping colors and line styles
print(f, 'Fc2_lin_vs_recpos_ps2', '-dpng', '-r300', '-opengl');

