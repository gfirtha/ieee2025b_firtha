clear
close all


rSSD = 2;         % SSD radius
nSSD = 360*4;       % secondary source number
xRef = [0,0];
c = 343;
dx = 8*0.0025;
taper = 0.25;

phi = (0:nSSD-1)'/nSSD'*2*pi + pi;
dphi = mean(diff(phi));
x0 = rSSD*[cos(phi) sin(phi)];
n0 = [ -cos(phi) -sin(phi) ];

rhoG = sqrt(sum(x0.^2,2));
dl = mean(rhoG*2*pi/size(x0,1));
%
x = (-1.25*rSSD:dx:1.25*rSSD);
y = (-1.25*rSSD:dx:1.25*rSSD);
[X,Y] = meshgrid(x,y);
Xfield = [X(:),Y(:)];
field = zeros(size(Xfield,1),1);
R = sqrt(  (Xfield(:,1)-x0(:,1)').^2 + (Xfield(:,2)-x0(:,2)').^2  );

for n = 1 : size(Xfield,1)
    n
    xSource = Xfield(n,:);
    Rs = norm(xSource);
    rhoP = sqrt(sum( (x0-xSource).^2,2 ));

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
    win0 = tukeywin(length(find(win==1)),taper);
    win0 = [win0;zeros(length(find(win~=1)),1)];
    win0 = circshift(win0,max_ix-round(length(find(win==1))/2)-1);
    win = win.*win0;
    rRef = rhoG;
    dref = sqrt(rRef.*rhoP./(F*rRef+rhoP));
    Px0 = 1/(4*pi)./rhoP;
    A_d = sqrt(8*pi).*abs(win.*dref.*khn.*Px0*dl);

    xRec = xRef;
    [x0Stat,xRefStat] = get_stat_pos(xRec, xSource, rSSD,F);

    [~,stat_ix] = min(sum((x0-x0Stat).^2,2));
    win_A = A_d./sqrt( sum((x0-xRec).^2,2) )/4/pi;
    win_A = win_A./win_A(stat_ix);
    win_A = circshift(win_A, -stat_ix+1);
    Leff = ( [sum( win_A(2:end/2) ), sum( win_A(end/2+1:end) )] + win_A(1)/2 )  *dl;

    rG = norm(xRec - x0Stat);
    rP = norm(xSource - x0Stat);
    d = rP*rG./(F*rP + rG);
    kc = pi/2*d./Leff.^2./khn(stat_ix).^2;
    fc = kc*c/2/pi;
    fc_vec(n) = max(fc,[],2);
end

%%
f0 = 1e3;
k0 = 2*pi*f0/c;
R0 = sqrt( (Xfield(:,1)-xSource(1)).^2 + (Xfield(:,2)-xSource(2)).^2  );
Pref = 1/4/pi*exp(-1i*k0*R0)./R0;

%%
crit_idcs = [find(win0>0,1,'first'), find(win0==1,1,'first') , find(win0==1,1,'last')  ,find(win0>0,1,'last')];

cbar = false;
t = 1e3;

LW = 1;
LW2 = 2;
ftsize = 13;


if cbar
    pos = [ 0.12 0.165 0.81 0.81 ];
    f = figure('Units','points','Position',[150,150,280,220]);
else
    pos = [ 0.16 0.165 0.81 0.81 ];
    f = figure('Units','points','Position',[150,150,220,220]);
end
ax1 = axes('Units','normalized','Position',pos(1,:));
pc1 = pcolor( x,y, reshape(fc_vec, length(y), length(x)));
pc1.FaceAlpha = 1;
shading interp;
xlabel('$x_{\mathrm{r}}$ [m]', 'FontSize',ftsize)
ylabel('$y_{\mathrm{r}}$ [m]', 'FontSize',ftsize)
clim(ax1, [-1, 1]*1/4/pi./norm(xSource(1))); % Skála a hullámnak
fc((sum((Xfield-xSource).*Xfield,2)<0)&(sqrt(sum(Xfield.^2,2))>rSSD) | Xfield(:,1)>xSource(1)) = nan;
hold on


clim([50,1e3])
set(gca, 'ColorScale', 'log');
xlabel('$x_r$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y_r$ [m]', 'FontSize',ftsize,'Interpreter','latex')

% --- Add Contour Lines ---
hold(ax1, 'on');
levels = [60, 100, 500, 1000];
[C, h] = contour(ax1, x, y, reshape(fc_vec, length(y), length(x)), levels, ...
    'LineColor', [1 1 1], ... % Dark gray for visibility
    'LineWidth', 1.0, ...
    'ShowText', 'on', ...           % Optional: labels the lines with the frequency
    'LabelSpacing', 200);           % Adjust text spacing
clabel(C, h, 'FontSize', 8, 'Interpreter', 'latex', 'Color', [1 1 1]);
% -------------------------
q = 15;
draw_ssd( f , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.04 );


%xlim(xl);
%ylim(yl);

if cbar
    % 1. Add colorbar to ax2
    cb = colorbar(ax1);
    cb.Label.String = '[Hz]';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = ftsize;
    
    % 2. Force MATLAB to calculate the new layout immediately
    drawnow; 
    
    % 3. Use 'normalized' units for both to ensure a match
    set(ax1, 'Units', 'normalized');
    set(ax2, 'Units', 'normalized');
    
    % 4. Synchronize the positions
    newPos = get(ax2, 'Position'); 
    set(ax1, 'Position', newPos);
    
    % 5. Optional: If using 'axis equal', force it on both to prevent drift
    axis(ax1, 'equal');
    axis(ax2, 'equal');
end
% Set the fonts/interpreters for both axes consistently
allAxesInFigure = findall(f,'type','axes');
set(allAxesInFigure, 'TickLabelInterpreter', 'latex'); % Use latex for the numbers
set(allAxesInFigure, 'FontName', 'Times New Roman', 'FontSize', ftsize);
% Ensure the figure background is preserved
set(f, 'InvertHardcopy', 'off', 'Color', 'w');
set(gca,'FontName','Times New Roman');

set(gcf,'PaperPositionMode','auto');
% Print with OpenGL to preserve overlapping colors and line styles
%print(f, 'Fc_circ_fps2', '-dpng', '-r300', '-opengl');

