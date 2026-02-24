clear
close all

%% Parameters
Lssd = 2;                 % Half-length of integration domain in x-direction
dxSSD = 1e-3;
taper = 0.25;
freq = logspace(0,log10(10e3),256)'; % Frequency vector (logarithmic scale)
c = 343;                    % Speed of sound in air (m/s)
f0  = 1000;
xRef = [0,1.5];             % Reference position
xSource = [0,-1e3];           % Source position
x0 = (-Lssd:dxSSD:Lssd)';   % SSD positions
y0 = 0;                     % Fixed y-position of the SSD line
win = tukeywin(length(x0),taper); % Tukey window applied to SSD

%% Pre-calculate Geometry and Driving Functions
F = sign(y0-xSource(2));
rS = sqrt( (x0-xSource(1)).^2 + (y0-xSource(2)).^2 );
dref = 1./(abs(xRef(2)./(xRef(2)-xSource(2))))./rS;
SPA_corr = sqrt(2*pi./abs(dref)).*exp(-1i*pi/4*sign(dref));
khn_all = (y0-xSource(2))./rS;
Dx0_prop = 2/(4*pi)*win.*SPA_corr.*khn_all./rS;

%% Receiver Grid Setup
dxMesh = 1e-2;
x = (-2.5:dxMesh:2.5)';
y = (-1:dxMesh:2.5)';
[X,Y] = meshgrid(x,y);
xRec = [X(:),Y(:)];
numRec = size(xRec,1);

k0 = 2*pi*f0/c;
R0 = sqrt((X-xSource(1)).^2 + (Y - xSource(2)).^2 );
P_target = 1./4/pi.*exp(-1i*k0*R0)./R0;
P_target = P_target / (1/4/pi/norm( xRef-xSource ));

% Calculate stationary positions for all receivers
x0Stat = xSource(1)-(xRec(:,1)-xSource(1))./(xRec(:,2)-xSource(2)).*xSource(2);

%% Vectorized Cutoff Frequency Calculation
% 1. Distances from all SSD elements to all receivers (Matrix: numRec x numSSD)
% Use broadcasting to avoid explicit loops
dist_rec_ssd = sqrt((xRec(:,1) - x0').^2 + (xRec(:,2) - y0).^2); 

% 2. Stationary point mapping
[~, ix_stat_all] = min(abs(x0Stat - x0'), [], 2);

% 3. Integrand Amplitude Matrix A0
% A0(ri, si) = 1 / (dist_to_source * dist_to_receiver)
A0_mat = (1 ./ rS') .* (1 ./ dist_rec_ssd);

% 4. Normalize by stationary point amplitude
stat_indices = sub2ind(size(A0_mat), (1:numRec)', ix_stat_all);
A0_stat_vals = A0_mat(stat_indices);
A0norm_mat = A0_mat ./ A0_stat_vals;

% 5. Vectorized Integration for Effective Length (Leff)
% Using cumulative sums to calculate integrals on both sides of stationary point
A0_cumsum = cumsum(A0norm_mat, 2);
total_sum = A0_cumsum(:, end);

% Sum from index 1 to ix_stat-1
sum_left = zeros(numRec, 1);
idx_not_first = ix_stat_all > 1;
sum_left(idx_not_first) = A0_cumsum(sub2ind(size(A0_mat), find(idx_not_first), ix_stat_all(idx_not_first)-1));

% Sum from index ix_stat+1 to end
sum_right = total_sum - A0_cumsum(stat_indices);

% Apply the Trapezoidal correction (A0(stat)/2 + sums) * dx
Leff1 = (0.5 + sum_left) * dxSSD;
Leff2 = (0.5 + sum_right) * dxSSD;
Leff_min = min(Leff1, Leff2);

% 6. Physical Geometry at Stationary Points
rhoP_stat = rS(ix_stat_all);
rhoG_stat = sqrt((xRec(:,1) - x0(ix_stat_all)).^2 + (xRec(:,2) - y0).^2);
d_rec_vec = (rhoP_stat .* rhoG_stat) ./ (rhoP_stat + F * rhoG_stat);

% 7. Final Cutoff Frequency calculation
kc_vec = (pi/2 * d_rec_vec) ./ (Leff_min.^2 .* khn_all(ix_stat_all).^2);
fc_raw = abs((kc_vec) * c / 2 / pi);

% 8. Apply Validity Masks
valid_mask = ~isnan(x0Stat) & (x0Stat >= -Lssd & x0Stat <= Lssd) & (xRec(:,2) > y0);
fc = fc_raw;
fc(~valid_mask) = nan;
fc_reshaped = reshape(fc, length(y), length(x));

%% Visualization
ftsize = 13;
pos = [ 0.15 0.175 0.8 0.8 ];
f = figure('Units','points','Position',[150,150,250,200]);

% Layer 1: Target Wave Field
ax1 = axes('Units','normalized','Position',pos);
pc1 = pcolor(x,y,real(P_target));
set(pc1, 'FaceAlpha', 0.5);
shading interp
axis equal tight
xlabel('$x$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y$ [m]', 'FontSize',ftsize,'Interpreter','latex')
hold on

% Layer 2: Cutoff Frequency Color Map
ax2 = axes('Units','normalized','Position',pos);
pc2 = pcolor(x, y, fc_reshaped);
shading interp
axis equal tight
clim([50, 1000]);
set(ax2, 'ColorScale', 'log', 'Visible', 'off');

% Layer 3: Contour Lines for Cutoff
hold(ax2, 'on');
levels = [50, 100, 200, 500];
[C, h] = contour(ax2, x, y, fc_reshaped, levels, ...
    'LineColor', [0.2 0.2 0.2], 'LineWidth', 1.0, 'ShowText', 'on');
clabel(C, h, 'FontSize', 8, 'Interpreter', 'latex', 'Color', [0.2 0.2 0.2]);

% Link and Clean up
linkaxes([ax1, ax2]); 
ax2.XTick = []; ax2.YTick = [];

% Draw SSD elements
q = 100;
ixs = find(x0 <= x(end) & x0 >= x(1));
ix_plot = ixs(1:q:end);
n0 = repmat([0,1], length(ix_plot), 1);
draw_ssd(f, [x0(ix_plot), zeros(size(ix_plot))], n0, 0.03);

set(findall(f,'type','axes'), 'FontSize', ftsize, 'FontName', 'Times New Roman');