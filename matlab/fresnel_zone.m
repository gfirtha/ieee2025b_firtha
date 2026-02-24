clear
close all

rSSD = 2;         % SSD radius
nSSD = 360*4;       % secondary source number
xSource = [3, 0 ];    % source position
xRec = [0,0;
    0, 1.5];
c = 343;
dx = 8*0.0025;
taper = 0.25;

phi = (0:nSSD-1)'/nSSD'*2*pi + pi;
dphi = mean(diff(phi));
x0 = rSSD*[cos(phi) sin(phi)];
n0 = [ -cos(phi) -sin(phi) ];
Rs = norm(xSource);
rhoP = sqrt(sum( (x0-xSource).^2,2 ));
rhoG = sqrt(sum(x0.^2,2));
dl = mean(rhoG*2*pi/size(x0,1));
%
x = (-5.25*rSSD:dx:1.1*rSSD);
y = (-1.1*rSSD:dx:1.1*rSSD);
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
Px0 = 1/(4*pi)./rhoP;
A_d = sqrt(8*pi).*abs(win.*dref.*khn.*Px0*dl);

[fc,A0,phi0] = get_cutoff_circular(x0, Xfield ,xSource, rSSD, F, A_d , dl, khn,c, win);
fc = max(fc,[],2);


%% Plotting with Filled Contours

%% Discrete Banded Plotting Fix
ftsize = 12;
levels = [50, 100, 200, 500, 1000, 2000];
FC_grid = reshape(fc, length(y), length(x));

% 1. Create a "Binned" version of your data
% This ensures that values between 100-500 get value 1, 500-1000 get value 2, etc.
binned_fc = zeros(size(FC_grid));
for i = 1:length(levels)
    binned_fc(FC_grid >= levels(i)) = i;
end

fig = figure('Color', 'w');
hold on;

% 2. Plot the binned data using contourf
% We plot levels 1 through 5
[C, h] = contourf(x, y, binned_fc, 0:length(levels), 'LineStyle', 'none');

% 3. Define the Exact Discrete Colormap
% Each row is the color for one of your frequency bands
myMap = [
    0.750, 0.750, 0.750;  % 0 - 100 Hz (White/Background)
    0.70, 0.85, 1.00;  % 100 - 500 Hz (Light Blue)
    1.00, 0.95, 0.70;  % 500 - 1000 Hz (Pale Yellow)
    1.00, 0.80, 0.50;  % 1000 - 2000 Hz (Orange)
    0.95, 0.60, 0.50;  % 2000 - 5000 Hz (Coral/Red)
    0.75, 0.60, 0.85   % > 5000 Hz (Purple)
];
colormap(myMap);
clim([0, length(levels)]); % Lock the color limits to our bins

% 4. Add the white lines and labels using the ORIGINAL fc data
[C2, h2] = contour(x, y, FC_grid, levels, 'LineColor', 'k', 'LineWidth', 1.2);
clabel(C2, h2, 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold');

axis equal tight;
grid on;
set(gca, 'Layer', 'top');

% 6. Custom Colorbar
cb = colorbar;
cb.Ticks = 0.5:1:length(levels)-0.5;
cb.TickLabels = {'<100', '100-500', '500-1k', '1k-2k', '2k-5k', '>5k'};
ylabel(cb, 'Cutoff Frequency $f_c$ [Hz]', 'Interpreter', 'latex');

q = 15;
draw_ssd( fig , x0(1:q:end,:) + n0(1:q:end,:)*0.02 , n0(1:q:end,:), 0.04 );


set(gca, 'ColorScale', 'log');
xlabel('$x_r$ [m]', 'FontSize',ftsize,'Interpreter','latex')
ylabel('$y_r$ [m]', 'FontSize',ftsize,'Interpreter','latex')


set(gca,'FontName','Times New Roman');
allAxesInFigure = findall(fig,'type','axes');
set(allAxesInFigure,'FontSize',ftsize);
set(gcf,'PaperPositionMode','auto');
