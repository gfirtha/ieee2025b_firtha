function [fc,A0, phi0] = get_cutoff_circular(x0, xRec,xSource, rSSD, F, A_d, dl, khn, c, win)


[x0Stat,xRefStat] = get_stat_pos(xRec, xSource, rSSD,F);

for ri = 1 :  size(xRec,1)
 %   ri
    [~,stat_ix(ri)] = min(sum((x0-x0Stat(ri,:)).^2,2));
    win_A = A_d./sqrt( sum((x0-xRec(ri,:)).^2,2) )/4/pi;
    win_A = win_A./win_A(stat_ix(ri));
    win_A = circshift(win_A, -stat_ix(ri)+1);
    Leff(ri,:) = ( [sum( win_A(2:end/2) ), sum( win_A(end/2+1:end) )] + win_A(1)/2 )  *dl;

    rG = norm(xRec(ri,:) - x0Stat(ri,:));
    rP = norm(xSource - x0Stat(ri,:));
    d = rP*rG./(F*rP + rG);
    kc(ri,:) = pi/2*d./Leff(ri,:).^2./khn(stat_ix(ri)).^2;
    fc(ri,:) = kc(ri,:)*c/2/pi;
    phi0(ri) = min(-Leff(ri,:).^2.*khn(stat_ix(ri)).^2./d/2,[],2);

end
A0 = win(stat_ix).*sqrt( dist(xRefStat-x0Stat)./dist(xRefStat-xSource).*...
    dist(xRec-xSource)./dist(xRec-x0Stat) );




end