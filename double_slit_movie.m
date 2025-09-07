clear;

[x y t psi psire psiim psimod v] = ...
sch_2d_adi(0.1/5, 7, 0.01, 1, [0.6, 0.5, 0.05, 0.05, -20, 0], 2, [0.4, 0.45, 0.55, 0.6, 10000000]);



vid = VideoWriter('double_slit.mp4', 'MPEG-4');
open(vid);

for n = 1:length(t)
    figure(1);
    clf;
    ylim([0, 1]);
    xlim([0, 1]);
    clf;
    hold on;
    pcolor(x, y, 10*squeeze(psimod(n, :, :)) + v/(1000*max(v, [], "all")));
    drawnow;
    writeVideo(vid, getframe(gcf));
end

close(vid);
