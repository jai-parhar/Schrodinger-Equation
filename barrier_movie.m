clear;

[x y t psi psire psiim psimod v] = ...
sch_2d_adi(0.1, 7, 0.05, 1, [0.5, 0.5, 0.1, 0.1, -5, 0], 1, [0.1, 0.2, 0.25, 0.75, 1000]);



vid = VideoWriter('barrier.mp4', 'MPEG-4');
open(vid);

for n = 1:length(t)
    figure(1);
    clf;
    ylim([0, 1]);
    xlim([0, 1]);
    clf;
    hold on;
    pcolor(x, y, squeeze(psimod(n, :, :)));
    rectangle("Position",[0.1, 0.25, 0.1, 0.5], "LineWidth",3, "EdgeColor",'r')
    drawnow;
    writeVideo(vid, getframe(gcf));
end

close(vid);
