[x M Y] = sim(150,-0.6,-0.1,-0.4,0.7,0.7,10);
[M_in Y_in eff pmat] = M6(x,M,Y);

% p-value matrix plot for original and shulffled
figure;imagesc(-log(pmat));colormap jet; colorbar;





