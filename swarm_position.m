% clc
clear all
Z = load('att_rep_position.dat');
N= 100;
iteration = 50000;
k = iteration/100+1;
Z1 = reshape(Z,N*3,k);
% figure
for i = k
    Z2 = Z1(:,i);
    Z3 = reshape(Z2,3,N);
    Z3 = Z3';
    scatter(Z3(:,1),Z3(:,2),60,mod(Z3(:,3),2*pi),'filled')
    xlabel('x','fontsize',20)
    ylabel('y','fontsize',20)
    caxis([0 2*pi]);
    colorbar('Ticks',[0,pi/2,pi,3*pi/2,2*pi],...
         'TickLabels',{'0','\pi/2','\pi','3\pi/2','2\pi'})
     set(gca,'linewidth',2,'Fontsize',20)
%     title('J = 1.0, K = 1.0, r = 0.2, N = 1000')
%     axis([-5 5 -5 5]);
    pause(0.1);
end
