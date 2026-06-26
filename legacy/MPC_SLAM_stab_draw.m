function MPC_SLAM_stab_draw (t,xx,xx1,u_cl,xs,N,rob_diam,xx_true,u_cl_SLAM)


set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;

%--------------------------------------------------------------------------
%-----------------------Simulate robots -----------------------------------
%--------------------------------------------------------------------------
x_r_1 = [];
y_r_1 = [];
x_r_1_true = [];
y_r_1_true = []; 


r = rob_diam/2;  % obstacle radius
ang=0:0.005:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);

figure(500)
% Animate the robot motion
%figure;%('Position',[200 200 1280 720]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

for k = 1:size(xx,2)
    h_t = 0.14; w_t=0.09; % triangle parameters
    
    x1 = xs(1); y1 = xs(2); th1 = xs(3); %ref
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot reference state
    hold on;
    x1 = xx(1,k,1); y1 = xx(2,k,1); th1 = xx(3,k,1);
    x_r_1 = [x_r_1 x1];
    y_r_1 = [y_r_1 y1];






    %true state
    x1_true = xx_true(1,k,1); y1_true = xx_true(2,k,1); th1_true = xx_true(3,k,1);    
    x1_tri_true = [ x1_true+h_t*cos(th1_true), x1_true+(w_t/2)*cos((pi/2)-th1_true), x1_true-(w_t/2)*cos((pi/2)-th1_true)];%,x1+(h_t/3)*cos(th1)];
    y1_tri_true = [ y1_true+h_t*sin(th1_true), y1_true-(w_t/2)*sin((pi/2)-th1_true), y1_true+(w_t/2)*sin((pi/2)-th1_true)];%,y1+(h_t/3)*sin(th1)];
    x_r_1_true = [x_r_1_true x1_true];
    y_r_1_true = [y_r_1_true y1_true];



    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];

    plot(x_r_1,y_r_1,'-r','linewidth',line_width);hold on % plot exhibited trajectory
    plot(x_r_1_true,y_r_1_true,'-k','linewidth',line_width);
    if k < size(xx,2) % plot prediction
         plot(xx1(1:N,1,k),xx1(1:N,2,k),'r--*')
%          plot(xx_true(1:N,1,k),xx_true(1:N,2,k),'k--*')
    end
    



    fill(x1_tri_true, y1_tri_true, 'k'); % plot robot position
    plot(x1_true+xp,y1_true+yp,'--k'); % plot robot circle

    fill(x1_tri, y1_tri, 'r'); % plot robot position
    plot(x1+xp,y1+yp,'--r'); % plot robot circle
    
   
    hold off
    %figure(500)
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    legend("reference position","estimated position","true position")
    axis([-0.2 1.8 -0.2 1.8])
    pause(0.1)
    box on;
    grid on
    %aviobj = addframe(aviobj,gcf);
    drawnow
    % for video generation
    F(k) = getframe(gcf); % to get the current frame
end



% close(gcf)
%viobj = close(aviobj)
%video = VideoWriter('exp.avi','Uncompressed AVI');

% video = VideoWriter('exp.avi','Motion JPEG AVI');
% video.FrameRate = 5;  % (frames per second) this number depends on the sampling time and the number of frames you have
% open(video)
% writeVideo(video,F)
% close (video)

figure(300)
subplot(211)
stairs(t,u_cl(:,1),'k','linewidth',1.5); axis([0 t(end) -0.35 0.75])
hold on
stairs(t,u_cl_SLAM(:,1),'r','linewidth',1.5); axis([0 t(end) -0.35 0.75])
ylabel('v (rad/s)')
grid on
subplot(212)
stairs(t,u_cl(:,2),'k','linewidth',1.5); axis([0 t(end) -1.5 1.5])
hold on
stairs(t,u_cl_SLAM(:,2),'r','linewidth',1.5); axis([0 t(end) -1.5 1.5])
xlabel('time (seconds)')
ylabel('\omega (rad/s)')
grid on
