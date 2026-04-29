% clear, clc, close all

nFrames = 101;
nPointsCircle = 101;

theta1 = x_sem(1,:); theta2 = x_sem(3,:); theta3 = x_sem(5,:); % Sem anti-unwinding
theta4 = x_com(1,:); theta5 = x_com(3,:); theta6 = x_com(5,:); % Com anti-unwinding

% Centros dos círculos
c1 = [-2 0]; c2 = [ 0 0]; c3 = [ 2 0];
c4 = [-2 -3]; c5 = [ 0 -3]; c6 = [ 2 -3];

figure
hold on
axis equal
axis off

thetaCircle = linspace(0,2*pi,nPointsCircle);

% Desenha os 3 círculos
plot(c1(1)+cos(thetaCircle), c1(2)+sin(thetaCircle),'k')
plot(c2(1)+cos(thetaCircle), c2(2)+sin(thetaCircle),'k')
plot(c3(1)+cos(thetaCircle), c3(2)+sin(thetaCircle),'k')

plot(c4(1)+cos(thetaCircle), c4(2)+sin(thetaCircle),'k')
plot(c5(1)+cos(thetaCircle), c5(2)+sin(thetaCircle),'k')
plot(c6(1)+cos(thetaCircle), c6(2)+sin(thetaCircle),'k')

xlim([-3.5 3.5])
ylim([-4.2 1.9])
set(gcf, 'Position', [100 100 1200 800]);

% Criar vídeo
writerObj = VideoWriter('Comparacaonew.mp4','MPEG-4');
writerObj.Quality = 100;   % máximo
writerObj.FrameRate = 6;
open(writerObj);

for i = 1:nFrames
    
    % Sem anti-unwinding
    xa1 = c1(1) + cos(theta1(i));    ya1 = c1(2) + sin(theta1(i)); % Agente 1
    xa2 = c2(1) + cos(theta2(i));    ya2 = c2(2) + sin(theta2(i)); % Agente 2
    xa3 = c3(1) + cos(theta3(i));    ya3 = c3(2) + sin(theta3(i)); % Agente 3

    % Com anti-unwinding
    xa4 = c4(1) + cos(theta4(i));    ya4 = c4(2) + sin(theta4(i)); % Agente 1
    xa5 = c5(1) + cos(theta5(i));    ya5 = c5(2) + sin(theta5(i)); % Agente 2
    xa6 = c6(1) + cos(theta6(i));    ya6 = c6(2) + sin(theta6(i)); % Agente 3
    
    % Sem anti-unwinding
    h1 = plot(xa1,ya1,'o','MarkerSize',8,'Color',[0.39,0.83,0.07],'MarkerFaceColor',[0.39,0.83,0.07]);
    h2 = plot(xa2,ya2,'o','MarkerSize',8,'Color',[0.07,0.62,1.00],'MarkerFaceColor',[0.07,0.62,1.00]);
    h3 = plot(xa3,ya3,'o','MarkerSize',8,'Color',[0.93,0.69,0.13],'MarkerFaceColor',[0.93,0.69,0.13]);

    hline1 = plot([c1(1) xa1],[c1(2) ya1],'Color',[0.39,0.83,0.07],'LineWidth',1.5);
    hline2 = plot([c2(1) xa2],[c2(2) ya2],'Color',[0.07,0.62,1.00],'LineWidth',1.5);
    hline3 = plot([c3(1) xa3],[c3(2) ya3],'Color',[0.93,0.69,0.13],'LineWidth',1.5);

    % Com anti-unwinding
    h4 = plot(xa4,ya4,'o','MarkerSize',8,'Color',[0.39,0.83,0.07],'MarkerFaceColor',[0.39,0.83,0.07]);
    h5 = plot(xa5,ya5,'o','MarkerSize',8,'Color',[0.07,0.62,1.00],'MarkerFaceColor',[0.07,0.62,1.00]);
    h6 = plot(xa6,ya6,'o','MarkerSize',8,'Color',[0.93,0.69,0.13],'MarkerFaceColor',[0.93,0.69,0.13]);

    hline4 = plot([c4(1) xa4],[c4(2) ya4],'Color',[0.39,0.83,0.07],'LineWidth',1.5);
    hline5 = plot([c5(1) xa5],[c5(2) ya5],'Color',[0.07,0.62,1.00],'LineWidth',1.5);
    hline6 = plot([c6(1) xa6],[c6(2) ya6],'Color',[0.93,0.69,0.13],'LineWidth',1.5);

    if i == 1;
    text(-0.7,1.2,'Sem anti-unwinding','Color','k','FontWeight','bold')
    text(-0.7,-1.8,'Com anti-unwinding','Color','k','FontWeight','bold')
    end
    F(i) = getframe;
    writeVideo(writerObj, F(i));
    
    delete(h1); delete(h2);delete(h3); delete(hline1); delete(hline2); delete(hline3)
    delete(h4);delete(h5);delete(h6);delete(hline4); delete(hline5); delete(hline6)

end

close(writerObj);