%Created by Olga Sosnovtseva 07.06.2022
%Logistic equation (population dynamics)
function logistic_eq
r=0.1;
K=2000;

% Set integration time, time linespace and initial conditions
tSpan = [0 200];
tTransient=0;
t = linspace(tTransient, tSpan(2), 1000);
yZero = [100];

% open new figure,  hold on to draw multiple images in one figure
 figure (1)
 hold on
%cycle for parameter (to see how form of the limit cycle changes)
    %run solver choose ode45, ode23s or smth else; deval it in region of
    %interest
%for r=0.01:0.01:0.1 % How does r parameter affect dynamics
%for yZero= 100:250:2500 %How do trajectories approach equilibrium from
%different initial conditions
    sol=ode45(@model, tSpan, yZero);
    y = deval(sol,t); 
    
%draw figures, plot (X,Y)
     plot(t,y(1,:));
      xlabel('t');
      ylabel('N');
 %     title('Logistic equation');
% end
 
%hold off % to stop drawing on this figure


%legend('N(0)=100','N(0)=350', 'N(0)=600');
%legend('r=0.01','r=0.02', 'r=0.03');
%legend('Location', 'SouthEast');

%figure (2)
% hold on
%Plot vector field. L is normalization factor
 [T,Y] = meshgrid(0:10:200,-50:50:2500);
 dY =r*Y.*(1-Y/K);
 dT=ones(size(dY));
 L=sqrt(dY.^2+dT.^2);
 quiver(T,Y,dT./L,dY./L,0.6,'r');
% 
% %quiver(T,Y,dT,dY,0.6,'r');
% 
 hold off


% model itself (diferential equations
function dydt = model(t,y) 
%    diferential equations
dydt=zeros(1,1);
dydt(1)=r*y(1)*(1-y(1)/K);
end

end