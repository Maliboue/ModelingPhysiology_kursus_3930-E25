%Created by Olga Sosnovtseva 09.06.2022
%van der Pol oscillator
% Some examples are here http://www.stat.harvard.edu/Research/matlab/test3.html

function vanderpol
% Set initial Mu 
Mu=-0.5; %for stable spiral
%Mu=0.0; %for center
%Mu=-5; %for stable node
%Mu=1.0; %for cycle
A=0;
%for chaos
% Mu=8.53;
% A=1.2; 
 W=2*pi/10;

% set integration time, time linespace and initial conditions
tSpan = [0 500];
%(to plot spiral, we set transient to 0,otherwise 490)
tTransient=0;
%tTransient=500
t = linspace(tTransient, tSpan(2), 10000);
yZero= [-1;2]; % for all tasks

% open new figure,  hold on to draw multiple images in one figure
figure
hold on
%cycle for mu (to see how form of the limit cycle changes)
%for Mu=0.01:0.5:4.0
for par=-1.0:0.2:0 %for spiral and node at different initial conditions
yZero(1)=par;
%for par=0:0.5:2 %for cycle at different initial conditions
%    yZero(2)=par;
%run solver choose ode45, ode23s or smth else; deval it in region of
    %interest
    sol=ode45(@model, tSpan, yZero);
    y = deval(sol,t); 
    %draw figures, plot (X,Y)
    plot(y(1,:),y(2,:));
%   plot(t,y(1,:));
 xlabel('x');
ylabel('y');
%     title('van der Pol oscillator. Chaos.');  
end
% hold off to stop drawing on this figure

%Plot vector field
[X,Y] = meshgrid(-3:0.25:3,-3:0.25:3);
dX = Y;
dY = Mu*(1 - X.^2).*Y - X;
L=sqrt(dX.^2+dY.^2);
quiver(X,Y,dX./L,dY./L,0.6,'k');


hold off

% model itself (diferential equations)
function dydt = model(t,y)
    dydt=zeros(2,1);
    xx=y(1);  
    yy=y(2);
    dydt(1)=yy;
    dydt(2)=Mu*(1-xx^2)*yy-xx + A*sin(W*t);
 
end

end
