%Created by Olga Sosnovtseva 10.08.20222
% Introduction to modeling. Template for Van der Pol equations (templateVDP).
% This template can be used for any model in the course.

function templateVDP
%parameters of the model
Mu=1000; %for stiff Van der Pol model, standard value for non-stiff model is 1

% set integration time, time linespace and initial conditions
tSpan = [0 3000]; %calculation time
tTransient=0; % to remove trainsient interval from plots
nPoints=50; % number of steps; the higher, the smoother the curve will be
t = linspace(tTransient, tSpan(2), nPoints);
yZero = [2;0]; % initial conditions (the number is elements depends on the number of ODEs)

% open new figure,  hold on to draw multiple images in one figure, e.g. "for loop"

figure
hold on

%run solver choose ode45, ode15s or smth else;

% If you use ODE function to solve ODE, t will be unevenly distributed which depends on a solver %disp(t)
% If you use deval function to evaluate the solution "sol" at particular "t" from linspace. "t" will be evenly distributed.

% You can do simulation for multiple values of parameters with "for loop"
% for Mu=0.01:0.5:4.0

%TASK 4
%t are unevenly distributed if you use ODE solver
%[t,y]=ode45(@model, tSpan, yZero);
[t,y]=ode15s(@model, tSpan, yZero);
 plot(t,y(:,1),'-o'); % the first column of the matrix 

%TASK 5
%t is evenly distributed when you add deval function
%sol=ode15s(@model, tSpan,yZero);
%  y = deval(sol,t); 
%  plot(t, y(1,:),'-o'); % all elements in the first row
% disp(t)
  
  %end %"for loop"

% Differential equations for the model (Van der Pol model as in https://en.wikipedia.org/wiki/Van_der_Pol_oscillator)
function dydt = model(~,y)% one may  model(t,y) if there is a funtion of t, e.g. sin(w*t)
    dydt=zeros(2,1);
    xx=y(1);  % xx and yy are variables of Van der Pol system. It is good practice to use the same names as in the sourse where this model comes from
    yy=y(2);  % y(1) and y(2) are elements of the vector y that will be used by a solver and where we write solutions
    dydt(1)=yy;
    dydt(2)=Mu*(1-xx^2)*yy-xx;
end    
end
