%Created by Olga Sosnovtseva 07.06.2022
%Enzyme kinetics: Michaelis-Menten equation
function MM
k1=1;
k_1=1;
k2=1;
n=1; 
Vmax=2;
K=2; 

% set integration time, time linespace and initial conditions
tSpan = [0 30];
tTransient=0;
t = linspace(tTransient, tSpan(2), 1000);

yZero=[10,1,0,0];
%yZero=[10,0.2,1.0,0]
%yZero=[5e-7,2e-7,0,0]
yaZero=[10,0];

% open new figure,  hold on to draw multiple images in one figure
figure(1)
hold on
%cycle for parameter (to see how form of the limit cycle changes)
    %run solver choose ode45, ode23s or smth else; deval it in region of
    %interest

    sol=ode45(@modelMM, tSpan, yZero);
    y = deval(sol,t); 
    
    %draw figures, plot (X,Y)
   plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:));
     xlabel('Time');
    ylabel('Concentration');
   title('Michaelis-Menten equation'); 
 %  end
hold off
 legend('S','E', 'ES', 'P');
 legend('Location', 'NorthEast');

% Approximation
figure(2)
hold on
    sol=ode45(@modelMMapprox, tSpan, yaZero);
    y = deval(sol,t); 
    
    %draw figures, plot (X,Y)
   plot(t,y(1,:),t,y(2,:));
     xlabel('Time');
    ylabel('Concentration');
   title('Michaelis-Menten equation'); 
hold off

%For Hill eq. we do not need to integrate
%since X-asis is concentration of S not time.
figure(3)
hold on
for n=[0.5,1,2,4]
SS=[0:0.1:8];
PP=(Vmax.*SS.^n)./(K.^n+SS.^n);
plot(SS,PP);
end
    xlabel('[S]');
    ylabel('dP/dt');
   title('Hill equation'); 
hold off
legend('n=0.5','n=1', 'n=2', 'n=4');
legend('Location', 'SouthEast');
 
% model itself (diferential equations
function dydt = modelMM(t,y) 
    % diferential equations
    dydt=zeros(4,1);
    S=y(1);
    E=y(2);
    ES=y(3);
    P=y(4);
    dydt(1)=-k1*E*S+k_1*ES;
    dydt(2)=-k1*E*S+k_1*ES+k2*ES;
    dydt(3)= k1*E*S-k_1*ES-k2*ES;
    dydt(4)= k2*ES;
end

function dydt = modelMMapprox(t,y) 
 %  diferential equations
     dydt=zeros(2,1);
     S=y(1);
     P=y(2);
     dydt(1)=-(Vmax*S)/(K+S);
     dydt(2)=(Vmax*S)/(K+S); 
 end

% function dydt = modelHill(t,y) 
%   %  diferential equations
%     dydt=zeros(4,1);
%     S=y(1);
%     E=y(2);
%     ES=y(3);
%     P=y(4);
%     dydt(1)=n*(-k1*E*S^n+k_1*ES);
%     dydt(2)=-k1*E*S^n+k_1*ES+k2*ES;
%     dydt(3)= k1*E*S^n-k_1*ES-k2*ES;
%     dydt(4)= k2*ES;
%     
% end

end

