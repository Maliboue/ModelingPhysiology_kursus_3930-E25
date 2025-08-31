%Olga Sosnovtseva, 13.06.2022
%CICR IP3-receptor based model
function calciumDynamics

IP3=1.5e-3;
A=4.0;
Kip3=0.65e-3;
Xact=4.0;
Kact=1.3e-4;
Xinact=4.0;
Kinact=3.5e-4;
Kca=0.0526;
Rip3=7.5;
Rleak=2.0;
Rserca=0.170979;
P=2.5;
Kserca=0.7e-4;
Tinact=6.25;
K=0.2;

% set integration time, time linespace and initial conditions
tSpan = [0 150];
tTransient=0;
t = linspace(tTransient, tSpan(2), 5000);
yZero = [1e-4;0];

%run solver choose ode45, ode23s or smth else; deval it in region of
%interest
figure(3)
hold on
% for Tinact=5:0.5:10 % Task 4, changes with inaction time 
%for IP3=1.0e-3:5e-4:10e-3 %for saturation % Task 3
for IP3=1.0e-3:1e-5:1.5e-3 % Task 2, bifurcation diagram
sol=ode23s(@model, tSpan, yZero);
y = deval(sol,t);
 YMin=real(min(y(1,:)));
 YMax=real(max(y(1,:)));
 [pks,locs] = findpeaks(y(1,:));
 plot(IP3,(YMax-YMin)/2,'bo' ); %for bifurcation diagram and saturation
%plot(Tinact,(YMax-YMin)/2,'ro' ); %for Task 4
%xlabel('[IP_3]');
ylabel('Amplitude [Ca^{2+}]');
%show period of last oscillation as dot. "(end)" - for last element in vector.
%plot(IP3,t(locs(end))-t(locs(end-1)),'mo');
%plot(Tinact,t(locs(end))-t(locs(end-1)),'mo');
 xlabel('[IP_3]');
%  xlabel('\tau_{inact}');
% ylabel('Period of [Ca^{2+} oscillations]');
% title('Saturation')
% title('Inactivation')
%title('Bifurcation diagram')
end
hold off

% Plot for Open probability
%ca2=0.00021

figure(1)
% hold on
%For P_open vs Ca (Task 1)
for IP3=6e-4:2e-4:2e-3
    ca2=[1e-5:1e-5:1e-2];

   fActSteady=(ca2.^Xact)./(ca2.^Xact+Kact.^Xact);
   fInactSteady=(ca2.^Xinact)./(ca2.^Xinact+Kinact^Xinact);
   pOpen=(fActSteady.*(1-fInactSteady))*(IP3^A)./(IP3^A+Kip3^A)*K;
semilogx(ca2,pOpen);
hold on
xlabel('[Ca^{2+}]');
ylabel('P_{open}');
end

figure(2)
 %For P_open vs IP3 (Task 1)
 IP3=[5e-4:2e-4:1e-2];
 ca2=0.0002;
   fActSteady=(ca2.^Xact)./(ca2.^Xact+Kact.^Xact);
   fInactSteady=(ca2.^Xinact)./(ca2.^Xinact+Kinact^Xinact);
   pOpen=(fActSteady.*(1-fInactSteady))*(IP3.^A)./(IP3.^A+Kip3^A)*K;
 plot(IP3,pOpen)
 xlabel('IP3');
 ylabel('P_{open}');
 title('Open probability')

    function dydt = model(t,y)
        dydt=zeros(2,1);  
        % diferential equations
        ca2=y(1);
        fInact=y(2);
        fActSteady=(ca2^Xact)/(ca2^Xact+Kact^Xact);
        fInactSteady=(ca2^Xinact)/(ca2^Xinact+Kinact^Xinact);
        pOpen=(fActSteady*(1-fInact))*(IP3^A)/(IP3^A+Kip3^A)*K;
                
        dydt(1)=Kca*(Rip3*pOpen+Rleak)-(Rserca*(ca2^P))/(ca2^P+Kserca^P);
        dydt(2)=(fInactSteady-fInact)/Tinact;
            
    end

end
