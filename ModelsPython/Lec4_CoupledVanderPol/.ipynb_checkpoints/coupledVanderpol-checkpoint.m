%Created by Olga Sosnovtseva 09.06.2022
%Coupled van der Pol oscillators
function coupledVanderpol
% Set of parameters 
lambda1=0.1;
lambda2=0.1;
w1=1.0;
w2=0.98; %for locking
%w2=0.85; %for suppression
B=0.0;

% set integration time, time linespace and initial conditions
tSpan = [0 2000];
%tTransient=0;
tTransient=1000; %B=0.11
pointsN=20000;
t = linspace(tTransient, tSpan(2), pointsN);
yZero = [-0.5;0.1;0.5;0.1];

% open new figure,  hold on to draw multiple images in one figure
figure (1)
%run solver choose ode45, ode23s or smth else; deval it in region of interest
    sol=ode45(@model, tSpan, yZero);
    y = deval(sol,t); 
    
%draw figures, plot (X,Y)
    plot(y(1,:),y(3,:));
%    xlabel('x1');
%    ylabel('x2');
%    xlim([1000,1100])
%    xlim([-0.1 0.1]);
%    ylim([-0.1 0.1]);
%    title('Time series, B=0.015')
  
figure(2)    
%FFT, fs is sampling frequency, nfft is number of points in time series
nfft = length(y(1,:));
    fs=pointsN/(tSpan(2)-tTransient);
    [pxx1,f1] = periodogram(y(1,:),[],nfft,fs);
    [pxx2,f2] = periodogram(y(3,:),[],nfft,fs);
    plot(f1,pxx1,f2,pxx2);
    xlim([0.12 0.2])
%    ylim([0 140])  %for suppression
    xlabel('f');
    ylabel('S(f)');

% model itself (diferential equations
function dydt = model(t,y)
dydt=zeros(4,1);   
  xx1=y(1);  
  yy1=y(2);
  xx2=y(3);  
  yy2=y(4);
    % diferential equations
    dydt(1)=yy1;
    dydt(2)=(lambda1-xx1^2)*yy1-w1^2*xx1+B*(yy2-yy1);
    dydt(3)=yy2;
    dydt(4)=(lambda2-xx2^2)*yy2-w2^2*xx2+B*(yy1-yy2);   
end

end
