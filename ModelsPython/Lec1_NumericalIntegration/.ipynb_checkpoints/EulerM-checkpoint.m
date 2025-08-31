%Created by Olga Sosnovtseva 10.08.2022
%Exercise on Euler method

function EulerM
y0 = 1; % initial condition
t0 = 0; % initial time
h = 1; % try: h = 0.01
tn = 4; % equal to: t0 + h*n, with n the number of steps
[t, y] = Eu(t0, y0, h, tn);
plot(t, y, 'b');
hold on

% exact solution (y = e^t):
    tt = t0:0.001:tn;
    yy = exp(tt);
    plot(tt, yy, 'r');
    legend('Euler', 'Exact');

% Euler method
function [t, y] = Eu(t0, y0, h, tn)
    fprintf('%10s%10s%10s%15s\n', 'i', 'ti', 'yi', 'f(yi,ti)');
    fprintf('%10d%+10.2f%+10.2f%+15.2f\n', 0, t0, y0, f(y0,t0));
    t = t0:h:tn;
    y = zeros(size(t));
    y(1) = y0;
    for i = 1:1:length(t)-1
        y(i+1) = y(i) + h*f(y(i),t(i));
        fprintf('%10d%+10.2f%+10.2f%+15.2f\n', i, t(i+1), y(i+1), f(y(i+1),t(i+1)));
    end
end

% Differential equation
function dydt = f(y,~)
    dydt = y;
end
end
