clc
clear
close all

syms theta2 theta2_0 L2 L3 K1 K2 K3 theta3_0 real


theta3 = asin((L2 / L3) * sin(theta2));


F_in = ...
    K1 * (theta2 - theta2_0) * ...
        (-cos(asin((L2 / L3) * sin(theta2))) / (L2 * sin(theta3 - theta2))) + ...
    K2 * (((theta2 - theta2_0) - (theta3 - theta3_0)) * ...
        (-cos(theta3) / (L2 * sin(theta3 - theta2))) + ...
        cos(theta2) / (L3 * sin(theta3 - theta2))) + ...
    K3 * ((theta3 - theta3_0) * -cos(theta2) / (L3 * sin(theta3)));

F_in = subs(F_in,theta3,asin((L2 / L3) * sin(theta2)));

disp(F_in);

%Real values
E = 8*10^9; %Youngs modulus of PLA 
l = 0.05; %length of compliant part
w = 0.05; %width of the compliant part
t = 0.005; %Thickeness
L2_real = 0.2;
L3_real = 0.25;
K1_real = (E*w*t^3)/(12*l);
K2_real = K1_real;
K3_real = K1_real;

theta2_0_real = pi/10;
theta3_0_real = 2*pi - pi/10;

F_in_real = subs(F_in,[K1 K2 K3 L2 L3 theta2_0 theta3_0], [K1_real K2_real K3_real L2_real L3_real theta2_0_real theta3_0_real]);

%Plotting

n = 101;
theta2_real = linspace(theta2_0_real,pi/2,n);
F_in_plot = zeros(n,1);

for i = 1:length(theta2_real)

    F_in_plot(i) = double(subs(F_in_real,theta2,theta2_real(i)));

end
theta2_real = theta2_real*180/pi;
F_in_plot = F_in_plot/1000;
figure(1)
plot(theta2_real,F_in_plot);
grid on
xlabel(['\theta_{2} [', char(176), ']'])
ylabel('Input force required [KN]')
xlim([0 100]);

W = trapz(theta2_real,F_in_plot/180*pi)*L2_real*2;
disp([num2str(W), ' KJ is the potential energy of the jumping robot'])
