%% Solving the nonlinear Schrodinger equation using Fourier spectral method
% differentiation matrix 
clear all;close all;clc
% defining the computational parameters
N =1024; L =40*pi;
h = 2*pi/N; x = h*(1:N);x = L*(x-pi)/pi; 
dt = 0.001;
%% Second order differentiation
column = [-pi^2/(3*h^2)-1/6 ...
-.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2 = (pi/L)^2*toeplitz(column);
 

Tmax=50; tplot =0.075;

plotgap= round(tplot/dt); dt = tplot/plotgap;
nplots = round(Tmax/tplot);
% defining the initial condition
%% The fundamental soliton and its instability
%a = 10^(-4);
%u = sech(x).^2 + a*cos(x); uold = sech(x-dt).^2 + a*cos(x -dt);
u = exp(-(1/2)*x.^2); uold = exp( - (1/2)*x.^2- dt);
%% Twisted modes
%u = 2*sech(2.*x).^2 + 2i*sech(2.*x).*tanh(2.*x)  + a*cos(x);
%uold  = 2*sech(2.*(x - dt)).^2 + 2i*sech(2.*(x - dt)).*tanh(2.*(x - dt)) + a*cos(x - dt) ;
udata = [u; zeros(nplots,N)]; t=0; tdata=0;
%% loop
%coeff_fourier = 2*pi/(2*L);
b =0; lambda = 0;
% timing the loop
t0 = cputime;
for k =1:nplots
    for n = 1:plotgap
        t = t+dt;
        v_tilda = D2*u';
        %unew = uold -2*dt*sin(t).*(1./(1+k*u.^2)).*coeff_fourier.*v_tilda';
        unew  = uold  + (2*dt/(1i + lambda)).*( v_tilda' + 6*u.^2 - 4*u + b.*(u - conj(u)));
        uold =u; u=unew;
    end
    udata(k+1,:)=u;tdata =[tdata;t];
end
runningtime = cputime-t0 

%% plotting the numerical results
figure(1)
surf(x, tdata, imag(udata), 'FaceColor','interp','EdgeColor','none');view( -30,60) 
%colormap cool
%set(gca, 'FontSize',15)
set(gca, 'YDir','reverse')
axis([-5 5 0 Tmax min(min(imag(udata))) max(max(imag(udata)))])

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex'); 
zlabel('$Im\: \psi(x,t)$', 'Interpreter', 'latex');

figure(2)
surf(x, tdata, real(udata), 'FaceColor','interp','EdgeColor','none');view( -30,60) 
%colormap cool
%set(gca, 'FontSize',15)
set(gca, 'YDir','reverse')
axis([-5 5 0 Tmax min(min(real(udata))) max(max(real(udata)))])
%axis([-5 5 0 Tmax -3 3])
%colormap winter


xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex'); 
zlabel('$Re\: \psi(x,t)$', 'Interpreter', 'latex');
%colormap winter

figure(3)
surf(x, tdata, abs(udata), 'FaceColor','interp','EdgeColor','none');view( -30,60) 
%colormap cool
%set(gca, 'FontSize',15)
set(gca, 'YDir','reverse')
axis([-5 5 0 Tmax min(min(abs(udata))) max(max(abs(udata)))])

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex'); 
zlabel('$ |\psi(x,t)|$', 'Interpreter', 'latex');

figure(4)
plot(x, abs(udata(1, :)))
hold on 

plot(x, abs(udata(196, :)), 'g')

plot(x, abs(udata(201, :)))

legend('$|\psi(x, t=0)|$','$|\psi(x, t\approx14.6)|$', '$|\psi(x, t=15)|$', 'Interpreter', 'latex')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$|\psi(x, t)|$', 'Interpreter', 'latex')
axis([-3*pi 3*pi -0.3 max(max(abs(udata)))])
%imagesc(udata)