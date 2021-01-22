%% Solution of the quadratic nonlinear Schrodinger equation via Fourier spectral methods
% differentiation matrix 
clear all;close all;clc
% defining the computational parameters
N = 180; L = 40*pi; dt = 0.001;
% computing the first order differentiation matrix via Fourier Method and
% the x value

[x,D]=fourier_matrix(N);x=(20*pi*x-20*pi)', D2 = D*D;
Tmax=50; tplot =0.15;

plotgap= round(tplot/dt); dt = tplot/plotgap;
nplots = round(Tmax/tplot);

% defining the initial conditions
%right wave 
u = exp(-x.^2); uold = exp( - x.^2- dt);
%left wave
%v = sech(x + 1).^2; vold = sech(x +1  - dt).^2;
% Main calculation
udata = [u; zeros(nplots,N)];%vdata = [v; zeros(nplots,N)]; 
t=0; tdata=0;


%% loop
coeff_fourier = (pi/L)^2;
% timing the loop
t0 = cputime;
%constants 
b = 0 ;lambda=0;

%% solving for the right wave
for k =1:nplots
    for n = 1:plotgap
        t = t+dt;
        v_tilda = D2*u';
        %unew = uold -2*dt*sin(t).*(1./(1+k*u.^2)).*coeff_fourier.*v_tilda';
        unew = uold + (2*dt/(1i + lambda )).*(coeff_fourier.*v_tilda' + 6*u.^2 - 4*u + b.*(u - conj(u)));
        uold =u; u=unew;
    end
    udata(k+1,:)=u; tdata =[tdata;t];
end

%% solving for the left wave

%for k =1:nplots
    %for n = 1:plotgap
        %t = t+dt;
        %w_tilda = D2*v';
        %unew = uold -2*dt*sin(t).*(1./(1+k*u.^2)).*coeff_fourier.*v_tilda';
        %vnew = vold + (2*dt/(1i + lambda )).*(coeff_fourier.*w_tilda' + 6*v.^2 - 4*v + b.*(v - conj(v)));
        %vold =v; v=vnew;
    %end
    %vdata(k+1,:)=v; 
%end

runningtime = cputime-t0 


%% plotting the numerical results

%surf(x,tdata,udata)

figure(1)
surf(x, tdata, imag(udata), 'FaceColor','interp','EdgeColor','none');view( -30,60) 
set(gca,'FontSize',15)
axis([-2*pi 2*pi 0 Tmax min(min(imag(udata))) max(max(imag(udata)))])

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex'); 
zlabel('$\psi(x,t)$', 'Interpreter', 'latex');

figure(2)
surf(x, tdata, real(udata), 'FaceColor','interp','EdgeColor','none');view( -30,60) 
set(gca,'FontSize',15)
axis([-2*pi 2*pi 0 Tmax min(min(real(udata))) max(max(real(udata)))])

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$t$', 'Interpreter', 'latex'); 
zlabel('$\psi(x,t)$', 'Interpreter', 'latex');

