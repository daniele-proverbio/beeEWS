%% Create other figures for the article

clc; clear; close all;

%% Fig 1: examples of dynamics

% No transition
x = linspace(1,100);
y = rand(1,length(x)) + 2;

figure()
plot(x,y,linewidth=3)
ylim([0,15])

% Smooth non-bifurcation
y2 = 100 + 0.04*rand(1,length(x)) + exp(-0.1*x);

figure()
plot(x,y2,linewidth=3)


% Abrupt non-bifurcation
y3 = zeros(1,length(x));
y3(1:50) = sqrt(100-x(1:50)) + 0.5*rand(1,50);
y3(51:end) = 5 + 0.5*rand(1,50);

plot(x,y3)

%% Transcritical
y4 = zeros(1,length(x));
y4(1:50) = 5 - 0.1*x(1:50) + 0.5*rand(1,50);
y4(51:end) = 0.3*rand(1,50);
y4(y4<0) = 0;

plot(x,y4)


%% Bifurcation plots for the models

% Parameters (values used in the respective publications)
L = 2000;
alfa_1 = 0.25;
alfa_2 = 0.25;
sigma = 0.75;
w = 27000;

v = 5000;
phi = 1/9;
gamma_B = 0.018;
gamma_A = 0.007;
b = 500;
c = 0.1;


% --- Khoury et al, 2011
m_c = L / (2*w*(alfa_1 - L/w)) * (alfa_1+sigma+ sqrt((alfa_1-sigma)^2 + 4*sigma*L/w) );
m = linspace(0.1,m_c,1000);
mm = linspace(m_c,0.55,1000); 

J = 0.5*((alfa_1./m - sigma./m - 1) + sqrt((alfa_1./m - sigma./m - 1).^2 + 4*alfa_1./m));
JJ = 0.5*((alfa_1./mm - sigma./mm - 1) + sqrt((alfa_1./mm - sigma./mm - 1).^2 + 4*alfa_1./mm));
FH = (L./m).*(1+J)./J - w; 
FHFH = (L./mm).*(1+JJ)./JJ - w; 

figure()
hold on
plot(m, FH, linewidth = 1.5, Color='black')
plot(mm, FHFH, '--', linewidth = 1.5, Color='black')
plot(m, 0*m, '--', linewidth = 1.5, Color='black')
plot(mm, 0*mm, linewidth = 1.5, Color='black')
xlabel('m', fontsize = 20)
ylabel('$\hat{H}+\hat{F}$', Interpreter='latex', fontsize = 20)
ylim([-10000,40000])


% --- Khoury et al., 2013
m_c2 = (c-gamma_A)*phi/gamma_B;
m2 = linspace(0.3,0.55,5000);

P = c/gamma_A - 1 - (gamma_B * m2)/(gamma_A * phi);
f = b*sqrt(alfa_2 ./ ((sigma./(P+1)) + (m2./P) - alfa_1) -1);
Fhat = (L*f.^2) ./ (m2 .* (f.^2 + b^2));

HBF = Fhat .* (P + m2/phi +1);
Hhat = Fhat .* P;
idx = find(HBF < 0, 1, 'first');

figure()
hold on
plot(m2(1:idx), HBF(1:idx), linewidth = 1.5, Color='black')
plot(m2(idx:end), HBF(idx:end), '--', linewidth = 1.5, Color='black')
plot(m2(1:idx), 0*m2(1:idx), '--', linewidth = 1.5, Color='black')
plot(m2(idx:end), 0*m2(idx:end), linewidth = 1.5, Color='black')
xlabel('m', fontsize = 20)
ylabel('$\hat{H}+\hat{F}+\hat{B}$', Interpreter='latex', fontsize = 20)
ylim([-10000,40000])



% --- Breda et al., 2021

n = [0.1 0.25 0.4 0.55 0.7];
L1 = (1-n)*L;

colors = [0., 0.1, 0.2 ; 0.25, 0.35, 0.45; 0.35, 0.45, 0.55; 0.55, 0.65, 0.75; 0.7, 0.8, 0.9]; 

figure()
hold on

for i =1:length(n)

    m_c = L1(i) / (2*w*(alfa_1 - L1(i)/w)) * (alfa_1+sigma+ sqrt((alfa_1-sigma)^2 + 4*sigma*L1(i)/w) );
    m = linspace(0.01,m_c,1000);
    mm = linspace(m_c,0.55,1000); 
    
    J = 0.5*((alfa_1./m - sigma./m - 1) + sqrt((alfa_1./m - sigma./m - 1).^2 + 4*alfa_1./m));
    JJ = 0.5*((alfa_1./mm - sigma./mm - 1) + sqrt((alfa_1./mm - sigma./mm - 1).^2 + 4*alfa_1./mm));
    FH = (L1(i)./m).*(1+J)./J - w; 
    FHFH = (L1(i)./mm).*(1+JJ)./JJ - w; 
    
    plot(m, 0*m, '--', linewidth = 1.5, Color='black')
    plot(mm, 0*mm, linewidth = 1.5, Color='black')
    plot(m, FH, linewidth = 1.5, Color=colors(i,:))
    plot(mm, FHFH, '--', linewidth = 1.5, Color=colors(i,:))

end

xlabel('m', fontsize = 20)
ylabel('$\hat{H}+\hat{F}$', Interpreter='latex', fontsize = 20)
ylim([-20000,40000])
legend('', '', 'n=0.1', '', '', '',  'n=0.25', '', '', '', 'n=0.4', '', '', '', 'n=0.55', '', '', '', 'n=0.7', fontsize = 15)