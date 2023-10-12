%% Question 1B: - 0.5  

% Passive membrane equation

clearvars;
close all

% Biophysical parameters

C = 1;
El = -65;
Gl = 0.1;
Iapp = -0.5;
D = 0;

% Time definitions

Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;

% Square wave (Heaviside function)

ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initial conditions (for D = 0)

V = zeros(1,length(t));
V(1) = El;

% Computation of the solution (for D = 0)

for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) = V(1);

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'b','linewidth',0.7);
set(gca,'fontsize',19);
hold on

D = 0;
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) = V(1);

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'r','linewidth',0.5);
axis([0 Tmax -65 -20]);
set(gca,'fontsize',19);
xlabel('t')
ylabel('V');
legend('analytical');
hold off; 
%% Question 1B: 0.5  

% Passive membrane equation

clearvars;
close all

% Biophysical parameters

C = 1;
El = -65;
Gl = 0.1;
Iapp = 0.5;
D = 0;

% Time definitions

Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;

% Square wave (Heaviside function)

ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initial conditions (for D = 0)

V = zeros(1,length(t));
V(1) = 3.5527e-14;

% Computation of the solution (for D = 0)

for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) = V(1);

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'b','linewidth',0.7);
set(gca,'fontsize',19);
legend('numerical');
hold on

D = 0;
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) =  -60.0000;

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'r','linewidth',0.5);
axis([0 Tmax -65 -20]);
set(gca,'fontsize',19);
xlabel('t')
ylabel('V');
legend('numerical','analytical');
hold off; 

%% Question 1B: 0.5  

% Passive membrane equation

clearvars;
close all

% Biophysical parameters

C = 1;
El = -65;
Gl = 0.1;
Iapp = 0.5;
D = 0;

% Time definitions

Tmax = 1000;
dt = 0.1;
t = 0:dt:Tmax;

% Square wave (Heaviside function)

ti = dt;
tf = 10000;
H = zeros(1,length(t));
H(floor(ti/dt):floor(tf/dt))=1;

% Initial conditions (for D = 0)

V = zeros(1,length(t));
V(1) = 3.5527e-14;

% Computation of the solution (for D = 0)

for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) = V(1);

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'b','linewidth',0.7);
set(gca,'fontsize',19);
legend('numerical');
hold on

D = 0;
for j=1:length(t)-1
    kv1 = (-Gl*(V(j)-El)+Iapp*H(j))/C;
    av = V(j)+kv1*dt;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    V(j+1) = V(j) + (kv1+kv2)*dt/2;   
end

% Initial conditions (for D > 0)

Vn = zeros(1,length(t));
Vn(1) =  -60.0000;

% Computation of the solution (for D > 0)

for j=1:length(t)-1
    eta = randn;
    Vnaux = Vn(j)+sqrt(2*D*dt)*eta;
    kv1 = (-Gl*(Vnaux-El)+Iapp*H(j))/C;
    av = Vn(j)+kv1*dt;
    av = av + sqrt(2*D*dt)*eta;
    kv2 = (-Gl*(av-El)+Iapp*H(j+1))/C;
    Vn(j+1) = Vn(j) + (kv1+kv2)*dt/2;
    Vn(j+1) = Vn(j+1) + sqrt(2*D*dt)*eta;   
end
% Graph
plot(t,V,'r','linewidth',0.5);
axis([0 Tmax -65 -20]);
set(gca,'fontsize',19);
xlabel('t')
ylabel('V');
legend('numerical','analytical');
hold off; 

