clc, clearvars

%SETUP PARA SIMULAÇÃO
Fs = 20000; %Frequência de amostragem
Ts = 1/Fs; %Período de amostragem
Vo_ref_amp = 120*sqrt(2); %Tensão de pico na  Vref
freqHzSin = 50; %Frequência da senoide de referência em Hz
freqSin = freqHzSin*2*pi; %frequência da senoide de referência em rad
Vdc = 400; %Tensão no barramento CC considerada pelo modelo
Vdc_sim = Vdc*1; %Tensão no barramento CC efetivamente aplicada na simulação
t_final = 3*(1/freqHzSin); %tempo final da simulação em amostras

s = tf('s'); %definir s para função de transferência
z = tf('z'); %definir z para função de transferência discreta

% P = 10; %ERRO DE MODELAGEM NOS COMPONENTES
P = 0; %ERRO DE MODELAGEM NOS COMPONENTES
% Vnoise = 2; %RUÍDO NA MEDIDA DA SAÍDA
Vnoise = 0; %RUÍDO NA MEDIDA DA SAÍDA

%define os horizontes
N = 4; %HORIZONTE DE PREDIÇÃO
Nu = N; %HORIZONTE DE CONTROLE
%--------
N1 = 2; %COMEÇO DO HORIZONTE DE PREDIÇÃO
N2 = N+N1-1; %FIM DO HORIZONTE DE PREDIÇÃO 
Nu1 = 1; %COMEÇO DO HORIZONTE DE CONTROLE
Nu2 = Nu+Nu1-1; %FIM DO HORIZONTE DE CONTROLE 

%COMPONENTES
L = 2e-3; %INDUTOR UTILIZADO NO MODELO
C = 50e-6; %CAPACITOR UTILIZADO NO MODELO
R = 60; %RESISTOR UTILIZADO NO MODELO
Lsim = L*(1-P/100); %INDUTOR UTILIZADO NA SIMULAÇÃO
Csim = C*(1-P/100); %CAPACITOR UTILIZADO NA SIMULAÇÃO
Rsim = R*(1-P/100); %RESISTOR UTILIZADO NA SIMULAÇÃO

%ABC 2 ALPHA BETA
Clarke = (2/3)*[1 -0.5 -0.5; 0 sqrt(3)/2 -sqrt(3)/2];
ABC_TF = (Vdc/(L*C))/(s^2 + (1/(R*C))*s + 1/(L*C));
G_xy = Clarke*ABC_TF;
Gd_xy = c2d(G_xy,Ts,'zoh');
% step(XY_TF);

num11 = Gd_xy(1,1).num{1};
num12 = Gd_xy(1,2).num{1};
num13 = Gd_xy(1,3).num{1};
num21 = Gd_xy(2,1).num{1};
num22 = Gd_xy(2,2).num{1};
num23 = Gd_xy(2,3).num{1};
den11 = Gd_xy(1,1).den{1};
den12 = Gd_xy(1,2).den{1};
den13 = Gd_xy(1,3).den{1};
den21 = Gd_xy(2,1).den{1};
den22 = Gd_xy(2,2).den{1};
den23 = Gd_xy(2,3).den{1};

A1 = conv(conv(den11, den12), den13);
A2 = conv(conv(den21, den22), den23);
if length(A2)<length(A1)
    A2 = [A2 zeros(1,length(A1)-length(A2))];
elseif length(A1)<length(A2)
    A1 = [A1 zeros(1,length(A2)-length(A1))];
end
A = [A1; A2];

B1 = deconv(conv(A1,num11),den11);
B1 = B1(2:end);
B2 = deconv(conv(A1,num12),den12);
B2 = B2(2:end);
B3 = deconv(conv(A1,num13),den13);
B3 = B3(2:end);
B4 = deconv(conv(A2,num21),den21);
B4 = B4(2:end);
B5 = deconv(conv(A2,num22),den22);
B5 = B5(2:end);
B6 = deconv(conv(A2,num23),den23);
B6 = B6(2:end);
B = [B1; B2; B3; B4; B5; B6];

dimA = length(A1);
dimB = length(B1);

vectors = [ 0   0   0;
            0   0   1;
            0   1   0;
            0   1   1;
            1   0   0;
            1   0   1;
            1   1   0;
            1   1   1];

%Diofantina
F1 = zeros(N2,dimB);
E1 = zeros(1,N2);
R_A = [1 zeros(1,dimA-1)] - A1;
E1(1) = 1;
F1(1,:) = R_A(2:end);
for i=2:N2
    E1(i) = F1(i-1,1);
    F1(i,:) = [F1(i-1,2:end) zeros(1,dimA-dimB)] - F1(i-1,1)*A1(2:end);
end

F2 = zeros(N2,dimB);
E2 = zeros(1,N2);
R_A = [1 zeros(1,dimA-1)] - A2;
E2(1) = 1;
F2(1,:) = R_A(2:end);
for i=2:N2
    E2(i) = F2(i-1,1);
    F2(i,:) = [F2(i-1,2:end) zeros(1,dimA-dimB)] - F2(i-1,1)*A2(2:end);
end

%Encontra resposta livre e forçada
G1 = zeros(N2,Nu);
G2 = zeros(N2,Nu);
G3 = zeros(N2,Nu);
G4 = zeros(N2,Nu);
G5 = zeros(N2,Nu);
G6 = zeros(N2,Nu);
f_aux1 = zeros(N2,(dimA-1)+(3*(dimB-1)));
f_aux2 = zeros(N2,(dimA-1)+(3*(dimB-1)));
for i=1:N2
    EB1 = conv(E1(1:i),B1);
    EB2 = conv(E1(1:i),B2);
    EB3 = conv(E1(1:i),B3);
    G1(i,i:-1:1) = EB1(1:i);
    G2(i,i:-1:1) = EB2(1:i);
    G3(i,i:-1:1) = EB3(1:i);
    f_aux1(i,:) = [F1(i,:) EB1(end-(dimB-2):end) EB2(end-(dimB-2):end) EB3(end-(dimB-2):end)];

    EB4 = conv(E1(1:i),B4);
    EB5 = conv(E1(1:i),B5);
    EB6 = conv(E1(1:i),B6);
    G4(i,i:-1:1) = EB4(1:i);
    G5(i,i:-1:1) = EB5(1:i);
    G6(i,i:-1:1) = EB6(1:i);
    f_aux2(i,:) = [F2(i,:) EB4(end-(dimB-2):end) EB5(end-(dimB-2):end) EB6(end-(dimB-2):end)];
end