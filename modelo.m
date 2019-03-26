function rRMSE=modelo

clear; clc;
exp_O2=csvread('exp_O2.csv');
exp_P =csvread('exp_P.csv');
exp_Xt=csvread('exp_Xt.csv');
exp_S=csvread('exp_S.csv');
exp_N=csvread('exp_N.csv');
exp_QO2X=csvread('exp_QO2X.csv');
save('Variaveis','exp_O2','exp_P','exp_Xt','exp_S','exp_N','exp_QO2X')

%SensGlobal
% [t,Xr,P,C]=modelo(mu_Xrmax, Yxs, Ks, Yxn, Kn, Yps, mu_Pmax, Kp,Yxo, Ypo, mo, Csat, theta)

%OtimParametros
% SSE=modelo(v)
%  Kpn=v(1);Kps=v(2);Kxn=v(3);Kxs=v(4);Ypo=v(5);Yps=v(6);Yxn=v(7);Yxo=v(8);Yxs=v(9);mu_Pmax=v(10);mu_Xrmax=v(11);theta=v(12);

Salvar=1; Otimizacao=1;
Xri=0.35; Pi=0.21;    Si=19.4; Ni=2.02;
%%%%%%%%%%%%%%%%%%%%%%%% Cinética microbiana

Csat=0.00571;   % [gO2/L] Concentração de saturação
Kpn=0.424;       % [g/L] Constante de afinidade glicerol-produto
Kps=0.285;      % [g/L] Constante de afinidade glicerol-produto
Kxn=0.075;      % [g/L] Constante de afinidade
Kxs=3.59270;           % [g/L] Constante de afinidade
Ypo=15;%2.65;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
Yps=0.4383;       % [gP/gS] Fator de conv. de glicerol em produto
Yxn=8.1707;        % [gX/gN] Fator de conv. de nitrogênio em biomassa
Yxo=1.825;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
Yxs=0.46;       % [gX/gS] Fator de conv. de substrato em biomassa
mo=0;  % [gO2/gXr/h] Consumo de O2 para manutenção celular.
ms=0;     % [gS/gX/h] Consumo de S para manutenção celular
mu_Pmax=0.135;%106;  % [gP/gXr/h] Velocidade máxima de produção
mu_Xrmax=0.1710; % [h^-1] Vel. específica de crescimento
theta=0.4609;      % [-] Expoente fração acumulada


PXmax=4;        % [-] Capacidade de acúmulo máxima

Ci=Csat*0.90;   % [gO2/L] Concentração inicial

% ==== AGITAÇÃO E MISTURA ====
V=4;            % [L] Volume inicial de meio
Nrot=450;       % [rpm] Velocidade de agitação
Q=0.2;          % [lpm] Vazão de aeração(1 atm @ 20°C)

Fs=0; Sm=670.54;
Fn=0; Nm=10.61;

% 


% 

% 


ncond=18;
tf=[15.1 20.5 22.15 22.9 26.1 28.3 29.4 30.7 39.2 40.66 42.5 46 55 58.8 62 63.667 75 85];


tspan=linspace(0,tf(1),tf(1)*4+1);
par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Kxs Kxn Kps Kpn PXmax theta mo ms Csat Nrot Q Fs Sm Fn Nm];
[t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
EDO=[t,y];
kLa(1,1)=0; kLa(1,2)=agitacao(Nrot,Q,V);
kLa(2,1)=tf(1); kLa(2,2)=kLa(1,2);

for j=2:ncond
    Xri=y(end,1);Si=y(end,2);Ni=y(end,3);Pi=y(end,4);Ci=y(end,5);
    switch j
        case 2 % Casos para j=t-1
            Nrot=525;
        case 3
            Nrot=600;
        case 4
            Si=24.2;
        case 5
            Nrot=700;
        case 6
            Nrot=900; Q=1;
        case 7
            Q=3;
        case 8
            Q=5;Si=22;
        case 9
            Q=4;Fs=16*10^-3;Fn=Fs;
        case 10
            Q=2;Fs=18*10^-3;Fn=Fs;
        case 11
            Ni=0.11;Si=Si+5.1;
        case 12
            Ni=0.11;Q=3;Si=Si+6; 
        case 13
            Ni=0.15;Fs=35.2*10^-3; Fn=Fs; Si=Si+10;
        case 14
            Fn=0;Fs=Fn;
        case 15
            Q=1;
        case 16
            Q=0.8; 
        case 17
            Fs=42.2*10^-3; Fn=Fs; Nm=16.61;
    end      
 
       
tspan=linspace(tf(j-1)+0.01,tf(j),(tf(j)-tf(j-1))*4+1);
par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Kxs Kxn Kps Kpn PXmax theta mo ms Csat Nrot Q Fs Sm Fn Nm];
[t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
EDO=[EDO;t,y];
kLa(j*2-1,1)=tf(j-1); kLa(j*2-1,2)=agitacao(Nrot,Q,V);
kLa(j*2,1)=tf(j); kLa(j*2,2)=agitacao(Nrot,Q,V);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===== CÁLCULOS DENTRO DA EDO =====
P=zeros(1,length(EDO)); muP=zeros(1,length(EDO));
 S=zeros(1,length(EDO)); muS=zeros(1,length(EDO));
 N=zeros(1,length(EDO)); QO2=zeros(1,length(EDO));
 C=zeros(1,length(EDO));muXr=zeros(1,length(EDO));
Xt=zeros(1,length(EDO));  Xr=zeros(1,length(EDO));
QO2X=zeros(1,length(EDO));

for i=1:length(EDO)
    t(i)=EDO(i,1); Xr(i)=EDO(i,2); S(i)=EDO(i,3); 
    N(i)=EDO(i,4); P(i)=EDO(i,5); C(i)=EDO(i,6);
    
    muXr(i)=mu_Xrmax*(S(i)/(Kxs+S(i)))*(N(i)/(Kxn+N(i)));
    
    muP(i)=mu_Pmax*(1-N(i)/(N(i)+Kpn))*(S(i)/(S(i)+Kps))*(1-((P(i)/Xr(i))/PXmax)^theta);

    muS(i)=muXr(i)/Yxs+muP(i)/Yps+ms;
    
    QO2(i)=muXr(i)/Yxo+muP(i)/Ypo+mo;
    
    Xt(i)=Xr(i)+P(i);
    QO2X(i)=QO2(i)*Xr(i);
end
%%%%%%%%%%%%%%%%%%%%%% OTIMIZAÇÃO DE PARÂMETROS
if Otimizacao>0
    n_func_obj=6;
    C=C(:); P=P(:); Xt=Xt(:); S=S(:); N=N(:); QO2X=QO2X(:); 
    SSE=zeros(1,n_func_obj); Ef=SSE; Ef_den=SSE;rRMSE=SSE;
    dados_modelo=unique([t,C,P,Xt,S,N,QO2X],'rows');
    dados_modelo=unique(dados_modelo,'rows');
    load('Variaveis');
    
    w=[length(exp_O2) length(exp_P) length(exp_Xt) length(exp_S) length(exp_N) length(exp_QO2X)];
    for j=1:n_func_obj
        switch j
            case 1 
                dados_exp=exp_O2;
            case 2 
                dados_exp=exp_P;
            case 3
                dados_exp=exp_Xt;
            case 4
                dados_exp=exp_S;
            case 5
                dados_exp=exp_N;
            case 6
                dados_exp=exp_QO2X;
        end
        modelo_interp=zeros(1,w(j)); 
        media_dados_exp=mean(dados_exp(:,2));
        for i=1:w(j)
            modelo_interp(i)=interp1(dados_modelo(:,1),dados_modelo(:,j+1),dados_exp(i,1));
            SSE(j)=SSE(j)+(modelo_interp(i)-dados_exp(i,2))^2; % Função multiobjetivo
            Ef_den(j)=Ef_den(j)+(modelo_interp(i)-media_dados_exp)^2;
        end
        Ef(j)=1-(SSE(j))/(Ef_den(j));
        rRMSE(j)=((SSE(j)/w(j))^0.5)/media_dados_exp;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ef
% ================ Resultados =================
if Salvar>0 % Salvar dados em planilha xls
    muXr=muXr(:);muS=muS(:);muP=muP(:);Xt=Xt(:);QO2=QO2(:);
    Dados=[EDO Xt muXr muS muP QO2];
    xlswrite('Modelo.xlsx',Dados,1,'A6');
    xlswrite('Modelo.xlsx',kLa,1,'O6');
    xlswrite('Modelo.xlsx',par,1,'B3');
    winopen('Modelo.xlsx')
end

end

function dy=balancos(t,y,par)
 persistent VDados
% ==== Redefinição dos parâmetros ==== 
Yxs=par(1); Yxn=par(2); Yps=par(3); Yxo=par(4); Ypo=par(5); 
mu_Xrmax=par(6); mu_Pmax=par(7); 
Kxs=par(8); Kxn=par(9); Kps=par(10); Kpn=par(11); PXmax=par(12); theta=par(13);
mo=par(14); ms=par(15); Csat=par(16); Nrot=par(17); Q=par(18); 
Fs=par(19); Sm=par(20); Fn=par(21); Nm=par(22);

dy=zeros(5,1); Xr=y(1); S=y(2); N=y(3); P=y(4); C=y(5);
% =========================

muXr=mu_Xrmax*(S/(Kxs+S))*(N/(Kxn+N));
muP=mu_Pmax*(S/(S+Kps))*(1-N/(N+Kpn))*(1-((P/Xr)/PXmax)^theta);
muS=muXr/Yxs + muP/Yps+ms;
QO2=muXr/Yxo + muP/Ypo+mo;

if t<0.5 VDados=csvread('exp_Volume.csv');end
dvdt = interp1(VDados(:,1),VDados(:,2),t);
   V = interp1(VDados(:,1),VDados(:,3),t);

kLa=agitacao(Nrot,Q,V);

dy(1)=  Xr*muXr - dvdt/V*Xr;                 % Biomassa residual
dy(2)= -Xr*muS - dvdt/V*S + Fs*Sm/V;                   % Glicerol
dy(3)= -Xr*muXr/Yxn - dvdt/V*N + Fn*Nm/V;              % Nitrogênio
dy(4)=  Xr*muP - dvdt/V*P;                   % PHB
dy(5)=  kLa*(Csat-C) - QO2*Xr - dvdt*C;      % Oxigênio dissolvido
end

function kLa=agitacao(Nrot,Q,V)
% Geometria do biorreator
Dt=0.18;  % [m] Diâmetro do tanque
Di=0.059; % [m] Diâmetro do agitador
Wi=0.018; % [m] Largura do agitador

% Geral
g=9.81;   % [m/s^2] Gravidade
Np=5.5;   % [-] Número de potência
rho=1000; % [kg/m^3] Densidade

% Conversões
Nrot=Nrot/60;      % [min^-1] to [s^-1]
Q=Q*303.15/298.15; % Correção do volume pela mudança de T
Q=Q/1000/60;       % [lpm] to [m^3/s]
Vliq=V/1000;       % [L]   to [m^3]

% =============================================

% =========== Transferência de pot. ===========
% Como deltaC=Di, P2=1.5*P1 (Hudcova, 1989)
P=Np*rho*(Nrot^3)*(Di^5)*1.5; % [W]

As=(pi*Dt^2)/4; % [m^2] Área superficia do líquido
vs=Q/As;        % [m/s] Velocidade superficial

% ==== Pot. transferida sob aeração
Pg=1.224*(P^2*Nrot*Di^3/(Q^0.56))^0.432; % [W] Abardi (1988)
% ==== 
%Pg=P*0.10*(Q/(Nrot*Vliq))^(-1/4)*(Nrot^2*Di^4/(g*Wi*Vliq^(2/3)))^(-1/5); % [W] Hughmark (1980)

% ============================================

% =================== kLa ====================
% ==== Van't Riet (1979)
%c=0.002; alpha=0.7; beta=0.36;
% ==== Meu
c=0.012; alpha=0.308; beta=0.165;

kLa=c*(Pg/Vliq)^alpha*(vs)^beta*3600; % [h^-1]
end








% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CULTIVO D
% Xri=0.25; Si=24.17; Ni=1.23;
% ncond=5;
% tf=[12 18.6 23 27 48];
% 
%     switch j
%         case 2 % Caso para t=j-1
%             Nrot=600;
%         case 3
%             Q=0.4;
%         case 4
%             Nrot=750;
%         case 5
%             Si=30.02;
%     end     


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CULTIVO G
% Xri=0.50; Pi=0.05; Si=20.14; Ni=1.13;
% ncond=10;
% tf=[15.1 22.7 24 24.75 32.5 33 38.5 44 65.2 82];
% 
%     switch j
%         case 2 % Casos para t=j-1
%             Nrot=600;
%         case 3
%             Nrot=800;
%         case 4
%             Si=28.11;
%         case 5
%             Q=0.6;
%         case 6 
%             Si=35.07;
%         case 7
%             Q=0.4;
%         case 8
%             Si=41.35;
%         case 9
%             Q=0.2; Nrot=650;
%         case 10
%             Nrot=600;
%     end     

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CULTIVO H
% Xri=0.37; Pi=0.13;    Si=22.93; Ni=0.86;
% Fs=0; Sm=709;
% Fn=0; Nm=16;
% 
% ncond=17;
% tf=[10.3 16.7 21.5 23.7 25.4 30.21 33.1 35.5 38.81 42.8 45.5 54.18 55.51 60.6 68.95 70 72];
% 
% 
%     switch j
%         case 2 % Casos para j=t-1
%             Fs=3.58*10^-3;
%         case 3
%             Nrot=600;
%         case 4
%             Q=0.2; Fs=5.3*10^-3;
%         case 5
%             Fn=0;Fs=15.1*10^-3;
%         case 6
%             Nrot=750;
%         case 7
%             Ni=0.11; 
%         case 8
%             Ni=0.11; Fn=8*10^-3;
%         case 9
%             Fs=0;
%         case 10
%             Nrot=800;Fs=39.1*10^-3; Ni=0.11;
%         case 11
%             Nrot=850; Fs=7.9*10^-3;
%         case 12
%             Q=0.6; Fs=15.9*10^-3; Ni=0.11;
%         case 13
%             Q=1.8;Fs=26.9*10^-3;Ni=0.11;
%         case 14
%             Ni=0.11;
%         case 15
%             Q=2.2; Nrot=925; Fs=2.9*10^-3; Ni=0.11;
%         case 16
%             Fs=0;
%         case 17
%             Fn=0;
%     end 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CULTIVO I
% Xri=0.35; Pi=0.21;    Si=19.4; Ni=2.02;
% ncond=18;
% tf=[15.1 20.5 22.15 22.9 26.1 28.3 29.4 30.7 39.2 40.66 42.5 46 55 58.8 62 63.667 75 85];
% 
% Fs=0; Sm=670.54;
% Fn=0; Nm=10.61;
% 
%     switch j
%         case 2 % Casos para j=t-1
%             Nrot=525;
%         case 3
%             Nrot=600;
%         case 4
%             Si=24.2;
%         case 5
%             Nrot=700;
%         case 6
%             Nrot=900; Q=1;
%         case 7
%             Q=3;
%         case 8
%             Q=5;Si=22;
%         case 9
%             Q=4;Fs=16*10^-3;Fn=Fs;
%         case 10
%             Q=2;Fs=18*10^-3;Fn=Fs;
%         case 11
%             Ni=0.11;Si=Si+5.1;
%         case 12
%             Ni=0.11;Q=3;Si=Si+6; 
%         case 13
%             Ni=0.15;Fs=35.2*10^-3; Fn=Fs; Si=Si+10;
%         case 14
%             Fn=0;Fs=Fn;
%         case 15
%             Q=1;
%         case 16
%             Q=0.8; 
%         case 17
%             Fs=42.2*10^-3; Fn=Fs; Nm=16.61;
%     end   







% 
% mu_Xrmax=0.158; % [h^-1] Vel. específica de crescimento
% 
% % Glicerol = S
% Yxs=0.43;       % [gX/gS] Fator de conv. de substrato em biomassa
% Kxs=2.68;           % [g/L] Constante de afinidade
% ms=2*10^-4;     % [gS/gX/h] Consumo de S para manutenção celular
% % Fonte de nitrogênio = N
% Yxn=7.8;        % [gX/gN] Fator de conv. de nitrogênio em biomassa
% Kxn=0.189;      % [g/L] Constante de afinidade
% 
% % PHB = P
% Yps=0.42;       % [gP/gS] Fator de conv. de glicerol em produto
% mu_Pmax=0.097;  % [gP/gXr/h] Velocidade máxima de produção
% Kps=1.73;      % [g/L] Constante de afinidade glicerol-produto
% Kpn=0.04;       % [g/L] Constante de afinidade glicerol-produto
% PXmax=4;        % [-] Capacidade de acúmulo máxima
% theta=0.8;      % [-] Expoente fração acumulada
% 
% Pi=Xri*0.1;     % [g/L] Concentração inicial de PHB
% 
% %Oxigênio dissolvido = C
% Yxo=2.28;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
% Ypo=2.54;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
% mo=2.33*10^-4;  % [gO2/gXr/h] Consumo de O2 para manutenção celular.
% Csat=0.00571;   % [gO2/L] Concentração de saturação
% Ci=Csat*0.90;   % [gO2/L] Concentração inicial
