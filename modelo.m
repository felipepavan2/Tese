function [t,Xr,P,C]=modelo(Csat, Kpn, Kps, Kxn, Kxs, Ypo, Yps, Yxn, Yxo, Yxs, mo, ms, mu_Pmax, mu_Xrmax, theta)
% 

%  clc; clear;

%SensGlobal
% [t,Xr,P,C]=modelo(mu_Xrmax, Yxs, Ks, Yxn, Kn, Yps, mu_Pmax, Kp,Yxo, Ypo, mo, Csat, theta)

%OtimParametros
% SSE=modelo(v)

% Kxn=x(1);Kps=x(2);Kxs=x(3);
% Ypo=x(4);Yxn=x(5);Yxo=x(6);Yxs=x(7);
% mu_Xrmax=x(8); mu_Pmax=x(9);theta=x(10);
% antes SSE  C=0.0008   P=43.6410   Xt=44.2664

Xri=0.32; Si=32; Ni=1;
Salvar=0; Otimizacao=0;

%%%%%%%%%%%%%%%%%%%%%%%% Cinética microbiana

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
% Kps=3.496;      % [g/L] Constante de afinidade glicerol-produto
% Kpn=0.04;       % [g/L] Constante de afinidade glicerol-produto
PXmax=4;        % [-] Capacidade de acúmulo máxima
% theta=0.8;      % [-] Expoente fração acumulada

Pi=Xri*0.1;     % [g/L] Concentração inicial de PHB

%Oxigênio dissolvido = C
% Yxo=1.93;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
% Ypo=3.34;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
% mo=2.33*10^-4;  % [gO2/gXr/h] Consumo de O2 para manutenção celular.
% Csat=0.00571;   % [gO2/L] Concentração de saturação
Ci=Csat*0.90;   % [gO2/L] Concentração inicial

% ==== AGITAÇÃO E MISTURA ====
V=4;            % [L] Volume inicial de meio
Nrot=450;       % [rpm] Velocidade de agitação
Q=0.2;          % [lpm] Vazão de aeração(1 atm @ 20°C)



ncond=1;
tf=[60];

tspan=linspace(0,tf(1),tf(1)*4+1);
par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Kxs Kxn Kps Kpn PXmax theta mo ms Csat Nrot Q];
[t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
EDO=[t,y];
kLa(1,1)=0; kLa(1,2)=agitacao(Nrot,Q,V);
kLa(2,1)=tf(1); kLa(2,2)=kLa(1,2);

% for j=2:ncond
%     Xri=y(end,1);Si=y(end,2);Ni=y(end,3);Pi=y(end,4);Ci=y(end,5);
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
% tspan=linspace(tf(j-1),tf(j),(tf(j)-tf(j-1))*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Kxs Kxn Kps Kpn PXmax theta mo ms Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[EDO;t,y];
% kLa(j*2-1,1)=tf(j-1); kLa(j*2-1,2)=agitacao(Nrot,Q,V);
% kLa(j*2,1)=tf(j); kLa(j*2,2)=agitacao(Nrot,Q,V);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ===== CÁLCULOS DENTRO DA EDO =====
P=zeros(1,length(EDO)); muP=zeros(1,length(EDO));
 S=zeros(1,length(EDO)); muS=zeros(1,length(EDO));
 N=zeros(1,length(EDO)); QO2=zeros(1,length(EDO));
 C=zeros(1,length(EDO));muXr=zeros(1,length(EDO));
Xt=zeros(1,length(EDO));  Xr=zeros(1,length(EDO));

for i=1:length(EDO)
    t(i)=EDO(i,1); Xr(i)=EDO(i,2); S(i)=EDO(i,3); 
    N(i)=EDO(i,4); P(i)=EDO(i,5); C(i)=EDO(i,6);
    
    muXr(i)=mu_Xrmax*(S(i)/(Kxs+S(i)))*(N(i)/(Kxn+N(i)));
    
    muP(i)=mu_Pmax*(1-N(i)/(N(i)+Kpn))*(S(i)/(S(i)+Kps))*(1-((P(i)/Xr(i))/PXmax)^theta);

    muS(i)=muXr(i)/Yxs+muP(i)/Yps+ms;
    
    QO2(i)=muXr(i)/Yxo+muP(i)/Ypo+mo;
    
    Xt(i)=Xr(i)+P(i);
end

Xr=Xr(:);P=P(:);C=C(:);S=S(:);N=N(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OTIMIZAÇÃO DE PARÂMETROS
if Otimizacao>0
    C=C(:); P=P(:); Xt=Xt(:); SSE=[0 0 0];
    dados_modelo=unique([t,C,P,Xt],'rows');
    exp_O2=csvread('exp_O2.csv');
    exp_P =csvread('exp_P.csv');
    exp_Xt=csvread('exp_Xt.csv');
% load('Variaveis','exp_O2','exp_P','exp_Xt');
    w=[length(exp_O2) length(exp_P) length(exp_Xt)];
    for j=1:3
        switch j
            case 1
                dados_exp=exp_O2;
                modelo_interp=zeros(1,w(j));
            case 2
                dados_exp=exp_P;
                modelo_interp=zeros(1,w(j));
            case 3
                dados_exp=exp_Xt;
                modelo_interp=zeros(1,w(j));
        end
        for i=1:w(j)
            modelo_interp(i)=interp1(dados_modelo(:,1),dados_modelo(:,j+1),dados_exp(i,1));
            SSE(j)=SSE(j)+(modelo_interp(i)-dados_exp(i,2))^2; % Função multiobjetivo
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================ Resultados =================
if Salvar>0 % Salvar dados em planilha xls
    muXr=muXr(:);muS=muS(:);muP=muP(:);Xt=Xt(:);QO2=QO2(:);
    Dados=[EDO Xt muXr muS muP QO2];
    xlswrite('Modelo.xlsx',Dados,1,'A4');
    xlswrite('Modelo.xlsx',kLa,1,'O4');
end
% disp('Fim.');
% disp(SSE);

end

function dy=balancos(t,y,par)
% persistent VDados
% ==== Redefinição dos parâmetros ==== 
Yxs=par(1); Yxn=par(2); Yps=par(3); Yxo=par(4); Ypo=par(5); 
mu_Xrmax=par(6); mu_Pmax=par(7); 
Kxs=par(8); Kxn=par(9); Kps=par(10); Kpn=par(11); PXmax=par(12); theta=par(13);
mo=par(14); ms=par(15); Csat=par(16); Nrot=par(17); Q=par(18);

dy=zeros(5,1); Xr=y(1); S=y(2); N=y(3); P=y(4); C=y(5);
% =========================

muXr=mu_Xrmax*(S/(Kxs+S))*(N/(Kxn+N));
muP=mu_Pmax*(S/(S+Kps))*(1-N/(N+Kpn))*(1-((P/Xr)/PXmax)^theta);
muS=muXr/Yxs + muP/Yps+ms;
QO2=muXr/Yxo + muP/Ypo+mo;

%if t<0.5 VDados=csvread('exp_Volume.csv');end
dvdt = 0;%interp1(VDados(:,1),VDados(:,2),t);
   V = 4;%interp1(VDados(:,1),VDados(:,3),t);


% if t<1
%     dvdt=0;V=Volume(1,3);
% elseif t==82
%     dvdt=0;[a,b]=size(Volume);
%     V=Volume(a,1);
% else
% Vpos=find(abs(Volume(:,1)-t) < 0.05);
% dvdt=Volume(Vpos,2);dvdt(1,:)=[];
% V=Volume(Vpos,3);V(1,:)=[];
% end

kLa=agitacao(Nrot,Q,V);

dy(1)=  Xr*muXr - dvdt/V*Xr;                 % Biomassa residual
dy(2)= -Xr*muS - dvdt/V*S;                   % Glicerol
dy(3)= -Xr*muXr/Yxn - dvdt/V*N;              % Nitrogênio
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
% Salvar=1; Otimizacao=0;
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Cinética microbiana
% 
% mu_Xrmax=0.156; % [h^-1] Vel. específica de crescimento
% 
% % Glicerol = S
% Yxs=0.43;       % [gX/gS] Fator de conv. de substrato em biomassa
% Ks=2;           % [g/L] Constante de afinidade
% 
% % Fonte de nitrogênio = N
% Yxn=7.8;        % [gX/gN] Fator de conv. de nitrogênio em biomassa
% Kn=0.189;       % [g/L] Constante de afinidade
% 
% % PHB = P
% Yps=0.42;       % [gP/gS] Fator de conv. de glicerol em produto
% mu_Pmax=0.097;  % [gP/gXr/h] Velocidade máxima de produção
% Kp=3.496;       % [g/L] Constante de afinidade glicerol-produção
% PXmax=4;        % [-] Capacidade de acúmulo máxima
% theta=0.8;      % [-] Expoente fração acumulada
% 
% 
% Pi=Xri*0.1;     % [g/L] Concentração inicial de PHB
% 
% %Oxigênio dissolvido = C
% Yxo=1.473; %2.7016;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
% Ypo=2.303; %1.9169;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
% mo=1.5*10^-05;  % [gO2/gXr/h] Consumo de O2 para manutenção celular.
% Csat=0.00571;   % [gO2/L] Concentração de saturação
% Ci=Csat*0.90;   % [gO2/L] Concentração inicial
% 
% % ==== AGITAÇÃO E MISTURA ====
% V=4;            % [L] Volume inicial de meio
% Nrot=450;       % [rpm] Velocidade de agitação
% Q=0.2;          % [lpm] Vazão de aeração(1 atm @ 20°C)
% 
% 
% 
% 
% ncond=5;
% tf=[12 18.6 23 27 48];
% 
% tspan=linspace(0,tf(1),tf(1)*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[t,y];
% kLa(1,1)=0; kLa(1,2)=agitacao(Nrot,Q,V);
% kLa(2,1)=tf(1); kLa(2,2)=kLa(1,2);
% 
% for j=2:ncond
%     Xri=y(end,1);Si=y(end,2);Ni=y(end,3);Pi=y(end,4);Ci=y(end,5);
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
% tspan=linspace(tf(j-1),tf(j),(tf(j)-tf(j-1))*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[EDO;t,y];
% kLa(j*2-1,1)=tf(j-1); kLa(j*2-1,2)=agitacao(Nrot,Q,V);
% kLa(j*2,1)=tf(j); kLa(j*2,2)=agitacao(Nrot,Q,V);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % ===== CÁLCULOS DENTRO DA EDO =====








% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CULTIVO E
% mu_Xrmax=0.156; % [h^-1] Vel. específica de crescimento
% 
% % Glicerol = S
% Yxs=0.43;       % [gX/gS] Fator de conv. de substrato em biomassa
% Ks=2;           % [g/L] Constante de afinidade
% 
% % Fonte de nitrogênio = N
% Yxn=7.8;        % [gX/gN] Fator de conv. de nitrogênio em biomassa
% Kn=0.189;       % [g/L] Constante de afinidade
% 
% % PHB = P
% Yps=0.42;       % [gP/gS] Fator de conv. de glicerol em produto
% mu_Pmax=0.097;  % [gP/gXr/h] Velocidade máxima de produção
% Kp=3.496;       % [g/L] Constante de afinidade glicerol-produção
% PXmax=4;        % [-] Capacidade de acúmulo máxima
% theta=0.8;      % [-] Expoente fração acumulada
% 
% 
% Pi=Xri*0.1;     % [g/L] Concentração inicial de PHB
% 
% %Oxigênio dissolvido = C
% Yxo=2.7016;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
% Ypo=1.9169;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
% mo=1.5*10^-05;  % [gO2/gXr/h] Consumo de O2 para manutenção celular.
% Csat=0.00571;   % [gO2/L] Concentração de saturação
% Ci=Csat*0.90;   % [gO2/L] Concentração inicial
% 
% % ==== AGITAÇÃO E MISTURA ====
% V=4;            % [L] Volume inicial de meio
% Nrot=450;       % [rpm] Velocidade de agitação
% Q=0.2;          % [lpm] Vazão de aeração(1 atm @ 20°C)
% 
% 
% 
% ncond=10;
% tf=[15.1 22.7 24 24.75 32.5 33 38.5 44 65.2 82];
% 
% tspan=linspace(0,tf(1),tf(1)*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[t,y];
% kLa(1,1)=0; kLa(1,2)=agitacao(Nrot,Q,V);
% kLa(2,1)=tf(1); kLa(2,2)=kLa(1,2);
% 
% for j=2:ncond
%     Xri=y(end,1);Si=y(end,2);Ni=y(end,3);Pi=y(end,4);Ci=y(end,5);
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
% tspan=linspace(tf(j-1),tf(j),(tf(j)-tf(j-1))*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[EDO;t,y];
% kLa(j*2-1,1)=tf(j-1); kLa(j*2-1,2)=agitacao(Nrot,Q,V);
% kLa(j*2,1)=tf(j); kLa(j*2,2)=agitacao(Nrot,Q,V);
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % ===== CÁLCULOS DENTRO DA EDO =====












%back-up - SensGlobal


% function modelo
% % [t,Xr,P,C]=modelo(mu_Xrmax, Yxs, Ks, Yxn, Kn, Yps, mu_Pmax, Kp,Yxo, Ypo, mo,Csat)
% Xri=0.32; Si=32; Ni=1;
% Salvar=0; Otimizacao=0;
% 
% %%%%%%%%%%%%%%%%%%%%%%%% Cinética microbiana
% 
% mu_Xrmax=0.115; % [h^-1] Vel. específica de crescimento
% 
% % Glicerol = S
% Yxs=0.46;       % [gX/gS] Fator de conv. de substrato em biomassa
% Ks=2;           % [g/L] Constante de afinidade
% 
% % Fonte de nitrogênio = N
% Yxn=6.9;%7.8;        % [gX/gN] Fator de conv. de nitrogênio em biomassa
% Kn=0.189;       % [g/L] Constante de afinidade
% 
% % PHB = P
% Yps=0.46;       % [gP/gS] Fator de conv. de glicerol em produto
% mu_Pmax=0.097;  % [gP/gXr/h] Velocidade máxima de produção
% Kp=3.496;       % [g/L] Constante de afinidade glicerol-produção
% PXmax=4;        % [-] Capacidade de acúmulo máxima
% theta=0.8;      % [-] Expoente fração acumulada
% 
% 
% Pi=Xri*0.1;     % [g/L] Concentração inicial de PHB
% 
% %Oxigênio dissolvido = C
% Yxo=1.473;%2.7016;     % [gXr/gO2] Fator de conversão de oxigênio em célula %1.473; 
% Ypo=2.2;%1.9169;     % [gP/gO2] Fator de conversão de oxigênio em PHB %2.303;
% mo=1.505*10^-05;% [gO2/gXr/h] Consumo de O2 para manutenção celular.
% Csat=0.0069;    % [gO2/L] Concentração de saturação
% Ci=Csat*0.90;   % [gO2/L] Concentração inicial
% 
% % ==== AGITAÇÃO E MISTURA ====
% V=4;            % [L] Volume inicial de meio
% Nrot=450;       % [rpm] Velocidade de agitação
% Q=0.2;          % [lpm] Vazão de aeração(1 atm @ 20°C)
% 
% 
% 
% ncond=1;
% tf=35;
% 
% tspan=linspace(0,tf(1),tf(1)*4+1);
% par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% EDO=[t,y];
% % kLa(1,1)=0; kLa(1,2)=agitacao(Nrot,Q,V);
% % kLa(2,1)=tf(1); kLa(2,2)=kLa(1,2);
% 
% % for j=2:ncond
% %     Xri=y(end,1);Si=y(end,2);Ni=y(end,3);Pi=y(end,4);Ci=y(end,5);
% %     switch j
% %         case 2
% %             Nrot=600;
% %         case 3
% %             Q=0.4; Nrot=750;
% %         case 4
% %             Si=27.33;
% %     end     
% % tspan=linspace(tf(j-1),tf(j),(tf(j)-tf(j-1))*4+1);
% % par=[Yxs Yxn Yps Yxo Ypo mu_Xrmax mu_Pmax Ks Kn Kp PXmax theta mo Csat Nrot Q];
% % [t,y]=ode45(@balancos,tspan,[Xri Si Ni Pi Ci],[],par);
% % EDO=[EDO;t,y];
% % kLa(j*2-1,1)=tf(j-1); kLa(j*2-1,2)=agitacao(Nrot,Q,V);
% % kLa(j*2,1)=tf(j); kLa(j*2,2)=agitacao(Nrot,Q,V);
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % ===== CÁLCULOS DENTRO DA EDO =====
% P=zeros(1,length(EDO)); muP=zeros(1,length(EDO));
%  S=zeros(1,length(EDO)); muS=zeros(1,length(EDO));
%  N=zeros(1,length(EDO)); QO2=zeros(1,length(EDO));
%  C=zeros(1,length(EDO));muXr=zeros(1,length(EDO));
% Xt=zeros(1,length(EDO));  Xr=zeros(1,length(EDO));
% 
% for i=1:length(EDO)
%     t(i)=EDO(i,1); Xr(i)=EDO(i,2); S(i)=EDO(i,3); 
%     N(i)=EDO(i,4); P(i)=EDO(i,5); C(i)=EDO(i,6);
%     
%     muXr(i)=mu_Xrmax*(S(i)/(Ks+S(i)))*(N(i)/(Kn+N(i)));
%     
%     muP(i)=mu_Pmax*(1-N(i)/(N(i)+Kn))*(S(i)/(S(i)+Kp))*(1-((P(i)/Xr(i))/PXmax)^theta);
% 
%     muS(i)=muXr(i)/Yxs+muP(i)/Yps;
%     
%     QO2(i)=muXr(i)/Yxo+muP(i)/Ypo+mo;
%     
%     Xt(i)=Xr(i)+P(i);
% end
% 
% Xr=Xr(:);P=P(:);C=C(:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%% OTIMIZAÇÃO DE PARÂMETROS
% if Otimizacao>0
%     C=C(:); P=P(:); Xt=Xt(:); SSE=[0 0 0];
%     dados_modelo=unique([t,C,P,Xt],'rows');
%     exp_O2=csvread('exp_O2.csv');
%     exp_P =csvread('exp_P.csv');
%     exp_Xt=csvread('exp_Xt.csv');
%     w=[length(exp_O2) length(exp_P) length(exp_Xt)];
%     for j=1:3
%         switch j
%             case 1
%                 dados_exp=exp_O2;
%                 modelo_interp=zeros(1,w(j));
%             case 2
%                 dados_exp=exp_P;
%                 modelo_interp=zeros(1,w(j));
%             case 3
%                 dados_exp=exp_Xt;
%                 modelo_interp=zeros(1,w(j));
%         end
%         for i=1:w(j)
%             modelo_interp(i)=interp1(dados_modelo(:,1),dados_modelo(:,j+1),dados_exp(i,1));
%             SSE(j)=SSE(j)+(modelo_interp(i)-dados_exp(i,2))^2; % Função multiobjetivo
%         end
%     end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % ================ Resultados =================
% if Salvar>0 % Salvar dados em planilha xls
%     muXr=muXr(:);muS=muS(:);muP=muP(:);Xt=Xt(:);QO2=QO2(:);
%     Dados=[EDO Xt muXr muS muP QO2];
%     xlswrite('TransfO2.xlsx',Dados,1,'A4');
%     xlswrite('TransfO2.xlsx',kLa,1,'O4');
% end
% % disp('Fim.');
% end
