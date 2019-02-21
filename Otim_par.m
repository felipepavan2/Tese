function Otim_par()
clc; clear;
tic;
% 
% % v=[-0.1745 -0.5095 0.5095 -0.0945 -0.2568];
% v=[-0.0787 -3.4152 37.8897];
% 
% options = optimset('PlotFcns',@optimplotfval);
% 
% v_otim=fminsearch(@modelo,v,options);
% 
% assignin('base','v',v_otim); % Adicionar ao workspace
    exp_O2=csvread('exp_O2.csv');
    exp_P =csvread('exp_P.csv');
    exp_Xt=csvread('exp_Xt.csv');
    VDados=csvread('exp_Volume.csv');
    
% for t=0:0.1:82
%     i=round(t*10+1);
%      dvt(i) = interp1(VDados(:,1),VDados(:,2),t);
%        V(i) = interp1(VDados(:,1),VDados(:,3),t);
%      Volume(i,:)=[t,dvt(i),V(i)];
% end

    save('Variaveis','exp_O2','exp_P','exp_Xt','VDados')
    
   
v=[0.3   1   0.5  7.8 0.41 0.09 0.15  1];
lb=[0.05 0.5 0.05 6.3 0.4  0.05 0.125 0.01];
ub=[0.3  3   2.5  8.3 0.5  0.1  0.17  2];

meta=[0 0 0];
peso=[1 1 1];
options = optimoptions('fgoalattain','Display','iter','FiniteDifferenceStepSize',1e-4);
 options.MaxIterations = 1000;
 options.MaxFunctionEvaluations = 800;
 options.PlotFcns =@optimplotfval;
[v_otim,fval]=fgoalattain(@(v)modelo(v),v,meta,peso,[],[],[],[],lb,ub,[],options)
assignin('base','v',v_otim); % Adicionar ao workspace
toc;
end

% v_otim =    0.1660    5.0000    0.4039    1.0000    6.9864    1.8565    0.4075    0.0898    0.1580   -0.5292
