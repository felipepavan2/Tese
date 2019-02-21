clc;clear;tic;

% amostra=csvread('teste.csv');
% Y=csvread('teste2.csv');
% [N a]=size(amostra);
% Rx=zeros(12,1);
% Ry=zeros(12,1);
% for j=1:a
%     for i=1:N
%         x(:,1)=amostra(i,:);
%         for k=1:a
%             y(k,1)=Y(i,1);
%         end
%         if i==1
%             z=amostra(i,:);
%             z(:,j)=[];
%             z=[1 z];
%             z=z';
%         end
%         [b,bint,rx] = regress(x,z);
%         [b,bint,ry] = regress(y,z);
%         
%         Rx=[Rx rx];
%         Ry=[Ry ry];
%     end
%     
% end
% Rx=Rx';Ry=Ry';
% xlswrite('Teste.xlsx',Rx,1,'A1');
% xlswrite('Teste.xlsx',Ry,2,'A1');

% 
% 
%   z=[1 156 420 335 823 964 769 68 422 336 440 510];
% x=[354 156 420 335 823 964 769 68 422 336 440 510];    
% y=[785 785 785 785 785 785 785 785 785 785 785 785];
% x=x(:);z=z(:);y=y(:);
% 
% [b,bint,rx] = regress(x,z)
% [b,bint,ry] = regress(y,z)

amostra=csvread('SensGlobal5000.csv');
N=length(amostra);

toc;
Y=[];

for i=1:N % Executa o modelo com os valores gerados por LHS
    Csat=amostra(i,1); Kpn=amostra(i,2); Kps=amostra(i,3); Kxn=amostra(i,4); Kxs=amostra(i,5); 
    Ypo=amostra(i,6);Yps=amostra(i,7); Yxn=amostra(i,8);Yxo=amostra(i,9);Yxs=amostra(i,10);
    mo=amostra(i,11);ms=amostra(i,12);mu_Pmax=amostra(i,13);mu_Xrmax=amostra(i,14); theta=amostra(i,15);
    
   [t,Xr,P,C]=modelo(Csat, Kpn, Kps, Kxn, Kxs, Ypo, Yps, Yxn, Yxo, Yxs, mo, ms, mu_Pmax, mu_Xrmax, theta);
   
%    if i==1
%        Y1=[Xr];Y2=[P];Y3=[C];Y4=[S];Y5=[N];
%        
%    else
%        Y1=[Y1,Xr];Y2=[Y2,P];Y3=[Y3,C];Y4=[Y4, S];Y5=[Y5, N];
%    end

Y=[Y,Xr,P,C];

end

% xlswrite('Yoi.xlsx',Y1,1,'A1'); xlswrite('Yoi.xlsx',Y2,2,'A1'); xlswrite('Yoi.xlsx',Y3,3,'A1');
% xlswrite('Yoi.xlsx',Y4,4,'A1'); xlswrite('Yoi.xlsx',Y5,5,'A1');
toc;
x=zeros(N,1); y=zeros(N,1);[a,b]=size(Y);PRCC=zeros(a,2);pvalue=zeros(a,2);
for i=1:3               % Para cada output
    for j=1:15          % Para cada parâmetro
        for m=1:a       % Para cada passo de tempo
            for k=1:N   % Para cada execução
                x(k,1)=amostra(k,j);
                if k==1
                    y(k,1)=Y(m,k+i-1);
                    z=amostra;
                    z(:,j)=[];
                else
                    y(k,1)=Y(m,k*3-3+i);
                    if k==N
                        [rho,p]=partialcorr(x,y,z,'Type','Spearman');
                        PRCC(m,j)=rho;
                         pvalue(m,j)=p;
                    end
                end
            end
        end
    end
    xlswrite('SensGlobal5000.xlsx',t,i  ,'A3');
    xlswrite('SensGlobal5000.xlsx',PRCC,i,'B3'); 
end



% [a b]=size(amostra);
% for j=1:b
%     [s,i]=sort(amostra(:,j));
%     r(i,j)=[1:a]';
% end
% xlswrite('SensGlobalPostos.xlsx',r,1,'A3');
% 
% 
% [a b]=size(Y);
% for j=1:b
%     [s,i]=sort(Y(:,j));
%     r(i,j)=[1:a]';
% end
% xlswrite('SensGlobalPostos.xlsx',r,2,'A3');

% n=length(x);
% [s,i]=sort(x);                                                
% r(i,1)=[1:n]';


toc;