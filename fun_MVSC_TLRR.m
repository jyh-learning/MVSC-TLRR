function [L,E,L1,L2,L3]=fun_MVSC_TLRR(X,para)
[~,n,v]=size(X);
w1=para.w1;
w2=1-w1;
lambda=para.lambda;
alpha=para.alpha; 


%% init 
tol = 1e-4; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;


L = zeros(n,n,v);
L1 = L;
L2 = L;
L3 = L;
E=L;

Lambda11=L;
Lambda12=L;
Lambda13=L;


Lambda21=L;
Lambda22=L;
Lambda23=L;



for iter = 1 : max_iter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % update L1 or the Frontal slcies
    Temp=0.5*(X+L-E+(Lambda11+Lambda21)/mu);
    for i=1:v
        Tempvector=0.5*(Temp(:,:,i)+Temp(:,:,i)');
        L1(:,:,i)=prox_nuclear(Tempvector,w1/(2*mu));
    end
    % update L2 or the Horizontal slcies
    Temp=0.5*(X+L-E+(Lambda12+Lambda22)/mu);

    for i=1:n
        Tempvector=Temp(:,i,:);
        Tempvector=reshape(Tempvector,n,v);
        L2(:,i,:)=reshape(prox_nuclear(Tempvector,w2/(2*mu)),n,1,v);
    end
    Temp=0.5*(X+L-E+(Lambda13+Lambda23)/mu);
    for i=1:n
        Tempvector=Temp(:,i,:);
        Tempvector=reshape(Tempvector,n,v);
        L3(:,i,:)=reshape(prox_l21(Tempvector,(w2*alpha)/(2*mu)),n,1,v);
    end

    % update E
    Temp=(1/3)*(3*X-L1-L2-L3+(Lambda11+Lambda12+Lambda13)/mu);
    E=Temp/(1+(2*lambda)/(3*mu));
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   % update L and E
    L=(1/3)*(L1+L2+L3-(Lambda21+Lambda22+Lambda23)/mu);

    
  
    d11=X-L1-E;
    d12=X-L2-E;
    d13=X-L3-E;

    d21=(L-L1);
    d22=(L-L2);
    d23=(L-L3);


    chg = max([ max(abs(d11(:))),max(abs(d12(:))),max(abs(d13(:))),max(abs(d21(:))),max(abs(d22(:))),max(abs(d23(:))) ]);
    
    disp(['#iteration: ',num2str(iter),'. chg: ',num2str(chg)]);
    if chg < tol
        break;
    end 
    
    

    Lambda11=Lambda11+mu*d11;
    Lambda12=Lambda12+mu*d12;
    Lambda13=Lambda13+mu*d13;
    
    Lambda21=Lambda21+mu*d21;
    Lambda22=Lambda22+mu*d22;
    Lambda23=Lambda23+mu*d23;
    
    mu = min(rho*mu,max_mu);    
end