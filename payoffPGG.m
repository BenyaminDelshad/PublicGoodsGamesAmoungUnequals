function [pi,coop,Mat]=payoffPGG(p1,p2,p3,evec,rvec,Xset,nPlayer)

% pi,coop ... payoffs and cooperation rates of the three players
% p1,p2,p3 ... players' strategies. 
% Old Version: The strategy entries p_{ijk} give the probability the player cooperates
%%% Update: p1,p2,p3 are two dimentional arrays give the probability the
%%% player choose cooperation rate x_i from Xset, 8 * 2 for two cooperation
%%% rate options, 27 * 3 for three, 64 * 4 for four and ...
% given the players 1,2,3 chose the actions i,j,k in {C,D} in previous
% round. 
% evec=(e1,e2,e3) ... individual endowments
% rvec=(r1,r2,r3) ... individual multiplication factors

%% Parameters and preparations
pi=zeros(1,nPlayer); coop=zeros(1,nPlayer); 
% OS=[1 1 1; 1 1 0; 1 0 1; 1 0 0; 0 1 1; 0 1 0; 0 0 1; 0 0 0]; % Possible outcomes of one-shot game
OS=zeros(length(Xset)^nPlayer,nPlayer);
index = 1;
for i=1:length(Xset)
    for j=1:length(Xset)
        for k=1:length(Xset)
            OS(index,:) = [Xset(i), Xset(j), Xset(k)];
            OS2(index,:) = [i,j,k]; % keep index instead of exact value because we need it!
            index = index + 1;
        end
    end
end
nOS=size(OS,1); % number of outcomes
PayOS=zeros(nOS,nPlayer); % Possible one-shot payoffs for the three players
yvec=evec.*rvec; % players' absolute contributions
for i=1:nOS
    s=OS(i,:);
    PayOS(i,1)=s*yvec'/nPlayer+(1-s(1))*evec(1);
    PayOS(i,2)=s*yvec'/nPlayer+(1-s(2))*evec(2); 
    PayOS(i,3)=s*yvec'/nPlayer+(1-s(3))*evec(3);
end


%% Constructing the transition matrix Mat
Mat=ones(nOS,nOS);
for iOld=1:nOS
    so=OS2(iOld,:); 
    pC1=p1(iOld,:); pC2=p2(iOld,:); pC3=p3(iOld,:); 
    
    for iNew=1:nOS
        sn=OS2(iNew,:);
        %pr1=pC1*sn(1)+(1-pC1)*(1-sn(1)); 
        %pr2=pC2*sn(2)+(1-pC2)*(1-sn(2)); 
        %pr3=pC3*sn(3)+(1-pC3)*(1-sn(3));
        pr1=pC1(sn(1));
        pr2=pC2(sn(2));
        pr3=pC3(sn(3));
        Mat(iOld,iNew)=pr1*pr2*pr3; 
    end
end

v=null(Mat'-eye(nOS)); v=v/sum(v);
pi=PayOS'*v; coop=OS'*v; 
dM = det(Mat);
end