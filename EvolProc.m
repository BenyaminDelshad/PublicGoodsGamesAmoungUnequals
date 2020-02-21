function [x1T,x2T,x3T,AvCoop,AvPi,nInv,evec,rvec,Xset,s,nGen]=EvolProc(evec,rvec,Xset,s,nGen)

% Old version: Structure xiT: 1st-8th column: mem-1 strategy, 9th column: coopRate, 10th
% column: payoff

%%% Update: last comments is only for the case that we have only two
%%% possible action(size of Xset = 2) for example Defect(0) or
%%% cooperate(1), but now, In general form, we can have more actions for
%%% each player to do! so XiT have |Xset|^4 + 2 coloumns. two last columns 
%%% are for coopRate and payoff and others are for mem-1 strategy!

%%Description of one row of XiT in 3-plaayer mode: player i have strategy p_t(|Xset|^3,|Xset|)
%%which we describe it as one dinemtional |Xset|^4 array and after that
%%coopRate and after that, payoff.

% evec ... (e1,e2,e3) distribution of endowments
% rvec ... (r1,r2,r3) players' productivities 
% s ... selection strength, nGen ... number of strategy updates

%% Setting up all objects
C=clock; rng(C(5)*60+C(6)); % Setting up the random number generator
% initilize starting strategies for all players, in this case, all of them
% choose first action of Xset with probability 1.
nPlayer=3; % But this code is not simply upgradable for using other values instead of 3 yet.
p1=zeros(length(Xset)^nPlayer,length(Xset)); p2=p1; p3=p1;
for i=1:(length(Xset)^nPlayer)
    p1(i,1) = 1;
end
p2 = p1;
p3 = p1;

% end initialing p1, p2, p3
x1T=zeros(nGen,length(Xset)^(nPlayer + 1) + 2); x2T=x1T; x3T=x1T;
p_tmp = DimentionReducer(p1);
for i=1:length(Xset)^(nPlayer + 1)
    x1T(1,i) = p_tmp(i);
end
x2T = x1T;
x3T = x1T;

[pi,coop]=payoffPGG(p1,p2,p3,evec,rvec,Xset,nPlayer);
% bug report: xTi instead of xiT ...!
x1T(1,length(Xset)^(nPlayer + 1) + 1:length(Xset)^(nPlayer + 1) + 2)=[coop(1), pi(1)]; 
x2T(1,length(Xset)^(nPlayer + 1) + 1:length(Xset)^(nPlayer + 1) + 2)=[coop(2), pi(2)]; 
x3T(1,length(Xset)^(nPlayer + 1) + 1:length(Xset)^(nPlayer + 1) + 2)=[coop(3), pi(3)];

nInv=zeros(1,nPlayer); 

%% Evolutionary process
for i=1:nGen-1
    % Randomly choosing a player
    zR=rand(1); iPop=floor(3*zR)+1; 
    pMut=rand(length(Xset)^nPlayer,length(Xset)); % Memory-1
    for j=1:length(Xset)^nPlayer
        pMut(j,:) = pMut(j,:) / sum(pMut(j,:)); % sum of probabilities must be 1!
    end
    % Player 1 explores a new strategy
    if iPop==1 
        [pi,coop]=payoffPGG(pMut,p2,p3,evec,rvec,Xset,nPlayer);
        rho=1/(1+exp(-s*(pi(1)-x1T(i,end))));
        rn=rand(1);
        if rn<rho
            p1=pMut;
            
            x1T(i+1,:)=[DimentionReducer(p1),coop(1),pi(1)];
            x2T(i+1,:)=[DimentionReducer(p2),coop(2),pi(2)];
            x3T(i+1,:)=[DimentionReducer(p3),coop(3),pi(3)];
            nInv(1)=nInv(1)+1;
        else
            x1T(i+1,:)=x1T(i,:);
            x2T(i+1,:)=x2T(i,:);
            x3T(i+1,:)=x3T(i,:); 
        end
    
    % Player 2 explores a new strategy
    elseif iPop==2 
        [pi,coop]=payoffPGG(p1,pMut,p3,evec,rvec,Xset,nPlayer);
        rho=1/(1+exp(-s*(pi(2)-x2T(i,end))));
        rn=rand(1);
        if rn<rho
            p2=pMut;
            x1T(i+1,:)=[DimentionReducer(p1),coop(1),pi(1)];
            x2T(i+1,:)=[DimentionReducer(p2),coop(2),pi(2)];
            x3T(i+1,:)=[DimentionReducer(p3),coop(3),pi(3)];
            nInv(2)=nInv(2)+1;
        else
            x1T(i+1,:)=x1T(i,:);
            x2T(i+1,:)=x2T(i,:);
            x3T(i+1,:)=x3T(i,:);
        end
        
    % Player 3 explores a new strategy
    else
        [pi,coop]=payoffPGG(p1,p2,pMut,evec,rvec,Xset,nPlayer);
        rho=1/(1+exp(-s*(pi(3)-x3T(i,end))));
        rn=rand(1);
        if rn<rho
            p3=pMut;
            x1T(i+1,:)=[DimentionReducer(p1),coop(1),pi(1)];
            x2T(i+1,:)=[DimentionReducer(p2),coop(2),pi(2)];
            x3T(i+1,:)=[DimentionReducer(p3),coop(3),pi(3)];
            nInv(3)=nInv(3)+1;
        else
            x1T(i+1,:)=x1T(i,:);
            x2T(i+1,:)=x2T(i,:);
            x3T(i+1,:)=x3T(i,:);
        end
        
    end
end

AvCoop=[mean(x1T(:,end-1)), mean(x2T(:,end-1)), mean(x3T(:,end-1))]; 
AvPi=[mean(x1T(:,end)), mean(x2T(:,end)), mean(x3T(:,end))]; 
end