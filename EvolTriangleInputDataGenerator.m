%This code must give EVec as an Input and generate
%EVec, Coop ,Pay ,rvec ,nGen and s
function [EVec,Coop,Pay,rvec,nGen,s]=EvolTriangleInputDataGenerator(EVec, rvec,Xset1,Xset2,Xset3,s,nGen)
%Initializing and defining
Coop = EVec;
Pay = EVec;
sz = size(EVec);
nQuery = sz(1);
%fill Coop and Pay Matrices
for i=1:nQuery
    [~,~,~,AvCoop,AvPi,~,~,rvec,Xset1,Xset2,Xset3,s,nGen]=EvolProc(EVec(i,:),rvec,Xset1,Xset2,Xset3,s,nGen);
    Coop(i,:) = AvCoop;
    Pay(i,:) = AvPi;
end
end
