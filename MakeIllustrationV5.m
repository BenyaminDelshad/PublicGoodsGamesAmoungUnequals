% You must change input file address every time you want to use this code
% on line 67.
function[sum_pay]= MakeIllustrationV5();

%% Fixing the Parameters
global lw ms fsA fsL fsH fsT fname le bo wi he dx dy col wi2 he2 dx2 b_factor

b_factor = 8/3;
lw=3; ms=10; fsT=10; fsA=12; fsL=12; fsH=16; fname='Arial'; 
le=0.076; bo=b_factor*0.105; wi=0.165; he=wi*b_factor*90/60; wi2=wi*0.9; he2=he*0.9; dx2=wi*0.05; dx=0.07; dy=b_factor*0.4; % change dy from 0.365 to 0.55
col=[0 0.447 0.741;
    0.85 0.325 0.098;
    0.494 0.184 0.556;
    0.929 0.694 0.125;
    0.466 0.674 0.188];

%% Plotting the general graphical structure 
figure('Position',[300,300,1200,800 / b_factor]); % 300 - 300 - 1200 - 750 
ax1=axes('Position',[le-0.35*dx 0 3*wi+3*dx 1]); hold on
axis([-0.01 1.01 0 1]); dxf=(wi+dx*0.4)/(3*wi+2*dx); delf=0.27*dx/(3*wi+2*dx);
grey1=242/256*ones(1,3);  
blue1=[79,129,189]/256; 
y1=0.005; y2=0.37; lwP=2;%change y2 from 0.735 to 0.94!
rd=0.01; t1=pi:0.01:3*pi/2; t2=3*pi/2:0.01:2*pi; t3=0:0.01:pi/2; t4=pi/2:0.01:pi; 
bxx=[rd+rd*cos(t1),dxf-rd+rd*cos(t2),dxf-rd+rd*cos(t3),rd+rd*cos(t4)]; 
bxy=[y1+rd+rd*sin(t1),y1+rd+rd*sin(t2),y2-rd+rd*sin(t3),y2-rd+rd*sin(t4)]; 
%fill(bxx,bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);
%fill(bxx+dxf+delf,bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);
%fill(bxx+2*(dxf+delf),bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);

%fill(bxx+3*(dxf+delf),bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);

% y1=0.76; y2=0.94; lwP=2;
% rd=0.01; t1=pi:0.01:3*pi/2; t2=3*pi/2:0.01:2*pi; t3=0:0.01:pi/2; t4=pi/2:0.01:pi; 
% bxx=[rd+rd*cos(t1),dxf-rd+rd*cos(t2),dxf-rd+rd*cos(t3),rd+rd*cos(t4)]; 
% bxy=[y1+rd+rd*sin(t1),y1+rd+rd*sin(t2),y2-rd+rd*sin(t3),y2-rd+rd*sin(t4)]; 
% fill(bxx,bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);
% fill(bxx+dxf+delf,bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);
% fill(bxx+2*(dxf+delf),bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);

%fill(bxx+3*(dxf+delf),bxy,grey1,'EdgeColor',blue1,'LineWidth',lwP);

S1=[0.065-0.028 0.420]; red=[0.8 0 0];  msS=8; lwS=2;
%plot(S1(1),S1(2),'o','Color',red,'MarkerSize',msS,'LineWidth',lwS); 
%text(S1(1)+0.015,S1(2)-0.0075,{'Extreme endowment inequality','prevents cooperation '},'Color',red,...
%    'FontName',fname,'FontSize',fsT,'VerticalAlignment','middle'); 
% S2=S1+[2*(wi+dx) 0]; 
% plot(S2(1),S2(2),'o','MarkerSize',msS,'LineWidth',lwS,'Color',col(3,:));
% plot(S2(1),S2(2),'x','MarkerSize',msS,'Color',col(3,:));
% text(S2(1)+0.015,S2(2)-0.0075,{'Complete endowment equality','prevents cooperation '},'Color',col(3,:),...
%     'FontName',fname,'FontSize',fsT,'VerticalAlignment','middle'); 


he1=0.03; wi1=3*wi+3*dx; k=wi1;
rd=0.006; t=0:0.01:2*pi; red=[0.8 0 0];  msS=8; 
axis off

% PlotGameSetup(1,'a'); 
% PlotGameSetup(2,'b'); 
% PlotGameSetup(3,'c'); 

%% Plotting the static predictions
delta1=0.8; delta2=0.3; delta3=0.35;
%Input_file = 'Inputs/InputDataXsetSize3_11.mat';
%Input_file = 'Inputs/InputDataXsetSize3and2_10.mat';
%Input_file = 'TempInputData3.mat';
%Input_file = 'Inputs/InputData_exact_Fig2_1MnGen.mat';
Input_file = 'Inputs/InputDataXsetSize2_15_ExactThreshold_2.mat';

load(Input_file);
%rvec1=[2 2 2]; rvec2=[1.1 1.5 2.9]; 
rvec1 = rvec;
rvec2 = rvec;
%plotEvo2(1,delta1,rvec1,'a',Input_file); 
% plotEvo2(2,delta2,rvec2,'b',Input_file);
% plotEvo2(3,delta3,rvec2,'c',Input_file);
%plotStaticLegend(); 

%% Plotting the evolutionary predictions
delta1=1; delta2=1; delta3=1;
plotEvo(0,delta1,rvec1,'a',Input_file); 
plotEvo(1,delta1,rvec1,'b',Input_file); 
plotEvo(2,delta2,rvec2,'c',Input_file);
plotEvo(3,delta2,rvec2,'d',Input_file);

end


function PlotGameSetup(nr,lett); 
global lw ms fsA fsL fsH fsT fname le bo wi he dx dy col

ax1=axes('Position',[le+(wi+dx)*(nr-1) bo+2*dy wi he]);
axis([-1.05 1.05 -sqrt(3)*0.1 sqrt(3)*1.1]); hold on

if nr==1, head1='Head!'; 
elseif nr==2, head1='Asymmetric and linear'; 
elseif nr==3, head1='Symmetric and nonlinear'; 
end
text(0,1.08,head1,'FontSize',fsA,'FontName',fname,...
    'HorizontalAlignment','center','FontWeight','bold');
text(0,0.9,'public goods game','FontSize',fsA,'FontName',fname,...
        'HorizontalAlignment','center','FontWeight','bold');
axis off
end


function plotEvo(nr,delta,rvec,lett,datafile); 
global lw ms fsA fsT fsL fsH fname le bo wi he dx col wi2 he2 dx2 b_factor

ax1=axes('Position',[le+(wi+dx)*(nr)+dx2 bo wi2 he2]); 
axis([-1.05 1.05 -sqrt(3)*0.1 sqrt(3)*1.1]); hold on
edges=[-1 0; 1 0; 0 sqrt(3)]; 
x1=edges(1,:); x2=edges(2,:); x3=edges(3,:);
%rvec
load(datafile);
sumpi=sum(Pay,2);
% group payoff previously but will changes in next lines
% to Cooperation rates of each player!

if(nr == 1)
    sumpi = Coop(:,1); %Cooperation rates of Player 1!
    for i=1:length(sumpi)
        if (sumpi(i) < 0)
            sumpi(i)=0;
        end
    end
end
if(nr == 2)
    sumpi = Coop(:,2); %Cooperation rates of Player 2!
    for i=1:length(sumpi)
        if (sumpi(i) < 0)
            sumpi(i)=0;
        end
    end
end
if(nr == 3)
    sumpi = Coop(:,3); %Cooperation rates of Player 3!
    for i=1:length(sumpi)
        if (sumpi(i) < 0)
            sumpi(i)=0;
        end
    end
end
%sumpi
PiMax=max(sumpi);
PiMin=min(sumpi);

%sumpi
%PiMin
%PiMax

%rvec
if nr==0
text(-1.9,sqrt(3)/2,'Evolutionary analysis','FontSize',fsL,'Color',col(5,:),...
    'FontName',fname,'Rotation',90,'HorizontalAlignment','center','FontWeight','bold');
end


%% Plotting the contour plot
dxx=0.48; dyy=0.45; ss=0.1; epsi=0.03;



for i=size(EVec,1):-1:1;
    %EVec(i,:)
    %debug infos to find why we have white areas?!
%     tmp_val = sum(Pay(i,:));
%     sum_pay(i) = tmp_val;
%     if (tmp_val<1.1)
%         EVec(i,:)
%         tmp_val
%         i
%     end
    E1=EVec(i,1:3); E2=(1-epsi)*E1+epsi*(1-E1); 
    e=E2(1)*x1+E2(2)*x2+E2(3)*x3; 
    f1=e(1)+ss*[-dxx dxx dxx-dxx*(E1(1)==0) -dxx+dxx*(E1(2)==0)]; f2=e(2)+ss*[-dyy -dyy dyy dyy]-0.02;
    cl=getcolor(sumpi(i),nr, PiMin, PiMax); 
%     for index=1:4
%         if (f1(index)==-0.236 && f2(index)==0.6382)
%             E1
%         end
%     end
    %[i, cl]
    
    fill(f1,f2,cl,'LineStyle','none'); 
    
end

%% Plotting the triangle
plot([edges(:,1); edges(1,1)],[edges(:,2); edges(1,2)],'k','LineWidth',lw); 
plot(edges(:,1),edges(:,2),'ko','MarkerSize',ms,'MarkerFaceColor','k','LineWidth',lw); 
plot(mean(edges(:,1)),mean(edges(:,2)),'kx','MarkerSize',ms,'MarkerFaceColor','k');

%% Labels
dyL=0.26; 
text(-1,2,lett,'FontSize',fsH,'FontName',fname,'FontWeight','bold');
text(x1(1)+0.175,x1(2)-dyL-0.035,{'Full endowment','to player 1'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
text(x2(1)-0.175,x2(2)-dyL-0.035,{'Full endowment','to player 2'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
text(x3(1),x3(2)+dyL+b_factor*0.035,{'Full endowment','to player 3'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
axis off

%% Color legend
%if nr==1 
%    xmax=PiMax; xmin=1; xtl={'1.00','1.25','1.50','1.75'};
%else
%    xmax=PiMax; xmin=1.2; xtl={'1.2','1.6','2.0','2.4'}; 
%end
xmax=PiMax; xmin=PiMin;
if(nr ~= 0)
    xmax = 1;
    xmin = 0;
end
xtl={round(xmin,2),round(xmin + ((xmax-xmin)/3),2),round(xmax - ((xmax-xmin)/3),2),round(xmax,2)};

ax2=axes('Position',[le+(wi+dx)*(nr) bo-b_factor*0.055 wi 0.01*b_factor],'XTick',0:1/3:1,...
    'XTickLabel',xtl,'FontSize',fsT,'FontName',fname,'YTick',[]); 
box(ax2,'on'); hold on
axis([-0.05 1.05 -1 1]); 

ss=0.01; 
for x=0:ss:1; 
    cl=getcolor((1-x)*xmin+x*xmax,nr, PiMin, PiMax);
    fill(x+ss*[-0.55 0.55 0.55 -0.55],[-0.5 -0.5 0.5 0.5],cl,'LineStyle','none');
end
if(nr == 0)
    xlabel('Group Payoff','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
if(nr == 1)
    xlabel('Player 1 Contributions','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
if(nr == 2)
    xlabel('Player 2 Contributions','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
if(nr == 3)
    xlabel('Player 3 Contributions','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
end



function plotEvo2(nr,delta,rvec,lett, datafile);  % Prevously plotStatic!
global lw ms fsA fsT fsL fsH fname le bo wi he dx dy col wi2 he2 dx2

ax1=axes('Position',[le+(wi+dx)*(nr-1)+dx2 bo+dy wi2 he2]); 
axis([-1.05 1.05 -sqrt(3)*0.1 sqrt(3)*1.1]); hold on
edges=[-1 0; 1 0; 0 sqrt(3)]; 
x1=edges(1,:); x2=edges(2,:); x3=edges(3,:);

% if nr==1
% text(-1.9,sqrt(3)/2,'Equilibrium analysis','FontSize',fsL,'Color',col(1,:),...
%     'FontName',fname,'Rotation',90,'HorizontalAlignment','center','FontWeight','bold');
% end
fill([edges(:,1); edges(1,1)],[edges(:,2); edges(1,2)],'w','LineWidth',lw,'EdgeColor','k'); 


% Copy plotEvo code part to edit here:

load(datafile);
sumpi=sum(Pay,2);
if(nr == 2)
    sumpi = Coop(:,1); %Cooperation rates of Player 1!
    for i=1:length(sumpi)
        if (sumpi(i) < 0)
            sumpi(i)=0;
        end
    end
end
if(nr == 3)
    sumpi = Coop(:,2); %Cooperation rates of Player 2!
    for i=1:length(sumpi)
        if (sumpi(i) < 0)
            sumpi(i)=0;
        end
    end
end

PiMax=max(sumpi);
PiMin=min(sumpi);

if nr==1
text(-1.9,sqrt(3)/2,'Evolutionary analysis','FontSize',fsL,'Color',col(5,:),...
    'FontName',fname,'Rotation',90,'HorizontalAlignment','center','FontWeight','bold');
end


%% Plotting the contour plot
dxx=0.48; dyy=0.45; ss=0.1; epsi=0.03;



for i=size(EVec,1):-1:1;

    E1=EVec(i,1:3); E2=(1-epsi)*E1+epsi*(1-E1); 
    e=E2(1)*x1+E2(2)*x2+E2(3)*x3; 
    f1=e(1)+ss*[-dxx dxx dxx-dxx*(E1(1)==0) -dxx+dxx*(E1(2)==0)]; f2=e(2)+ss*[-dyy -dyy dyy dyy]-0.02;
    cl=getcolor(sumpi(i),nr, PiMin, PiMax); 
    
    fill(f1,f2,cl,'LineStyle','none'); 
    
end

%% Plotting the triangle
plot([edges(:,1); edges(1,1)],[edges(:,2); edges(1,2)],'k','LineWidth',lw); 
plot(edges(:,1),edges(:,2),'ko','MarkerSize',ms,'MarkerFaceColor','k','LineWidth',lw); 
plot(mean(edges(:,1)),mean(edges(:,2)),'kx','MarkerSize',ms,'MarkerFaceColor','k');

%% Labels
dyL=0.26; 
text(-1,2,lett,'FontSize',fsH,'FontName',fname,'FontWeight','bold');
text(x1(1)+0.175,x1(2)-dyL,{'Full endowment','to player 1'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
text(x2(1)-0.175,x2(2)-dyL,{'Full endowment','to player 2'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
text(x3(1),x3(2)+dyL+0.035,{'Full endowment','to player 3'},'FontSize',fsT,...
    'FontName',fname,'HorizontalAlignment','center'); 
axis off

%% Color legend

xmax=PiMax; xmin=PiMin; xtl={round(xmin,2),round(xmin + ((xmax-xmin)/3),2),round(xmax - ((xmax-xmin)/3),2),round(xmax,2)};

ax2=axes('Position',[le+(wi+dx)*(nr-1) (bo+dy)-0.05 wi 0.01],'XTick',0:1/3:1,...
    'XTickLabel',xtl,'FontSize',fsT,'FontName',fname,'YTick',[]); 
box(ax2,'on'); hold on
axis([-0.05 1.05 -1 1]); 

ss=0.01; 
for x=0:ss:1; 
    cl=getcolor((1-x)*xmin+x*xmax,nr, PiMin, PiMax);
    fill(x+ss*[-0.55 0.55 0.55 -0.55],[-0.5 -0.5 0.5 0.5],cl,'LineStyle','none');
end
if(nr == 1)
    xlabel('Group payoff','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
if(nr == 2)
    xlabel('Player 1 Contributions','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end
if(nr == 3)
    xlabel('Player 2 Contributions','FontSize',fsT,'FontName',fname,'Position',[0.5,-6]);
end


% End copy plotEvo part to edit here! start main code of this part in
% prevous code:

%% Calculating the thresholds for the players' endowments
% if nr==3
%     M1=(x1(1)+x2(1)+x3(1))/3; M2=(x1(2)+x2(2)+x3(2))/3;
%     p1=[7/11 4/11 0]; p2=[5/13 5/13 3/13]; p3=[4/11 7/11 0]; 
%     P1=p1(1)*x1+p1(2)*x2+p1(3)*x3; P2=p2(1)*x1+p2(2)*x2+p2(3)*x3; P3=p3(1)*x1+p3(2)*x2+p3(3)*x3; 
%     fill([P1(1) P2(1) P3(1)],[P1(2) P2(2) P3(2)],col(1,:),'LineStyle','none');
%     %plot([P2(1) M1],[P2(2) M2],'k--','LineWidth',2);
%     
%     p1=[7/11 0 4/11]; p2=[5/13 3/13 5/13]; p3=[4/11 0 7/11]; 
%     P1=p1(1)*x1+p1(2)*x2+p1(3)*x3; P2=p2(1)*x1+p2(2)*x2+p2(3)*x3; P3=p3(1)*x1+p3(2)*x2+p3(3)*x3; 
%     fill([P1(1) P2(1) P3(1)],[P1(2) P2(2) P3(2)],col(1,:),'LineStyle','none');
%     %plot([P2(1) M1],[P2(2) M2],'k--','LineWidth',2);
%     
%     p1=[0 7/11 4/11]; p2=[3/13 5/13 5/13]; p3=[0 4/11 7/11]; 
%     P1=p1(1)*x1+p1(2)*x2+p1(3)*x3; P2=p2(1)*x1+p2(2)*x2+p2(3)*x3; P3=p3(1)*x1+p3(2)*x2+p3(3)*x3; 
%     fill([P1(1) P2(1) P3(1)],[P1(2) P2(2) P3(2)],col(1,:),'LineStyle','none');
%     %plot([P2(1) M1],[P2(2) M2],'k--','LineWidth',2);
%     
%     %text(mean(edges([1,3],1))+0.12,mean(edges([1,3],2)-0.08),{'C. f.'},'FontSize',fsT,...
%     %'FontName',fname,'HorizontalAlignment','center','Color','w');
%     %text(mean(edges([2,3],1))-0.11,mean(edges([2,3],2)-0.08),{'C. f.'},'FontSize',fsT,...
%     %'FontName',fname,'HorizontalAlignment','center','Color','w');
%     %text(mean(edges([1,2],1)),mean(edges([1,2],2)+0.1),{'C. f.'},'FontSize',fsT,...
%     %'FontName',fname,'HorizontalAlignment','center','Color','w');
%    
% elseif nr<3
%     e12=delta*rvec(2)/(3-rvec(1)+delta*rvec(2)),
%     e13=delta*rvec(3)/(3-rvec(1)+delta*rvec(3)),
%     A10=x1; A12=e12*x1+(1-e12)*x2; A13=e13*x1+(1-e13)*x3;
% 
%     e21=delta*rvec(1)/(3-rvec(2)+delta*rvec(1)), 
%     e23=delta*rvec(3)/(3-rvec(2)+delta*rvec(3)), 
%     A20=x2; A21=(1-e21)*x1+e21*x2; A23=e23*x2+(1-e23)*x3;
% 
%     e31=delta*rvec(1)/(3-rvec(3)+delta*rvec(1)),
%     e32=delta*rvec(2)/(3-rvec(3)+delta*rvec(2)),
%     A30=x3; A31=(1-e31)*x1+e31*x3; A32=(1-e32)*x2+e32*x3;
% 
%     if nr==1
%         fill([A12(1),A21(1),A23(1),A32(1),A31(1),A13(1)],[A12(2),A21(2),A23(2),A32(2),A31(2),A13(2)],col(1,:),'LineStyle','none'); 
%     elseif nr==2
%         ka=0.61; A1x=ka*A13+(1-ka)*A12;
%         fill([A1x(1),A23(1),A32(1),A31(1),A13(1)],[A1x(2),A23(2),A32(2),A31(2),A13(2)],col(1,:),'LineStyle','none'); 
%     end
    %fill([A10(1),A12(1),A13(1)],[A10(2),A12(2),A13(2)],'w','LineStyle','none'); 
    %fill([A20(1), A21(1),A23(1)],[A20(2), A21(2),A23(2)],'w','LineStyle','none');
    %fill([A30(1),A31(1),A32(1)], [A30(2),A31(2),A32(2)],'w','LineStyle','none');

    %plot([A12(1),A13(1)],[A12(2),A13(2)],'k--','LineWidth',2);
    %plot([A21(1),A23(1)],[A21(2),A23(2)],'k--','LineWidth',2);
    %plot([A31(1),A32(1)],[A31(2),A32(2)],'k--','LineWidth',2);
%end

%% Plotting the triangle
% plot([edges(:,1); edges(1,1)],[edges(:,2); edges(1,2)],'k','LineWidth',lw); 
% plot(edges(:,1),edges(:,2),'ko','MarkerSize',ms,'MarkerFaceColor','k','LineWidth',lw);  

%% Plotting the symbols for no cooperation 
% laS=0.86; rd=0.06; t=0:0.01:2*pi; red=[0.8 0 0]; msS=8; lwS=2;
% S1=laS*x1+(1-laS)/2*(x2+x3); S2=laS*x2+(1-laS)/2*(x1+x3); S3=laS*x3+(1-laS)/2*(x2+x1); 
% plot(S1(1),S1(2),'o','Color',red,'MarkerSize',msS,'LineWidth',lwS); 
% plot(S2(1),S2(2),'o','Color',red,'MarkerSize',msS,'LineWidth',lwS);
% plot(S3(1),S3(2),'o','Color',red,'MarkerSize',msS,'LineWidth',lwS);
% plot(mean(edges(:,1)),mean(edges(:,2)),'kx','MarkerSize',ms,'MarkerFaceColor','k');
% if nr>1, plot(mean(edges(:,1)),mean(edges(:,2)),'o','MarkerSize',msS,'LineWidth',lwS,...
%         'Color',col(3,:),'MarkerFaceColor','w'); 
%     plot(mean(edges(:,1)),mean(edges(:,2)),'x','MarkerSize',msS,'Color',col(3,:));
% end


%% Labels
% text(-1,2,lett,'FontSize',fsH,'FontName',fname,'FontWeight','bold');
% dyL=0.26;
% text(x1(1)+0.175,x1(2)-dyL,{'Full endowment','to player 3'},'FontSize',fsT,...
%     'FontName',fname,'HorizontalAlignment','center'); 
% text(x2(1)-0.175,x2(2)-dyL,{'Full endowment','to player 2'},'FontSize',fsT,...
%     'FontName',fname,'HorizontalAlignment','center'); 
% text(x3(1),x3(2)+dyL+0.035,{'Full endowment','to player 1'},'FontSize',fsT,...
%     'FontName',fname,'HorizontalAlignment','center'); 
%  
% text(max(edges(:,1))*0.85,1.1,{'Cooperation','feasible'},'FontSize',fsT,...
%         'FontName',fname,'HorizontalAlignment','center','Color',col(1,:));
% fill(max(edges(:,1))*0.85+[-0.1 0.1 0.1 -0.1],[0.73 0.73 0.87 0.87],col(1,:),'LineStyle','none'); 
% 
% axis off
end




function cl=getcolor(pi,nr, PiMin, PiMax); 
global col
%PiMin=1;

if(nr ~= 0) % I am not sure that it was a correct way!
    PiMin = PiMin + 1;
    PiMax = PiMax + 1;
    PiMax = 2;
    pi = pi + 1;
end

%if nr==1 
%    PiMax=1.75; 
%elseif nr==2
%    PiMax=2.4;  
%elseif nr==3
%    PiMax=2.4;
%end
% or we can use 3 instead of PiMax to compare better all examples!
x=(log(pi)-log(PiMin))^2/(log(PiMax)-log(PiMin))^2;
cl=x*col(5,:)+(1-x)*[1 1 1]; 
end