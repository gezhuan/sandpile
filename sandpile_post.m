clear all
close all
stepp=2000
df         = dir('sandpile_*h5');


for i=stepp
    
    filenames  = df(i).name;
    num_str = regexp(filenames,'\d*\.?\d*','match');
    filenamebf = ['sandpile_bf_',num_str{1},'h5']; 

    % dfbf = dir('sandpile_bf_0*h5');
    % 
    % filenamebf = dfbf(i).name;
    % num_str = regexp(filenamebf,'\d*\.?\d*','match');
    % 
    % filenames = ['sandpile_',num_str{1},'h5'];
    % 

    Fnbf        = h5read(filenamebf,'/Normal');
    Ftbf        = h5read(filenamebf,'/Tangential');
    posbf       = h5read(filenamebf,'/Position');
    %Fbf         = h5read(filenamebf,'/Fn');
    ID1         = h5read(filenamebf,'/ID1');
    ID2         = h5read(filenamebf,'/ID2');
    Branch      = h5read(filenamebf,'/Branch');
    Angler      = h5read(filenames,'/PAngVelocity');
    M           = h5read(filenames,'/PAngMomentum');
    pos         = h5read(filenames,'/Position');
    
    PR          = h5read(filenames,'/Radius');   
    PTag        = h5read(filenames,'/PTag'); 
    FPT         = h5read(filenames,'/PUForce'); 
    vp          = h5read(filenames,'/PVelocity');
    EK          = h5read(filenames,'/PEkin');
    
    Fnx = Fnbf(1:3:end);
    Fny = Fnbf(2:3:end);
    Fnz = Fnbf(3:3:end);
    Ftx = Ftbf(1:3:end);
    Fty = Ftbf(2:3:end);
    Ftz = Ftbf(3:3:end);
    xp = pos(1:3:end);
    yp = pos(2:3:end);
    zp = pos(3:3:end); 

    Branchx = Branch(1:3:end);
    Branchy = Branch(2:3:end);
    Branchz = Branch(3:3:end);
    
    R=0.15;
    DR=2*R;
    DL=10*R; 
    Lz=7;
    Lx=9.2*2;
    Ly=9.2*2;
    numC=0;
    % 参数设置
    dr = 0.1;  % 间隔
    L = 1.0;   % 正方向和负方向的最大距离
    
    % 生成点
    x = -L:dr:L;
    
    % 如果你想确保 0 在正中间且点是对称的，可以用以下方式：
    N = floor(L / dr);
    x = (-N:N) * dr;
    
    % 显示结果
    disp(x);


    for j=1:ceil(Lz/DR)

        for k=-ceil(Ly/DR/2)*2:1:ceil(Ly/DR/2)*2
            for l=-ceil(Lx/DR/2)*2:1:ceil(Lx/DR/2)*2
                numC=numC+1;
                Rp(1,numC)=l*DR;
                Rp(2,numC)=k*DR;
                Rp(3,numC)=DL/2.0+j*DR;
            end
        end
    end



    id_bulk=find(PTag==-1);
    Xp=xp(id_bulk);
    Yp=yp(id_bulk);
    Zp=zp(id_bulk);
    id_bulk2=find(PTag==-15);
    ZBp=zp(id_bulk2);
    mean(Xp)
    %Centerx=11.249;
    %Centery=11.249;
    Centerx=0;
    Centery=0;
    Centerz=mean(ZBp)+0.1;
    %Rp(1,numC+1)=Centerx;
    %Rp(2,numC+1)=Centerz;

    for j=1:numC
        Vp(j)=0;
        nump(j)=0;
        Vs(j)=DL^3;
        for k=1:length(xp) 
            if PTag(k)==-1
                if abs(xp(k)-Rp(1,j))<=DL/2&&abs(yp(k)-Rp(2,j))<=DL/2&&abs(zp(k)-Rp(3,j))<=DL/2
                    Vp(j)=4*pi/3*PR(k)^3+Vp(j);
                    nump(j)=nump(j)+1;
                end
            end
        end

        Fcn=zeros(3);
        NUMContact(j)=0;
        for k = 1:length(ID2)
            Disx1=abs(xp(ID1(k)+1)-Rp(1,j));
            Disy1=abs(yp(ID1(k)+1)-Rp(2,j));
            Disz1=abs(zp(ID1(k)+1)-Rp(3,j));
            Disx2=abs(xp(ID2(k)+1)-Rp(1,j));
            Disy2=abs(yp(ID2(k)+1)-Rp(2,j));
            Disz2=abs(zp(ID2(k)+1)-Rp(3,j));
            testID1=ID1(k)+1;
            testID2=ID2(k)+1;
            if Disx1<=DL/2&&Disy1<=DL/2&&Disz1<=DL/2&&...
                    Disx2<=DL/2&&Disy2<=DL/2&&Disz2<=DL/2
                NUMContact(j)=NUMContact(j)+1;

                bfvector=[Branchx(k),Branchy(k),Branchz(k)];
                Fcn= ([Fnx(k);Fny(k);Fnz(k)]+[Ftx(k);Fty(k);Ftz(k)])*bfvector+Fcn;
            end
        end

        ES=Fcn/Vs(j);
        P(j)=abs(1/3*(ES(1,1)+ES(2,2)+ES(3,3)));
        tau(j)=sqrt(1/9*((ES(1,1)-P(j))^2+(ES(2,2)-P(j))^2+(ES(3,3)-P(j))^2+ES(1,2)^2+ES(1,3)^2+ES(2,3)^2 ...
        +ES(2,1)^2+ES(3,1)^2+ES(2,3)^2));
        mu(j)=abs(tau(j))/P(j);
        phi(j)=Vp(j)/Vs(j);

    end   
end

save('sandpile');

%%
clear all 
close all 
clc
load('sandpile');
nums=0;

ms=8;
lw=1.2;
numslayer1=0;
numslayer2=0;
numslayer3=0;
zlayer1=DL/2.0+DR*4.0;
zlayer2=DL/2.0+DR*8.0;
zlayer3=DL/2.0+DR*12.0;
for j=1:length(phi)
    if phi(j)>=0.2&&abs(Rp(2,j)-0.0)<=0.0001
        nums=nums+1;
        id_bulkp(nums)=j;
        RCX(nums)=Rp(1,j);
        RCZ(nums)=Rp(3,j);
        postphi(nums)=phi(j);
        postVs(nums)=Vs(j);
        postP(nums)=P(j);
        postmu(nums)=mu(j);
        postZc(nums)=NUMContact(j)/nump(j);
        if abs(Rp(3,j)-zlayer1)<=0.0001
            numslayer1=numslayer1+1;
            RCXlayer1(numslayer1)=Rp(1,j);
            Player1(numslayer1)=P(j);
            philayer1(numslayer1)=phi(j);
            Zclayer1(numslayer1)=NUMContact(j)/nump(j);
        end
        if abs(Rp(3,j)-zlayer2)<=0.0001
            numslayer2=numslayer2+1;
            RCXlayer2(numslayer2)=Rp(1,j);
            Player2(numslayer2)=P(j);
            philayer2(numslayer2)=phi(j);
            Zclayer2(numslayer2)=NUMContact(j)/nump(j);
        end
        if abs(Rp(3,j)-zlayer3)<=0.0001
            numslayer3=numslayer3+1;
            RCXlayer3(numslayer3)=Rp(1,j);
            Player3(numslayer3)=P(j);
            philayer3(numslayer3)=phi(j);
            Zclayer3(numslayer3)=NUMContact(j)/nump(j);
        end
    end
end

figure(1)
scatter(RCX,RCZ,36, postphi, 'filled' ,'Marker', 's')
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
title('$\phi$ Distribution','Interpreter','latex');
%plot(Rp(1,:),Rp(2,:))
caxis([0.5 0.65])
colorbar

figure(2)
scatter(RCX,RCZ,36, postmu, 'filled' ,'Marker', 's')
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
title('$\mu$ Distribution','Interpreter','latex');
%plot(Rp(1,:),Rp(2,:))
%caxis([0.45 0.7])
colorbar

figure(3)
scatter(RCX,RCZ,36, postP, 'filled' ,'Marker', 's')
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
title('$P$ Distribution','Interpreter','latex');
%plot(Rp(1,:),Rp(2,:))
%caxis([0.45 0.7])
colorbar

figure(4)
plot(RCXlayer1,Player1,'o','markersize',ms,'color','k','markerface','r','Linewidth',lw)
hold on
plot(RCXlayer2,Player2,'o','markersize',ms,'color','k','markerface','g','Linewidth',lw)
hold on
plot(RCXlayer3,Player3,'o','markersize',ms,'color','k','markerface','b','Linewidth',lw)
hold on
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('P', 'FontSize', 20);
legend('layer1 $z=2.85$','layer2 $z=3.6$','layer3 $z=4.35$','interpreter','latex','Location','southeast');

figure(5)
plot(RCXlayer1,philayer1,'o','markersize',ms,'color','k','markerface','r','Linewidth',lw)
hold on
plot(RCXlayer2,philayer2,'o','markersize',ms,'color','k','markerface','g','Linewidth',lw)
hold on
plot(RCXlayer3,philayer3,'o','markersize',ms,'color','k','markerface','b','Linewidth',lw)
hold on
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('\phi', 'FontSize', 20);
figure(6)
plot(RCXlayer1,Zclayer1,'o','markersize',ms,'color','k','markerface','r','Linewidth',lw)
hold on
plot(RCXlayer2,Zclayer2,'o','markersize',ms,'color','k','markerface','g','Linewidth',lw)
hold on
plot(RCXlayer3,Zclayer3,'o','markersize',ms,'color','k','markerface','b','Linewidth',lw)
hold on
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('Zc', 'FontSize', 20);
figure(7)
scatter(RCX,RCZ,36, postZc, 'filled' ,'Marker', 's')
set(gca,'FontSize',27,'FontName','Times New Roman','LineWidth',2);
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
title('$Z_{c}$ Distribution','Interpreter','latex');
%plot(Rp(1,:),Rp(2,:))
%caxis([0.45 0.7])
colorbar