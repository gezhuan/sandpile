
clear all
close all
stepp=114
for i=stepp
    
    
    dfbf = dir('ColumPRE7_bf_0*h5');

    filenamebf = dfbf(i).name;
    num_str = regexp(filenamebf,'\d*\.?\d*','match');

    filenames = ['ColumPRE7_',num_str{2},'h5'];
   

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
    DR=R/2;
    DZ=R/2;
    DZC=3*R;
    DRC=3*R;
    DZC0=3*R;
   
    DRC0=3*R;
    Lz=7;
    Rbase=9.2;
    numC=0;
    for j=1:ceil(Lz/DZ)
        for k=1:ceil(Rbase/DR)
            numC=numC+1;
            Rp(1,numC)=2*R+k*DR;
            Rp(2,numC)=2*R+j*DZ;
        end
    end
    Centerx=11.249;
    Centery=11.249;
    for j=1:numC
        Vp(j)=0;
        nump(j)=0;
        if j==1
            Vs(j)=pi*DRC^2*DZC;
        end
        if j>1
            Vs(j)=(pi*(Rp(1,j)+DRC/2)^2-pi*(Rp(1,j)-DRC/2)^2)*DZC;
        end
        if j==1
            for k=1:length(xp)
                if PTag(k)>=-10
                    DisR=(xp(k)-Centerx)^2+(yp(k)-Centery)^2;
                    DisZ=zp(k);
                    if DisR<DRC0^2&&DisZ<DZC0
                        Vp(j)=4*pi/3*PR(k)^3+Vp(j);
                        nump(j)=nump(j)+1;

                    end
                end
            end
            phi(j)=Vp(j)/Vs(j);
        end
        if j>1
            for k=1:length(xp) 
                if PTag(k)>=-10
                    DisR=(xp(k)-Centerx)^2+(yp(k)-Centery)^2;
                    DisZ=zp(k);
                    if DisR>=(Rp(1,j)-DRC/2)^2&&DisR<(Rp(1,j)+DRC/2)^2&&DisZ<Rp(2,j)+DZC/2&&DisZ>=Rp(2,j)-DZC/2
                        Vp(j)=4*pi/3*PR(k)^3+Vp(j);
                        nump(j)=nump(j)+1;
                    end
                end
            end
            phi(j)=Vp(j)/Vs(j);
        end
    end   
end

save('sandpile','xp','yp','zp','PR');
for j=1:length(phi)
    RCR(j)=Rp(1,j);
    RCZ(j)=Rp(2,j);
end

scatter(RCR,RCZ,36, phi, 'filled')
%plot(Rp(1,:),Rp(2,:))

colorbar
