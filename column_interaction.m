clc
clear
%{
X=input('Number of stripes: ');
for strip=1:X
    fprintf('Input for stripe %d\n\n', strip)
    Y(strip)=input('Number of rods in strip:\n');
    ag(strip)=input('area of a rod in strip:\n');
end
areas=Y.*ag;
a= input('Larger dimension of the column (inch):');
b= input('Smaller dimension of the column (inch):');
c= input('Reinforcement yield strength (ksi):');
x= input('28 days concrete cylinder strength (ksi):');
y= input('Concrete cover (inch):');
%}
X=4;
areas=[5*1.27 2*1.27 2*1.27 5*1.27];
c=60;
x=4;
y=2.5;
a=22;
b=22;

B1=0.85-0.05*((x*1000-4000)/1000);
if B1>0.85
    B1=0.85;
elseif B1<0.65
    B1=0.65;
else
    B1;
end

pector=[];
mector=[];

%Pure Compression Failure Point
pnc=0.85*(x)*a*b+c*sum(areas);
pector=[pector,pnc];
mector=[mector,0.0];

dis=((a-2*y)/(X-1));
areap=fliplr(areas);

%Balanced Failure
CB=87*(a-y)/(87+c);
AB=B1*CB;
CCC=0.85*x*AB*b;

NOCS=round(X*(CB-y)/(a-2*y));
NOTS=round(X*(a-CB-y)/(a-2*y));

for I=1:1:NOCS
    FSP(I)=87*(CB-y-dis*(I-1))/CB;
    if FSP(I)>c
        FSP(I)=c;
    elseif FSP(I)<-c
        FSP(I)=-c;
    else
        FSP(I);
    end

    COMP(I)=areas(I).*(FSP(I)-.85*x);
end

for IP=1:1:NOTS
    STRAIN(IP)=(.003/CB)*(a-CB-y-dis*(IP-1));
    FS(IP)=STRAIN(IP)*29000;
    
    if FS(IP)>c
        FS(IP)=c;
    elseif FS(IP)<-c
        FS(IP)=-c;
    else
        FS(IP);
    end
    
    TENS(IP)=FS(IP).*areap(IP);
end

PNB=ceil(CCC+sum(COMP)-sum(TENS));

for MOC=1:1:NOCS
    TENC(MOC)=COMP(MOC)*(a-y-y-dis*(MOC-1)-.5*a+y);
end

for MOT=1:1:NOTS
    TENT(MOT)=TENS(MOT)*(0.5*a-y-dis*(MOT-1));
end
MNB=ceil(CCC*(a-y-.5*AB-.5*a+y)+sum(TENC)+sum(TENT));
    
pector=[pector,PNB];
mector=[mector,MNB];

%Compression Zone Failure
pn1=PNB;
C1=CB;
%mn1=MNB;

while pn1>=PNB && pn1<=pnc 
    C1=C1+.001;
    no=(X*(C1-y)/(a-2*y));
    nosc=round(no);
    nost=round(X*(a-C1-y)/(a-2*y));
    
    if nosc<=0
        nosc=0;
    elseif nosc>X
        nosc=X;
    else
        nosc;
    end
    
    if nost<0
        nost=0;
    elseif nost>X
        nost=X;
    else
        nost;
    end
    
    %compression
    
    if nosc>=1
        for stc=1:1:nosc
            fsp(stc)=87*(C1-y-dis*(stc-1))/C1;
             if fsp(stc)>c
                fsp(stc)=c;
             elseif fsp(stc)<-c
                 fsp(stc)=-c;
             else
                 fsp(stc);
             end
        
            comp(stc)=areas(stc).*(fsp(stc)-.85*x);
            
         end
    else comp=0;
    end
    
    aaa=B1*C1;
    ccc=0.85*(x)*aaa*b;
    
    %tension
    if nost>=1
        for stt=1:1:nost
             strain(stt)=(0.003/C1)*(a-C1-y-dis*(stt-1));
             fs(stt)=strain(stt)*29000;
        
            if fs(stt)>c
                fs(stt)=c;
            elseif fs(stt)<-c
                fs(stt)=-c;
            else 
                fs(stt);
            end
        
            tens(stt)=fs(stt).*areap(stt);   
        end
        
    else tens=0;
    end
   
    pn1=ceil(ccc+sum(comp)-sum(tens))%final 
    
    if nosc>=1
        for moc=1:1:nosc
            momen(moc)=comp(moc)*(a-y-y-dis*(moc-1)-.5*a+y);
        end
    else momen=0;
    end
        
    if nost>=1
        for mot=1:1:nost
            moment(mot)=tens(mot)*(.5*a-y-dis*(mot-1));
        end
    else moment=0;
    end
        
    mn1=ceil(sum(momen)+sum(moment)+ccc*(a-y-.5*aaa-.5*a+y));
    
    pector=[pector,pn1];
    mector=[mector,mn1];
    
end


%Tension Zone Failure
pn2=PNB;
C2=CB;

while pn2>=0 && pn2<=PNB
    C2=C2-.001;
    no2=(X*(C2-y)/(a-2*y));
    nosc2=round(no2);
    nost2=round(X*(a-C2-y)/(a-2*y));
    
    if nosc2<=0
        nosc2=0;
    elseif nosc2>X
        nosc2=X;
    else
        nosc2;
    end
    
    if nost2<0
        nost2=0;
    elseif nost2>X
        nost2=X;
    else
        nost2;
    end
    
    %compression
    if nosc2>=1
        for stc2=1:1:nosc2
            fsp2(stc2)=87*(C2-y-dis*(stc2-1))/C2;
            if fsp2(stc2)>c
                fsp2(stc2)=c;
            elseif fsp2(stc2)<-c
                fsp2(stc2)=-c;
            else
                fsp2(stc2);
            end
        
        comp2(stc2)=areas(stc2).*(fsp2(stc2)-.85*x); 
        end
        
    else comp2=0;
    end
    aaa2=B1*C2;
    ccc2=0.85*(x)*aaa2*b;
    
    %tension
    if nost2>=1
        for stt2=1:1:nost2
            strain2(stt2)=(0.003/C2)*(a-C2-y-dis*(stt2-1));
        
            fs2(stt2)=strain2(stt2)*29000;
        
        if fs2(stt2)>c
            fs2(stt2)=c;
        elseif fs2(stt2)<=-c
            fs2(stt2)=-c;
        else
            fs2(stt2);
        end
        
        tens2(stt2)=fs2(stt2).*areap(stt2);
      
        end
    else tens2=0;
    end
    
    pn2=ceil(ccc2+sum(comp2)-sum(tens2))%final 
    
    if nosc2>=1
        for moc2=1:1:nosc2
        momen2(moc2)=comp2(moc2)*(a-y-y-dis*(moc2-1)-.5*a+y);
        end
    else momen2=0;
    end
    
    if nost2>=1
        for mot2=1:1:nost2
            moment2(mot2)=tens2(mot2)*(.5*a-y-dis*(mot2-1));
        end
    else moment2=0;
    end
        
    mn2=ceil(sum(momen2)+sum(moment2)+ccc2*(a-y-.5*aaa2-.5*a+y));
    
    pector=[pector,pn2];
    mector=[mector,mn2];
    
end

% Graph plotting
i=mector/12;
j=pector;

i(i<0)=NaN;
j(j<0)=NaN;

K=.65*i;
L=.65*j;
Z=.8*.65*pnc;
L(L>Z)=Z;

It=i';
Jt=j';

col=[It,Jt];

Kt=K';
Lt=L';

colp=[Kt,Lt];

length(col);

plot(i,j,'.')

hold on
plot(K,L,'r.')
grid on

legend('Nominal strength','ACI design strength')

xlabel('Mn, kips-ft')
ylabel('Pn, kips')
title('Column Strength Interaction Diagram')

hold off