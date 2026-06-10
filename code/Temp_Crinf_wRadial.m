function [T,xval,bestP,bestT,convo,minmet,structure]=Temp_Crinf_wRadial...
    (n,Tinmin,Tinmax,Pinmin,Pinmax,flow,Tset,metval,Cr) 

%{ 
This program is a modified version of Predict_Temp() to account for radial
heat loss. 

This program produces a temperature profile for a given reactor
configuration, maximum temperature, and set flow rate. It calls the
functions density(), heatcapacity(), thermalconductivity() and 
varyingpower(). 

Inputs:
    n: the index for given reactor congifuration, from 1-2^14-1, it is
    converted into a 14 digit binary number
    Tinmin: minimum inlet temperature 
    Tinmax: maximum inlet temperature 
    Pinmin: minimum total power
    Pinmax: maximum total power
    flow: set total flow rate in slpm 
    Tset: set maximum temperature 
    metval: maximum error metric 
    Cr: set ratio between high and low electrical conductivity sections 

Outputs:
    T: axial temperature profile in Celsius 
    xval: axial position array in m 
    bestP: total power input 
    bestT: inlet temperature 
    convo: CO2 conversion at outlet
    minmet: error metric of axial temperature profile
    strucure: array of 2's and 1's denoting reactor configuration 

%}
 
L=0.175; %m, length of reactor 
area=0.019^2*pi; %cross sectional area of reactor 

number=dec2bin(n); %binary number representing reactor configuration
structure=[]; %array of 1's and 0's representing high conductivity (1) and low conductivity (0) sections

for i=1:length(number)
    structure(i)=str2num(number(i)); %convert binary string entries to numbers
end

if length(number)<14 %if number does not have 14 digits, add 0's to the beginning 
    structure=[zeros(1,14-length(number)),structure];
end

structure=structure+1; %convert structure to array of 1's and 2's

%molar mass of reactants and products
M_CO2=44.09; %g/mol
M_H2=2.016; %g/mol
M_CO=28.01; %g/mol
M_H2O=18.0153; %g/mol
M_Ar=39.948; %g/mol 

%heat of reaction RWGS
Hrx=41.154; %kJ/mol 

%Constants 
P=101.325; %KPa atmosphere 
R=8.314*10^-3; %kJ/Kmol, gas constant 
CtoK=273.15;% convert celsius to kelvin 
Tinf=300;% K, temperature infinitely far away (for boundary conditions)

%Reaction Kinetics 
Eaf=150; %forward activation energy 
Af=9.673e8; % forward pre-exponential factor 
DelG0=37.5864; %kJ/mol, G0 in Gibbs free energy: DelG=G0+DelGp(T)
DelGp=-0.0349; %kJ/mol K
Ar=Af*exp(DelGp/R);% Reverse pre-exponential factor
Ear=-DelG0+Eaf; % reverse activation energy 
 


dS=L/length(structure); %m, length of single section of reactor 


nump=150;% number of points to iterate over
xval=linspace(0,L,nump); %m, axial position of reactor 
dx=xval(2)-xval(1); %m 


Q_stp=flow*0.001/60; %m^3/s, volumetric flow rate 
C=10713;% W/m, constant to define thermal conductivity: k=C/T
epsilon=0.3; %porosity
eta_R=1; %catalyst effectiveness factor 


Tdiff=100; %initiallizing variable for temeprature sweep
minmet=100; %initializing error value 
params=[]; %initializing to save parameters 



%coefficients of inlet and outlet boundary conditions
%{
h_x=A*exp(B(T-C)); 
%}

hincoef=[154.5 -0.0016 598.9; 134.1 -0.0027 592.6; 96.9 -0.0036 615.79;...
    79.7 -0.0047 575.8; 69.5 -0.0051 588.9]; 
houtcoef=[144.5 0.0048 540.5; 101.7 0.0019 547.4; 75.2 0.0015 536.8;...
    61.4 0.0002 526.3; 54.1 0.0002 573.68];

%U=1000; %radial loss value 
U_const=6.5*10^5; 



loop=1; %if it loops too many times for whatever reason, cut it off 


while Tdiff>10 && minmet>metval
    Tinvals=linspace(Tinmin,Tinmax,15);
    Pinvals=linspace(Pinmin,Pinmax,40);
    curve=1; %starting loop fresh 
    clear TCval hout Tkeep hout_calc params Tout

    for Tin=Tinvals
        for Pin=Pinvals

    
            fail=0; 

            h=hincoef(flow,1)*exp(hincoef(flow,2)*(Tin-hincoef(flow,3))); 
            pd=varyingpower(dS, Cr, area, xval, structure, Pin); %power dissipation 
            [C,U]=thermalconductivity(dS, xval, structure); 

            %end power density

            T(1)=Tin+CtoK; %now were in kelvin 

                rho=(1/5)*density("CO2",T(1))+(3/5)*density("H2",T(1))+...
                    (1/5)*density("Ar",T(1)); 

                cp=(1/5)*heatcapacity("CO2",T(1))+(3/5)*heatcapacity("H2",T(1))+...
                    (1/5)*heatcapacity("Ar",T(1));

                u=Q_stp/area*(T(1))/300;

                C_CO2(1)=(1/5)*P/(R*(T(1)));
                C_Ar(1)=(1/5)*P/(R*(T(1)));
                C_H2(1)=(3/5)*P/(R*(T(1)));
                C_CO(1)=0;
                C_H2O(1)=0; 

                Cond=[];
                Conv=[];
                Reac=[];
                Power=[]; 

 

                for i=2:nump-1 %i is the index of Tfit 

                    if i==2 %the left boundary 
            %the reaction rate 

                        rf=Af*C_CO2(i-1)*C_H2(i-1)*exp(-Eaf/(R*T(i-1))); %mol/m^3 s
                        rr=Ar*C_CO(i-1)*C_H2O(i-1)*exp(-Ear/(R*T(i-1))); %mol/m^3 s

                        r=-rf+rr; %consuming reactants is neg, consuming products is pos 
                        H=Hrx*1000+(heatcapacity("CO",T(i-1))*M_CO+heatcapacity("H2O",T(i-1))*M_H2O-...
                            heatcapacity("CO2",T(i-1))*M_CO2-heatcapacity("H2",T(i-1)))*(T(i-1)-Tinf); 

                       
                      


                        k=C(i-1)./(T(i-1));
                        U(i)=U(i)./(T(i-1)); %Radial heat loss term 
                        %units of each term in energy balacne is W/m^3 
                        T(i+1)=(h/k)*(T(i-1)-Tinf)*2*dx+T(i-1); %now we have i, i+1 ;

                        
                        T(i)=(-rho*cp*epsilon*u*(T(i+1)-T(i-1))/(2*dx)-...
                            pd(i-1)...
                            +eta_R*r*H...
                            +(k/dx^2)*(T(i-1)+T(i+1))...
                            -U(i)*Tinf)...
                            *1/(2*k/dx^2+U(i)); 
                        



                        C_CO2(i)=C_CO2(i-1)+(dx/u)*r-C_CO2(i-1)*(T(i)-T(i-1))/T(i-1); 
                        C_H2(i)=C_H2(i-1)+(dx/u)*r-C_H2(i-1)*(T(i)-T(i-1))/T(i-1); 
                        C_CO(i)=C_CO(i-1)-(dx/u)*r-C_CO(i-1)*(T(i)-T(i-1))/T(i-1); 
                        C_H2O(i)=C_H2O(i-1)-(dx/u)*r-C_H2O(i-1)*(T(i)-T(i-1))/T(i-1);
                        C_Ar(i)=C_Ar(i-1)-C_Ar(i-1)*(T(i)-T(i-1))/T(i-1);

                        rf=Af*C_CO2(i)*C_H2(i)*exp(-Eaf/(R*T(i)));
                        rr=Ar*C_CO(i)*C_H2O(i-1)*exp(-Ear/(R*T(i)));

                        r=-rf+rr; %consuming reactants is neg, consuming products is pos 

                        C_CO2(i+1)=C_CO2(i)+(dx/u)*r-C_CO2(i)*(T(i+1)-T(i))/T(i); 
                        C_H2(i+1)=C_H2(i)+(dx/u)*r-C_H2(i)*(T(i+1)-T(i))/T(i); 
                        C_CO(i+1)=C_CO(i)-(dx/u)*r-C_CO(i)*(T(i+1)-T(i))/T(i); 
                        C_H2O(i+1)=C_H2O(i)-(dx/u)*r-C_H2O(i)*(T(i+1)-T(i))/T(i-1);
                        C_Ar(i+1)=C_Ar(i)-C_Ar(i)*(T(i+1)-T(i))/T(i);

                        Cond=[Cond,k*(T(i+1)-2*T(i)+T(i-1))/(dx^2)];
                        Conv=[Conv,-rho*cp*epsilon*u*(T(i+1)-T(i-1))/(2*dx)];
                        Reac=[Reac,eta_R*r*H];
                        Power=[Power,pd(i)]; 

                        %so now we have T i-1, i, i+1, and C i-1, i, i+1 

    
                        Ctot=C_CO2(i+1)+C_H2(i+1)+C_CO(i+1)+C_H2O(i+1)+C_Ar(i+1); 


                        rho=(1/Ctot)*(C_CO2(i+1)*density("CO2",T(i+1))+C_H2(i+1)*density("H2",T(i+1))+...
                             C_CO(i+1)*density("CO",T(i+1))+C_H2O(i+1)*density("H2O",T(i+1))+...
                             C_Ar(i+1)*density("Ar",T(i+1))); 

                        cp=(1/Ctot)*(C_CO2(i+1)*heatcapacity("CO2",T(i+1))+C_H2(i+1)*heatcapacity("H2",T(i+1))+...
                            C_CO(i+1)*heatcapacity("CO",T(i+1))+C_H2O(i+1)*heatcapacity("H2O",T(i+1))+...
                            C_Ar(i+1)*heatcapacity("Ar",T(i+1))); 

                        u=Q_stp/area*(T(i+1))/Tinf;

                    else

                        rf=Af*C_CO2(i)*C_H2(i)*exp(-Eaf/(R*T(i))); %mol/m^3 s  
                        rr=Ar*C_CO(i)*C_H2O(i)*exp(-Ear/(R*T(i))); %mol/m^3 s  

                        r=-rf+rr; %consuming reactants is neg, consuming products is pos 
                        H=Hrx*1000+(heatcapacity("CO",T(i))*M_CO+heatcapacity("H2O",T(i))*M_H2O-...
                            heatcapacity("CO2",T(i))*M_CO2-heatcapacity("H2",T(i)))*(T(i)-Tinf); %J/mol 

                        coef=rho*cp*u*epsilon;
                        k=C(i)./(T(i));
                        U(i)=U(i)./(T(i));
                        T(i+1)=(-coef*T(i-1)/(2*dx)-k*T(i-1)/dx^2+k*2*T(i)/dx^2-...
                            pd(i-1)+eta_R*H*r+U(i)*(T(i)-Tinf))/(k/dx^2-coef/(2*dx));



                        C_CO2(i+1)=C_CO2(i)+(dx/u)*r-C_CO2(i)*(T(i+1)-T(i))/T(i); 
                        C_H2(i+1)=C_H2(i)+(dx/u)*r-C_H2(i)*(T(i+1)-T(i))/T(i); 
                        C_CO(i+1)=C_CO(i)-(dx/u)*r-C_CO(i)*(T(i+1)-T(i))/T(i); 
                        C_H2O(i+1)=C_H2O(i)-(dx/u)*r-C_H2O(i)*(T(i+1)-T(i))/T(i-1);
                        C_Ar(i+1)=C_Ar(i)-C_Ar(i)*(T(i+1)-T(i))/T(i); 

                        Cond=[Cond,k*(T(i+1)-2*T(i)+T(i-1))/(dx^2)];
                        Conv=[Conv,-rho*cp*epsilon*u*(T(i+1)-T(i-1))/(2*dx)];
                        Reac=[Reac,eta_R*r*H];
                        Power=[Power,pd(i)]; 


                        Ctot=C_CO2(i+1)+C_H2(i+1)+C_CO(i+1)+C_H2O(i+1)+C_Ar(i+1); 


                        rho=(1/Ctot)*(C_CO2(i+1)*density("CO2",T(i+1))+C_H2(i+1)*density("H2",T(i+1))+...
                            C_CO(i+1)*density("CO",T(i+1))+C_H2O(i+1)*density("H2O",T(i+1))+...
                            C_Ar(i+1)*density("Ar",T(i+1))); 

                        cp=(1/Ctot)*(C_CO2(i+1)*heatcapacity("CO2",T(i+1))+C_H2(i+1)*heatcapacity("H2",T(i+1))+...
                            C_CO(i+1)*heatcapacity("CO",T(i+1))+C_H2O(i+1)*heatcapacity("H2O",T(i+1))+...
                                C_Ar(i+1)*heatcapacity("Ar",T(i+1))); 

           
                       u=Q_stp/area*(T(i+1))/Tinf;

                    end
            
                    if T(i+1)<300 || T(i+1)>1000
                        fail=1; 
                        break
                    end 
                end %making temperature profile 

                if fail==0
                %now check it

                    Tkeep(curve,:)=T; 
                    Conversionfit(curve)=1-C_CO2(end)/C_Ar(end); 
                    hout(curve)=-k*(T(end)-T(end-2))/(2*dx)/(T(end)-Tinf);
                    hout_calc(curve)=houtcoef(flow,1)*exp(houtcoef(flow,2)*(T(end)-CtoK-houtcoef(flow,3)));              
                    Tout(curve)=max(T)-273.15; 
                    params(curve,:)=[Tin,Pin];
                    curve=curve+1;  
                    
                end
                clear T C_CO2 C_Ar C_CO C_H2O C_H2 Cond Conv Reac Power
               


        end %Pin loop
    end %Tin loop 

    if curve==1 %if nothing worked, 
        Tkeep(curve,:)=zeros(1,150);  
        convo=-1; 
        indmetric=1; 
        bestP=-1; 
        bestT=-1; 
        minmet=-1; 
        break %break the while loop 
    else %if it worked

    Deltah=hout-hout_calc;
    DeltaT=Tout-Tset; %pay attention to T out 



    met=abs(DeltaT)+abs(Deltah);
    indmetric=find(met==min(met));

    bestT=params(indmetric,1);
    bestP=params(indmetric,2);
    convo=Conversionfit(indmetric); 

    minmet=min(met); 

    Tmin_old=Tinmin; 
    Tmax_old=Tinmax;
    midT=(Tinmax+Tinmin)/2; 
    Tdiff=Tinmax-Tinmin; 

    if bestT<midT && bestT>Tinmin %bottom half
        Tinmax=midT;
    elseif bestT>midT && bestT<Tinmax
        Tinmin=midT; 
    elseif bestT==Tinmin
        Tinmin=Tmin_old-10;
        Tinmax=Tmin_old+10; 
    elseif bestT==Tinmax
        Tinmin=Tmax_old-10;
        Tinmax=Tmax_old+10;
    end

    end


    loop=loop+1; 
    
    if loop>5
        break
    end

end %end of while loop 

T=Tkeep(indmetric,:)-CtoK; 




end


%DENSITY 


function rho=density(gastype,T)
P=101.325e3; %Pa
R=8.314; %J/mol K

if gastype=="Ar"
    M=39.948; %g/mol 
elseif gastype=="H2"
    M=2.016; %g/mol
elseif gastype=="CO2"
    M=44.009; %g/mol
elseif gastype=="CO"
    M=28.01; %g/mol
else
    M=18.014; %g/mol 
end

rho=P*M/(R*T);

end

%HEAT CAPACITY

function cp=heatcapacity(gastype,T)
R=8.314; %J/mol K

if gastype=="Ar" 
   M=39.948; %g/mol 
   cp=5/2*R*1/M; 
elseif gastype=="H2"
    M=2.016; %g/mol
    alpha=R*3.057; 
    beta=R*2.677e-3; 
    gamma=R*-5.810e-6;
    delta=R*5.521e-9;
    epsilon=R*-1.812e-12;
    cp=(alpha+beta*T+gamma*T^2+delta*T^3+epsilon*T^4)/M; %J/(g K)
elseif gastype=="CO2"
    M=44.009; %g/mol
    alpha=R*2.401; 
    beta=R*8.735e-3; 
    gamma=R*-6.607e-6;
    delta=R*2.002e-9;
    epsilon=R*0e-12; 
    cp=(alpha+beta*T+gamma*T^2+delta*T^3+epsilon*T^4)/M; %J/(g K)
elseif gastype=="CO"
    M=28.01; %g/mol
    alpha=R*3.710; 
    beta=R*-1.619e-3; 
    gamma=R*3.692e-6;
    delta=R*-2.032e-9;
    epsilon=R*0.24e-12; 
    cp=(alpha+beta*T+gamma*T^2+delta*T^3+epsilon*T^4)/M; 
else %H2O
    M=18.014; 
    alpha=R*4.070; 
    beta=R*-1.108e-3; 
    gamma=R*4.125e-6;
    delta=R*-2.964e-9;
    epsilon=R*0.807e-12; 
    cp=(alpha+beta*T+gamma*T^2+delta*T^3+epsilon*T^4)/M;
end

end

%THERMAL CONDUCTIVITY
%function to calculate the effective thermal conductivity and appropriate
%radial loss term for infinite ratio structure. 

function [C,U]=thermalconductivity(dS, x, structure)
    C=NaN(1,length(x)); 
    U=NaN(1,length(x)); 

    for i=1:length(x) 
        if x(i)==x(end)
            index=length(structure);
        else
            index=floor(x(i)/dS)+1;
        end

        if structure(index)==1 %for infinity case, 1 and 2 are switched
            C(i)=10713;
            U(i)=6.8*10^5; 
        else
            C(i)=7650; 
            U(i)=6.8*10^5; 
        end
    end 
end


%VARYING POWER 

function P=varyingpower(dS, Cr, area, x, structure,Ptot)

bounds=0:dS:dS*length(structure); 
bounds=bounds-0.5*bounds(end); %shift so 0 is in the middle 
Lcoil=0.127; 
Rcoil=0.044; 

fun= @(x) ((x+Lcoil)./sqrt((x+Lcoil).^2+Rcoil^2)+...
    (Lcoil-x)./sqrt((Lcoil-x).^2+Rcoil^2)).^2; 

Bsum=0; 

for sec=1:length(structure)
    if structure(sec)==1
        Bsum=Bsum+integral(fun,bounds(sec),bounds(sec+1));
    else
        Bsum=Bsum+Cr*integral(fun,bounds(sec),bounds(sec+1)); 
    end
end
 
C=Ptot/(area*Bsum); 
 
xmid=x(end)/2; 
P=C*((x-xmid+Lcoil)./sqrt((x-xmid+Lcoil).^2+Rcoil^2)+...
    (Lcoil-x+xmid)./sqrt((Lcoil-x+xmid).^2+Rcoil^2)).^2; 

for i=1:length(x) 
    if x(i)==x(end)
        index=length(structure);
    else
        index=floor(x(i)/dS)+1;
    end

    if structure(index)==2
        P(i)=Cr*P(i); %scale the power of 2 structure
    end

end 

end

