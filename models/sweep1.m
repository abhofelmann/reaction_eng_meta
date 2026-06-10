clear all 
%{
This script runs each of the 2^14 possible reactor configurations. 
It calls the function Predict_Temp() and then saves the following
variables:
    Temp: array of modeled temeprature profile for each structure
    xval: axial position associated with Temp 
    bestP: array of total power dissipaiton for each structure
    bestT: array of inlet temperature for each structure
    convo: array of CO2 conversion for each structure
    minmet: array of solution metric for each structure
    struct: array of 1's and 2's representing each reactor configuration




%}
nvals=0:16383; %index of 2^14 structures  

flow=5; %set flow rate in slpm
Tmax=560; %maximum temperature withi reactor bed
Cr=2.34;  %AC impedance ratio between high and low conductivity sections


i=1; 
for n=nvals %iterate over each structure
    disp(i)
    Tinmin=100; %set min and max allowable temperature and power values
    Tinmax=400;
    Pinmin=90;
    Pinmax=150; 
    metval=10; %set large error metric
    [~,~,P,T,~,~,~]=...
        Predict_Temp(n,Tinmin,Tinmax,Pinmin,Pinmax,flow,Tmax,metval,Cr);
    Tinmin=T-10-2*rand(1); %refine min/max temperature and power 
    Tinmax=T+10+2*rand(1);
    Pinmin=P-10-2*rand(1);
    Pinmax=P+10+2*rand(1);
    metval=0.5; %set smaller metric 
    [Temp(i,:),xval,bestP(i),bestT(i),convo(i),minmet(i),struct(i,:)]=...
        Predict_Temp(n,Tinmin,Tinmax,Pinmin,Pinmax,flow,Tmax,metval,Cr);


    i=i+1; 
end

