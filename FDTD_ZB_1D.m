function [Ex,epsr_end]=FDTD_ZB_1D(number_of_duration_time,...
    dfu,NP,rate_delta,TN0,varingN0,L2,varingCenter,R0)
% ==============================================================================
% Copyright (c) [2024] [Zhou Bo]
% All Rights Reserved
%
% This code is the intellectual property of [Your Full Name]. The code
% or any portion of it may be used exclusively for academic and scientific
% research purposes. Any other use, including but not limited to commercial
% application or adaptation, redistribution, or incorporation into other
% software, is strictly prohibited without prior written consent from the
% copyright owner.
%
% By using this code, you agree to abide by these terms and acknowledge that
% this code is provided "AS IS" without warranty of any kind, either expressed
% or implied, including but not limited to the implied warranties of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% ==============================================================================
%%
c0=2.99792458e8;%speed of light in free space
mu=4.0*pi*1.0e-7;%permeability of free space
eps=1.0/(c0*c0*mu);%permittivity of free space
eta=sqrt(mu/eps);
freq=1.0e+9;%frequency of source excitation
lambda=c0/freq;%wavelength of source excitation
omega=2*pi*freq;%angle frequency
T0=1/freq;t0=0.5*T0;tao=t0/(sqrt(2)*pi*0.8);
%%
ttl=2;%time step factor
zl=128;%space step factor
dz=lambda/zl;dt=dz/c0/ttl;
%%
rmax=exp(-700);
if R0==-0.985
    upml=1;
else
    upml=0;
end
orderbc=4;
delbc=upml*dz;
sigmam=-log(rmax)*(orderbc+1.0)/(2.0*eta*delbc); 
sigfactor=sigmam/(dz*(delbc^orderbc)*(orderbc+1.0));
kmax=1;
kfactor=(kmax-1.0)/dz/(orderbc+1.0)/delbc^orderbc;
%%

varingN=varingN0;%Length of time-varying region = N wavelengths

Nzmax=2*upml+(varingN+0)*zl;%Spatial steps

TN=TN0;%number of period of source

it_startN=TN;%start time of time-varing
it_duraN=number_of_duration_time;%duration time of time-varing
it_waitN=0;%waiting time
time_cycle=it_startN+it_duraN+it_waitN;% Number of simulation cycles
Ntmax=zl*ttl*time_cycle;%time step
it_TN=zl*ttl*(time_cycle);
omega_dt=omega*dt;Nzmax_1=Nzmax+1;Nzmax_d1=Nzmax-1;
wavelenth=Nzmax*4;Nx_SF_max=Nzmax-0.5*zl;Nx_SF_min=upml+4;%Nx_SF_min is also the location of the source
%%
Ex=zeros(Nzmax+1,1);Hy=zeros(Nzmax,1);Jix=zeros(Nzmax+1,1);
epsr=ones(Nzmax+1,Ntmax+1);CA=zeros(Nzmax+1,1);CB=zeros(Nzmax+1,1);CP=zeros(Nzmax+1,1);CQ=zeros(Nzmax+1,1);
ezinc=zeros(1,wavelenth);hyinc=zeros(1,wavelenth-1);

position=Nzmax*varingCenter+0*zl;
varing_list=1+(position-L2/2*zl):(position+L2/2*zl);

rand_sin_sec=rand_step(zl*ttl*it_duraN+1,zl*ttl*NP,dfu,rate_delta);
epsr_time_list=[1*ones(1,zl*ttl*it_startN),   rand_sin_sec,   rand_sin_sec(end)*ones(1,zl*ttl*it_waitN)];
epsr(varing_list,:)=repmat(epsr_time_list,length(varing_list),1);%with time-varing
% epsr(varing_list,:)=repmat(1*ones(1,Ntmax+1),length(varing_list),1);%no time-varying

%%
% figure(1)
% subplot(2,1,1),plot(Z,Ex(1:Nzmax),'r');
% % axis([0 Nzmax/zl*lambda -1 1]);
% ylabel('EZ');
% subplot(2,1,2),plot(Z,Hy,'b');
% xlabel('z (meters)');ylabel('HY');
% rect=get(gcf,'Position');
% rect(1:2)=[0 0];pic_num=1;
%%
fh=dt/mu/dz;fe=dt/dz;epsr=epsr*eps;
ca=1;cb=dt/eps/dz;da=1;db=dt/mu/dz;
for it=1:Ntmax
    %%
    for kk=1:upml
        x1=(upml-kk+1)*dz;
        x2=(upml-kk)*dz;
        sigma=sigfactor*(x1^(orderbc+1)-x2^(orderbc+1));
        ki=1+kfactor*(x1^(orderbc+1)-x2^(orderbc+1));
        facm=(2*epsr(upml+1,it)*ki-sigma*dt);
        facp=(2*epsr(upml+1,it)*ki+sigma*dt);
        CA(kk+1)=facm/facp;CA(Nzmax-kk+1)=facm/facp;
        CB(kk+1)=2*dt*eps/epsr(upml+1,it)/facp;CB(Nzmax-kk+1)=2*dt*eps/epsr(Nzmax-upml,it)/facp;
        CP(kk+1)=facm/facp;CP(Nzmax-kk+1)=facm/facp;
        CQ(kk+1)=2*dt*eps/mu/facp;CQ(Nzmax-kk+1)=2*dt*eps/mu/facp;
    end
    %%
    %Introducing unidirectional Gaussian waves with TSFS. Note that TSFS...
    % cannot be used at the same time as Jix later.
    if it<it_TN
        ezinc(1)=sin(omega_dt*it)*exp(-(it*dt-t0)^2/(tao)^2);
        % ez_inc(1)=sin(omega_dt*it);
    elseif it<=it_TN+1
        ezinc(1)=0;
    end
    hyinc(1:wavelenth-1)=da*hyinc(1:wavelenth-1)-db*(ezinc(2:wavelenth)-ezinc(1:wavelenth-1));
    eziabc=ezinc(wavelenth-1);
    ezinc(wavelenth)=(eziabc+ezinc(wavelenth))/2;
    ezinc(2:wavelenth-1)=ca*ezinc(2:wavelenth-1)-cb*(hyinc(2:wavelenth-1)-hyinc(1:wavelenth-2));
    %%
    Ex(upml+2:Nzmax-upml)=epsr(upml+2:Nzmax-upml,it)./epsr(upml+2:Nzmax-upml,it+1).*Ex(upml+2:Nzmax-upml)...
        -fe./epsr(upml+2:Nzmax-upml,it+1).*(Hy(upml+2:Nzmax-upml)-Hy(upml+1:Nzmax-upml-1)-dz*Jix(upml+2:Nzmax-upml));
    Ex(2:upml+1)=CA(2:upml+1).*Ex(2:upml+1)-CB(2:upml+1).*(Hy(2:upml+1)-Hy(1:upml))/dz;
    Ex(Nzmax-upml+1:Nzmax)=CA(Nzmax-upml+1:Nzmax).*Ex(Nzmax-upml+1:Nzmax)...
        -CB(Nzmax-upml+1:Nzmax).*(Hy(Nzmax-upml+1:Nzmax)-Hy(Nzmax-upml:Nzmax-1))/dz;
    %%
    %Introduce bidirectional Gaussian waves with Jix. Note that Jix cannot...
    % be used at the same time as the preceding TSFS.
    % if it<it_TN
    %     Jix(Nx_SF_min)=sin(omega_dt*it)*exp(-(it*dt-t0)^2/(tao)^2);
    % elseif it<=it_TN+1
    %     Jix(Nx_SF_min)=0;
    % end
    %%
    % TFSF left ex in total field   
    Ex(Nx_SF_min) = Ex(Nx_SF_min)- cb*hyinc(Nx_SF_min-1) ;
    % TFSF right ex in total field  
    % Ex(Nx_SF_max) = Ex(Nx_SF_max)+cb*hy_inc(Nx_SF_max);
    %%
    Ex(Nzmax_1)=0;
    Ex(1)=0;
    %%
    Hy(upml+1:Nzmax-upml)=Hy(upml+1:Nzmax-upml)-fh*(Ex(upml+2:Nzmax-upml+1)-Ex(upml+1:Nzmax-upml));
    Hy(1:upml)=CP(2:upml+1).*Hy(1:upml)-CQ(2:upml+1).*(Ex(2:upml+1)-Ex(1:upml))/dz;
    Hy(Nzmax-upml+1:Nzmax)=CP(Nzmax-upml+1:Nzmax).*Hy(Nzmax-upml+1:Nzmax)...
        -CQ(Nzmax-upml+1:Nzmax).*(Ex(Nzmax-upml+2:Nzmax+1)-Ex(Nzmax-upml+1:Nzmax))/dz;
    % TFSF left hy in scattered field
    Hy(Nx_SF_min-1) = Hy(Nx_SF_min-1)- db*ezinc(Nx_SF_min);
    % TFSF right hy in scattered field
    % Hy(Nx_SF_max) = Hy(Nx_SF_max)+ db*ez_inc(Nx_SF_max); 
    %%
    % if mod(it,50)==0
    %     figure(1)
    %     rtime=num2str(floor(it*dt/T0),2);
    %     subplot(2,1,1),plot(Z,Ex(1:Nzmax),'r');
    %     axis([0 Nzmax/zl*lambda -0.35 0.35]);
    %     grid on;
    %     yticks(-1:0.1:1);
    %     xticks(dz*(0:zl:Nzmax));xticklabels(0:1:Nzmax/zl);
    %     title(['time = ',rtime,' T0']);
    %     ylabel('EZ');
    %     subplot(2,1,2),plot(Z,Hy(1:Nzmax),'b');
    %     axis([0 Nzmax/zl*lambda -inf inf]);
    %     grid on;
    %     yticks(-1:0.1:1);
    %     xticks(dz*(0:zl:Nzmax));xticklabels(0:1:Nzmax/zl);
    %     title(['time = ',rtime,' ns']);
    %     xlabel('z (meters)');ylabel('HY');
    %     M=getframe(gcf,rect);
    %     I=frame2im(M);
    %     [I,map]=rgb2ind(I,256);
    %     if pic_num == 1
    %     imwrite(I,map,'test.gif','gif','Loopcount',inf,'DelayTime',0.01);
    %     else
    %     imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.01);
    %     end
    %     pic_num = pic_num + 1;
    %     figure(6)
    %     plot(epsr(:,it));ylim([eps,2*eps])
    % end 

end

epsr_end=epsr(:,end);

    function [sec]=rand_step(Nt,Tt,delta,deleta_rate)
        sec=ones(1,Nt);
        cyc1=floor(Nt/Tt);
        for iu=1:cyc1
            temp1=rand*delta;
            % temp1=randzb;
            rate1=temp1/floor(Tt*deleta_rate);
            temp2=0.2+rand*(1-2*deleta_rate-0.2);
            % temp2=0.2+randzb*(1-2*deleta_rate-0.2);
            Tt1=floor(Tt*(1-2*deleta_rate-temp2));
            Tt2=floor(Tt*(deleta_rate));
            Tt4=Tt2;
            Tt3=Tt-Tt1-Tt2-Tt4;
            temp3=ones(1,Tt2);
            for ip=2:Tt2
                temp3(ip)=temp3(ip-1)+rate1;
            end
            sec(1+(iu-1)*Tt:iu*Tt)=[ones(1,Tt1),temp3,(1+temp1)*ones(1,Tt3),temp3(end:-1:1)];
       end
    end
    function [randzb0]=randzb()
        randzb0=rand/2+rand/2;
    end
end
