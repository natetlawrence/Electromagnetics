%--------------------------------------------------------------------------
%calculate scattering cross section with PW excitation
clear Q
Lambda=300:5:1000;
Radius=[1000];
epsilon={2 1};
isEz=1;

for ll=1:length(Lambda)
    [Q(:,ll)]=Cyl2D_ML_PW_Q_v12(Lambda(ll),Radius,epsilon,isEz);
end

h=figure;
plot(Lambda,Q(1,:),'linewidth',2)
xlabel('Wavelength','fontsize',18)
ylabel('Qscat','fontsize',18)
set(gca,'fontsize',18,'box','on')

%calculate scattering cross section with PW excitation
clear G
Lambda=300:5:1000;
Radius=[100];
epsilon={2 1};
source=[150 0];
SourceDir=3;

for ll=1:length(Lambda)
    [G(:,ll)]=Cyl2D_ML_PS_G_v12(Lambda(ll),Radius,epsilon,source,SourceDir);
end

h=figure;
hold on
plot(Lambda,(G(1,:)),'b','linewidth',2)
plot(Lambda,G(2,:),'r','linewidth',2)
plot(Lambda,G(3,:),'g','linewidth',2)
%plot(source_x-80,Q(:,3)/max(Q(:,3))*2,'--g','linewidth',2)
xlabel('Wavelength','fontsize',18)
ylabel('Gamma','fontsize',18)
legend({'G' 'Gr' 'Gnr'},'fontsize',18)
set(gca,'fontsize',18,'box','on')

Cyl2D_ML_PS_G_v12(Lambda(ll),Radius,epsilon,source,SourceDir);