%--------------------------------------------------------------------------
%Note

%This script is created by %Daniel Du from Proteomics and Metabolimics Core 
%Facility at MD Anderson Cancer Center. 

%This script is to show an example of using MIsoCor to correct fraction
%abundances of measured ions obtained by mass spectrometry.

%--------------------------------------------------------------------------

clear 
clc

%% perform correction
tracer = 'C'; %carbon tracer 
formula = 'C3H6O3'; %lactate
purity = 1; %purity is normally one if unspecified by the manufacturer
fam = [0.611 0.365 0.019 0.005]'; %fractional abundance of measured ion
[mdv,cm] = misocor(fam,formula,tracer,purity); %perform correction, one can also see how the correction matrix looks

%% plot out the results
h = bar([fam mdv]);
set(h,'edgecolor','w');
set(gca,'fontsize',12,'xtick',1:4,'xticklabel',{'C0','C1','C2','C3'});
xlabel('Lactate Isotopologue','fontsize',12); 
legend('FAM','MDV');
title('Example of MIsoCor','fontsize',14);

