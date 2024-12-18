function  smadfitffl(name,model,fit)
% tic
% clear all
close all
arInit;

load gename.mat;% load the names 

ind=find(ismember(gename,name));

if model==1;
    arLoadModel('Simple_mode_SMAD_ffl1');%load the model
elseif model==2;
    arLoadModel('Simple_mode_SMAD_ffl2');%load the model
elseif model==3;
    arLoadModel('Simple_mode_SMAD_ffl3');%load the model
elseif model==4;
    arLoadModel('Simple_mode_SMAD_ffl4');%load the model
elseif model==5;
    arLoadModel('Simple_mode_SMAD_ffl5');%load the model
elseif model==6;
    arLoadModel('Simple_mode_SMAD_ffl6');%load the model
elseif model==7;
    arLoadModel('Simple_mode_SMAD_ffl7');%load the model
elseif model==8;
    arLoadModel('Simple_mode_SMAD_ffl8');%load the model
end

arLoadData(['data_' gename{ind}]);%load the rest of the data

ar.config.fiterrors = 0;%no error model
arCompileAll;

arSetPars(['ks'],1,1,1,-3,3);%load the rest of the data
arSetPars(['ls'],-1,1,1,-5,1);%load the rest of the data
arSetPars(['rs'],-3,1,1,-5,2);%load the rest of the data
arSetPars(['hs'],0.001,1,1,-5,1);%load the rest of the data
arSetPars(['ks1'],1,1,1,-3,3);%load the rest of the data
arSetPars(['kact'],1,1,1,-3,3);%load the rest of the data
arSetPars(['b0'],-1,1,1,-5,1);%load the rest of the data
arSetPars(['syn'],-1,1,1,-5,1);%load the rest of the data
arSetPars(['rs1'],-3,1,1,-5,2);%load the rest of the data
arSetPars(['rs2'],-3,1,1,-5,2);%load the rest of the data
arSetPars(['rp1'],-3,1,1,-5,2);%load the rest of the data
arSetPars(['hs1'],0.001,1,1,-5,1);%load the rest of the data
arSetPars(['kp'],1,1,1,-3,3);%load the rest of the data
arSetPars(['deg'],-1,1,1,-5,1);%load the rest of the data
arSetPars(['rp'],-3,1,1,-5,2);%load the rest of the data
arSetPars(['hp'],0.001,1,1,-5,1);%load the rest of the data
arSetPars('sd_gnex',0.11,0,0,-5,5);

tic

arFit
toc

arFitLHS(fit)

%print and plot the data
% arPrint
% arPlot
% arPlotter

%get the fitted values
dr=ar.model.data;
ksolu=dr.pNum;

kname=dr.p; 

time_step=[{'45_high'},{'90_high'},{'180_high'},{'360_high'},{'720_high'},{'1440_high'},{'45_low'},{'90_low'},{'180_low'},{'360_low'},{'720_low'},{'1440_low'}];%the time values that are going to be used as data. should be in column
time_data=[{'45_highdata'},{'90_highdata'},{'180_highdata'},{'360_highdata'},{'720_highdata'},{'1440_highdata'},{'45_lowdata'},{'90_lowdata'},{'180_lowdata'},{'360_lowdata'},{'720_lowdata'},{'1440_lowdata'}];%the time values that are going to be used as data. should be in column

tab_title=[time_step,kname,'chisquare_high','chisquare_low','pvalue',time_data];

prd=[];
dta=[];
chis=[];

for i=2:7;
    prd(i-1,:)=[dr(i).tExp,dr(i).yExpSimu];
    dta(i-1,:)=[dr(i).tExp,dr(i).yExp];
    chis(i-1,:)=dr(i).chi2;
end
chis=sum(chis,1);

dta=sortrows(dta,1);
prd=sortrows(prd,1);

val_data=(dta(:,2:3))';
pred_activ=(prd(:,2:3))';

if model==1
    string_name=['fit_ffl_model1_' gename{ind}];
elseif model==2
    string_name=['fit_ffl_model2_' gename{ind}];
elseif model==3
    string_name=['fit_ffl_model3_' gename{ind}];
elseif model==4
    string_name=['fit_ffl_model4_' gename{ind}];
elseif model==5
    string_name=['fit_ffl_model5_' gename{ind}];
elseif model==6
    string_name=['fit_ffl_model6_' gename{ind}];
elseif model==7
    string_name=['fit_ffl_model7_' gename{ind}];
elseif model==8
    string_name=['fit_ffl_model8_' gename{ind}];
end

arChi2Test
pvl=ar.pval;
pzeit=[10,19,37,73,145,288];
time_prot=[{'45_highprot'},{'90_highprot'},{'180_highprot'},{'360_highprot'},{'720_highprot'},{'1440_highprot'},{'45_lowprot'},{'90_lowprot'},{'180_lowprot'},{'360_lowprot'},{'720_lowprot'},{'1440_lowprot'}];%the time values that are going to be used as data. should be in column

d_prot=ar.model.condition.xExpSimu(pzeit,[1,3]);
    
valus=array2table([pred_activ(1,:),pred_activ(2,:),ksolu,chis,pvl,val_data(1,:),val_data(2,:),d_prot(:,1)',d_prot(:,2)'],'VariableNames',[tab_title,time_prot]);

geneID=gename(ind);
genenm=table(geneID);

valus=[genenm,valus];

eval([string_name '=valus;']);

%save the fitted values
if model==1
    save(['model_fitting/fit_ffl_model1_' gename{ind} '.mat'],['fit_ffl_model1_' gename{ind}])
    arSave(['ffl_model1_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==2
    save(['model_fitting/fit_ffl_model2_' gename{ind} '.mat'],['fit_ffl_model2_' gename{ind}])
    arSave(['ffl_model2_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==3
    save(['model_fitting/fit_ffl_model3_' gename{ind} '.mat'],['fit_ffl_model3_' gename{ind}])
    arSave(['ffl_model3_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==4
    save(['model_fitting/fit_ffl_model4_' gename{ind} '.mat'],['fit_ffl_model4_' gename{ind}])
    arSave(['ffl_model4_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==5
    save(['model_fitting/fit_ffl_model5_' gename{ind} '.mat'],['fit_ffl_model5_' gename{ind}])
    arSave(['ffl_model5_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==6
    save(['model_fitting/fit_ffl_model6_' gename{ind} '.mat'],['fit_ffl_model6_' gename{ind}])
    arSave(['ffl_model6_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==7
    save(['model_fitting/fit_ffl_model7_' gename{ind} '.mat'],['fit_ffl_model7_' gename{ind}])
    arSave(['ffl_model7_' gename{ind} 'lhs_' num2str(fit) '.mat'])
elseif model==8
    save(['model_fitting/fit_ffl_model8_' gename{ind} '.mat'],['fit_ffl_model8_' gename{ind}])
    arSave(['ffl_model8_' gename{ind} 'lhs_' num2str(fit) '.mat'])
end

end
