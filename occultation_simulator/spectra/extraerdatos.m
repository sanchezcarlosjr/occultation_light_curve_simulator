%%Script para extraer datos espectrales de la estrella y multiplicarlos por
%%el filtro visible de TAOS-II.... se generan archivos .dat a partir de
%%llos archivos .h con dos columnas longitud de onda en amstrongs y flujo normalizado

lista=importdata('listah.txt');
%Fa=importdata('v_filter.h',',',6);
%Fb=importdata('v_filter.h',',',446);\
F=importdata('..\FiltCams\taos5564A.flt',' ',5);
indF=find(F.data(:,2)>0);

Fa.data=F.data(indF,1);
Fb.data=F.data(indF,2);

for k=1:length(lista)
    arch=char(lista(k));
    A=importdata(arch,',',7);
    B=importdata(arch,',',1905);
    figure(1)
    plot(A.data,B.data,Fa.data,Fb.data);
    ind=find(A.data>=min(Fa.data) & A.data<=max(Fa.data));
    res=Fb.data.*B.data(ind);
    figure(2)
    plot(Fa.data,res)
    dlmwrite(strcat(arch(1:length(arch)-2),'.dat'),[Fa.data res]);
end