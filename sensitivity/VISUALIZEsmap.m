clear all
close all

sens = load('sensitivitymatrix.dat');
%load sensitivitymatrix
fidp = fopen('PARexp29frA02optim.m','r');
npar = 93;
for n = 1:npar
   par = fgetl(fidp);
   fgetl(fidp);
   if ~ischar(par), break, end;
   par = strtok(par,' ');
if (par(1) == 'V' | par(1)=='v' | par(1)=='I')
switch par
case 'Vfgly'
labels{n} = 'GPa + GPb';
case 'Vffpglm'
labels{n} = 'PGLM';
case 'Vbbpgi'
labels{n} = 'PGI';
case 'Vffpfk'
labels{n} = 'PFK';
case 'Vffald'
labels{n} = 'ALD';
case 'Vfftpi'
labels{n} = 'TPI';
case 'Vbbg3pdh'
labels{n} = 'G3PDH';
case 'Vffgad'
labels{n} = 'GAPDH';
case 'Vbbpgk'
labels{n} = 'PGK';
case 'Vffpgm'
labels{n} = 'PGM';
case 'Vffen'
labels{n} = 'EN';
case 'Vffpk'
labels{n} = 'PK';
case 'Vffldh'
labels{n} = 'LDH';
case 'VforCK'
labels{n} = 'CK';
case 'Vfadk'
labels{n} = 'ADK';
case 'I'
labels{n} = par;
    otherwise
        labels{n} = par;
end
elseif strcmp(par,'Katp_ATPase')
    labels{n} = par;
else
labels{n} = ' ';
end
end
fclose(fidp);

index = strmatch('VmaxATPase',labels,'exact');
labels(index)=[];
sens(index,:)=[];
index = strmatch('Katp_ATPase',labels,'exact');
labels(index)=[];
sens(index,:)=[];

sensrc = size(sens);

iptsetpref('ImshowAxesVisible', 'on')

%%
figure(14);clf;

% Row Normalized
for i = 1:sensrc(1)
    visualr(i,:) = mat2gray(sens(i,:)); 
end

axes('Position',[0.375,0.1,0.2,0.8])
title('Row Normalized')
imshow(visualr);
colormap(gray);
set(gca,'XTick',1:sensrc(2),'TickDir','out','TickLength',...
[0.005 0.00],'FontSize',4,'FontWeight','normal','XTickLabel',...
{'1';'2';'3';'4';'5';'6';'7';'8';'9'},'YTick',1:sensrc(1), 'YTickLabel',labels)

% Column Normalized
for i = 1:sensrc(2)
    visualc(:,i) = mat2gray(sens(:,i)); 
end

axes('Position',[0.5,0.1,0.2,0.8])
title('Column Normalized')
imshow(visualc);
colormap(gray);
set(gca,'XTick',1:sensrc(2),'TickDir','out','TickLength',...
[0.005 0.00],'FontSize',4,'FontWeight','normal','XTickLabel',...
{'1';'2';'3';'4';'5';'6';'7';'8';'9'},'YTick',1:sensrc(1), 'YTickLabel',labels)

%print -f3 -deps -tiff 'figure14b.eps'
%print -f2 -dtiff -r1200 'figure14.tiff'

fig14 = figure(14)
resp14 = fig2plotly(fig14, 'filename', 'fig14', 'strip', false)
