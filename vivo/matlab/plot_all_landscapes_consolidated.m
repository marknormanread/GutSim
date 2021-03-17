%
% matlab script for calculating bacterial relative abundance surfaces given protein and carb consumption data.
%
addpath('../landscape')

% reads in <rows, columns>
data = xlsread('mouse-data-dup-removed.xlsx',3);

prot = data(:,13);
carb = data(:,14);
% taxa data are already relative abundances. 
firmP = data(:,19);
bactP = data(:,20);

clostC = data(:,21);
bacilC = data(:,22);
erysiC = data(:,23);
deferC = data(:,24);
verruC = data(:,25);

lachnF = data(:,26);
clostF = data(:,27);
ruminF = data(:,28);
eubacF = data(:,29);
rikenF = data(:,30);
bacteF = data(:,31);


calculate_draw_landscape(prot,carb,firmP,'Firmicutes phylum, relative abundance (%)');
view(2);
print('landscapes/p_Firmicutes','-dpng');

calculate_draw_landscape(prot,carb,bactP,'Bacteroidetes phylum, relative abundance (%)');
view(2);
print('landscapes/p_Bacteriodetes','-dpng');

calculate_draw_landscape(prot,carb,clostC,'Clostridia class, relative abundance (%)');
view(2);
print('landscapes/c_Clostridia','-dpng');

calculate_draw_landscape(prot,carb,bacilC,'Bacillae class, relative abundance (%)');
view(2);
print('landscapes/c_Bacillae','-dpng');

calculate_draw_landscape(prot,carb,erysiC,'Erysipelotrichia class, relative abundance (%)');
view(2);
print('landscapes/c_Erysipelotrichia','-dpng');

calculate_draw_landscape(prot,carb,deferC,'Deferribacteres class, relative abundance (%)');
view(2);
print('landscapes/c_Deferribacteres','-dpng');

calculate_draw_landscape(prot,carb,verruC,'Verrucomicrobia class, relative abundance (%)');
view(2);
print('landscapes/c_Verrucomicrobia','-dpng');

calculate_draw_landscape(prot,carb,lachnF,'Lachnospiraceae family, relative abundance (%)');
view(2);
print('landscapes/f_Lachnospiraceae','-dpng');

calculate_draw_landscape(prot,carb,clostF,'Clostridiaceae family, relative abundance (%)');
view(2);
print('landscapes/f_Clostredeaceae','-dpng');

calculate_draw_landscape(prot,carb,ruminF,'Ruminococcaceae family, relative abundance (%)');
view(2);
print('landscapes/f_Ruminococcaceae','-dpng');

calculate_draw_landscape(prot,carb,eubacF,'Eubacteriaceae family, relative abundance (%)');
view(2);
print('landscapes/f_Eubacteriaceae','-dpng');

calculate_draw_landscape(prot,carb,rikenF,'Rikenellaceae family, relative abundance (%)');
view(2);
print('landscapes/f_Rikenellaceae','-dpng');

calculate_draw_landscape(prot,carb,bacteF,'Bacteroidaceae family, relative abundance (%)');
view(2);
print('landscapes/f_bacteroidaceae','-dpng');

