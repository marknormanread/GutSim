%
% Works on the phylogenetic-detailed-dup-removed.xlsx file
%
% matlab script for calculating bacterial relative abundance surfaces given protein and carb consumption data.
%

clear;
addpath('../landscape')

% reads in <rows, columns>
file = 'phylogenetic-dup-STD-removed.xlsx'
[numGenus,txtGenus,rawGenus]  = xlsread(file,'Bacteria Genus');
[numFamily, txtFamily, rawFamily] = xlsread(file,'Bacteria Family');
[numClass, txtClass, rawClass] = xlsread(file,'Bacteria Class');

prot = numGenus(1:end-1,13);	% end-1 because genus columns have row of averages, for which prot and carb have
carb = numGenus(1:end-1,14);	% no value. They acquire NaN. 
fat =  numGenus(1:end-1,15)
% taxa data are already relative abundances. 
% firmP = data(:,19);
% bactP = data(:,20);

% clostC = data(:,21);
% bacilC = data(:,22);
% erysiC = data(:,23);
% deferC = data(:,24);
% verruC = data(:,25);

% lachnF = data(:,26);
% clostF = data(:,27);
% ruminF = data(:,28);
% eubacF = data(:,29);
% rikenF = data(:,30);
% bacteF = data(:,31);

drawDots = true;
stiffness = 0.0005;	% low values mean stiffer sheet. 

% Plot genus' relative abundances
% genus information starts at column 20
for col=20:660-1					% -1 so total reads not plotted. 
	genus = numGenus(1:end-1,col);	% ignore the row of averages. 
	avg = mean(genus);
	name = char(txtGenus(1,col))
	index = findstr(name, ' (genus)');		% remove this string from the name
	if ~isempty(index)
		name = name(1:index-1)
	end
	if avg > 0.5		% less than this isn't worth plotting. 
		calculate_draw_landscape(prot,carb,fat,genus,strcat(name, ' genus, relative abundance (%)'), drawDots, stiffness);
		view(2);
		filename = strcat('genus/g_',name);
		print(filename, '-dpng');
	end
end

% Plot family's relative abundances
for col=20:218-1					% -1 so total reads not plotted.  
	family = numGenus(1:end-1,col);	% ignore the row of averages. 
	avg = mean(family);
	name = char(txtFamily(1,col))
	index = findstr(name, ' (family)');		% remove this string from the name
	if ~isempty(index)
		name = name(1:index-1)
	end	
	if avg > 0.5		% less than this isn't worth plotting. 
		calculate_draw_landscape(prot,carb,fat,family,strcat(name, ' family, relative abundance (%)'), drawDots, stiffness);
		view(2);
		filename = strcat('family/f_',name);
		print(filename, '-dpng');
	end
end

% Plot class's relative abundances
for col=20:59-1						% -1 so total reads not plotted. 
	cls = numClass(1:end-1,col);	% ignore the row of averages. 
	avg = mean(cls);
	name = char(txtClass(1,col))
	index = findstr(name, ' (class)');		% remove this string from the name
	if ~isempty(index)
		name = name(1:index-1)
	end
	if avg > 0.5		% less than this isn't worth plotting. 
		calculate_draw_landscape(prot,carb,fat,cls,strcat(name, ' class, relative abundance (%)'), drawDots, stiffness);
		view(2);
		filename = strcat('class/c_',name);
		print(filename, '-dpng');
	end
end
