function addToolboxes(os,add)


if os==0,
   %dir2 = '/home/mbaxgma4/mlprojects/';
   %dir2 = '/media/malvarez/89602034657BF6C4/LastCopyFromManchesterMachinesMadeOnMay2011/mbaxgma4/mlprojects/';
%    dir3 = '/home/malvarez/Documents/mlprojects/';
   dir3 = '../../toolboxes/';
else
   dir3 = '../../toolboxes/'; 
%    dir3 = 'D:\Algoritmos_Doctorado\mlprojects2\mlprojects\';
end

if add
    func = str2func('addpath');
else
    func = str2func('rmpath');
end

% Updated in SVN yet
func(strcat(dir3,'netlab/NETLAB3p3'))
func(strcat(dir3,'mltools2/matlab'))
func(strcat(dir3,'kern2/matlab'))
func(strcat(dir3,'ndlutil2/matlab'))
func(strcat(dir3,'gp2/matlab'))
func(strcat(dir3,'optimi2/matlab'))
func(strcat(dir3,'datasets2/matlab'))
func(strcat(dir3,'multigp2/matlab'))
func(strcat(dir3,'switchedlfm/matlab'))
func(strcat(dir3,'sim/matlab'))




