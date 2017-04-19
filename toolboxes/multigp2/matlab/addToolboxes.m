function addToolboxes(os,add)


if os==0,
    %dir2 = '/home/anfelopera/Dropbox/JovenInvestigador2/codes/toolboxes/';
    dir2 = '/Users/usuarioutp/Dropbox/AndresFelipeLopezLopera/sdsim/codes/toolboxes/';
else
   dir2 = 'D:\Algoritmos_Doctorado\mlprojects2\mlprojects\';
end

if add
    func = str2func('addpath');
else
    func = str2func('rmpath');
end

% Updated in SVN yet
func(strcat(dir2,'netlab/NETLAB3p3'))
%func(strcat(dir2,'mocap/matlab'))
func(strcat(dir2,'mltools2/matlab'))
func(strcat(dir2,'kern2/matlab'))
func(strcat(dir2,'ndlutil2/matlab'))
func(strcat(dir2,'gp2/matlab'))
func(strcat(dir2,'optimi2/matlab'))
func(strcat(dir2,'datasets2/matlab'))
func(strcat(dir2,'switchedsim'))
func(strcat(dir2,'switchedlfm/matlab'))
%func(strcat(dir2,'drawing/matlab'))
func(strcat(dir2,'multigp2/matlab'))


