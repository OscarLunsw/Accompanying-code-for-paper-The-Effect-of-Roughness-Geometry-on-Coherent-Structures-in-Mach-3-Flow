%% Normalize velocity profiles

%%
Dia = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Diamonds\Boundary_profile.mat");
Dia_y = Dia.y;
Dia_u = Dia.u;
clear Dia
Dia = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Diamonds\DAQ_data\T625.mat");
Uinf = Dia.S.testConditions.trigger.Uinf_mps;
clear Dia
fact = Uinf/Dia_u(57);
Dia_u = Dia_u*fact;
%%
Squ = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Squares\Boundary_profile.mat");
Squ_y = Squ.y;
Squ_u = Squ.u;
clear Squ
Squ = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Squares\DAQ_data\T625.mat");
Uinf = Squ.S.testConditions.trigger.Uinf_mps;
clear Squ
fact = Uinf/Squ_u(57);
Squ_u = Squ_u*fact;

%%
Tri = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Tirangles\Boundary_profile.mat");
Tri_y = Tri.y;
Tri_u = Tri.u;
clear Tri
Tri = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Tirangles\DAQ_data\T625.mat");
Uinf = Tri.S.testConditions.trigger.Uinf_mps;
clear Tri
fact = Uinf/Tri_u(57);
Tri_u = Tri_u*fact;

%%
Smo = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Smooth\Boundary_profile.mat");
Smo_y = Smo.y;
Smo_u = Smo.u;
clear Smo
Smo = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Smooth\DAQ_data\T625.mat");
Uinf = Smo.S.testConditions.trigger.Uinf_mps;
clear Smo
fact = Uinf/Smo_u(57);
Smo_u = Smo_u*fact;

%%
Pin = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Pins\Boundary_profile.mat");
Pin_y = Pin.y;
Pin_u = Pin.u;
clear Pin
Pin = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Pins\DAQ_data\T625.mat");
Uinf = Pin.S.testConditions.trigger.Uinf_mps;
clear Pin
fact = Uinf/Pin_u(57);
Pin_u = Pin_u*fact;

%%
Pyr = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Pyramids\Boundary_profile.mat");
Pyr_y = Pyr.y;
Pyr_u = Pyr.u;
clear Pyr
Pyr = load("C:\Users\oscar\Documents\UNSW\Thesis\BOS\Final testing 23-07\Pyramids\DAQ_data\T625.mat");
Uinf = Pyr.S.testConditions.trigger.Uinf_mps;
clear Pyr
fact = Uinf/Pyr_u(57);
Pyr_u = Pyr_u*fact;

%% Plotting

plot(Pyr_y,Pyr_u,'-r' ...
    ,Pin_y,Pin_u,'-b' ...
    ,Smo_y,Smo_u,'-g' ...
    ,Tri_y,Tri_u,'-k' ...
    ,Squ_y,Squ_u,'-y' ...
    ,Dia_y,Dia_u,'-m')
legend("Pyramid","Pin","Smooth","Triangle","Square","Diamond")
set(gca,'Xscale','log')
ylabel('$u \, \mathrm{ (ms^{-1})}$','Interpreter','latex')
xlabel('$y \, \mathrm{(m)}$','Interpreter','latex')
xlim([0.0001,0.015])
%title("Mean Velocity Profiles Scaled to Freestream Velocity","FontName",'Times New Roman')
