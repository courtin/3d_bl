function airplane = POC_Input()
%Inputs to airplane data structure.
%NOTE: Vehicle geometry is defined in initialize_geometry.m Harmonizing is
%on the TODO list. 
    insq2msq = 0.00064516;


airplane.weight.W = 36 * 4.45;
Sref = 17.95*.3048^2;
bref = 13*.3048;
bh = 4.4*.3048;
bv = 2.3*.3048;


airplane.geometry.Wing.Sref = Sref;
airplane.geometry.Wing.breg = bref;
airplane.geometry.Htail.S = 639*insq2msq;
airplane.geometry.Htail.Swet = 2.1*airplane.geometry.Htail.S;
airplane.geometry.Htail.x_c_m = .3;
airplane.geometry.Htail.t_c_avg = .12;
airplane.geometry.Htail.c_ma = 1.05*.3048;
airplane.geometry.Htail.b = bh;
airplane.geometry.Htail.AR = bh^2/airplane.geometry.Htail.S;

airplane.geometry.Vtail.S = 480*insq2msq;
airplane.geometry.Vtail.Swet = 2.1*airplane.geometry.Vtail.S;
airplane.geometry.Vtail.x_c_m = .3;
airplane.geometry.Vtail.t_c_avg = .12;
airplane.geometry.Vtail.c_ma = 1.29*.3048;
airplane.geometry.Vtail.b = bv;
airplane.geometry.Vtail.AR = bv^2/airplane.geometry.Vtail.S;

airplane.geometry.Fuse.l = 5.68*.3048;
airplane.geometry.Fuse.fr = 5.68/1.79;
airplane.geometry.Fuse.Swet = 2.5*airplane.geometry.Wing.Sref;

airplane.sim.flight_condition = 'std day';
end

