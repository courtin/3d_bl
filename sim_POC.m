function [CL, CXt, CM, ai_t_avg, CL_t] = sim_POC(V, alpha, throttles, flaps, verbose, blow_center, blowing_model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sim_POC returns the aero force coefficients and tail flow condition for
%the 30% proof-of-concept vehicle.  It uses a vortex lattice method in
%combination with wind tunnel test data, QProp and XFOIL/empirical fits to
%generate the 3D corrects for a blown lift vehicle. 
%
%Dependencies:
%   Qprop must be installed at the location shown in the qprop_path
%   variable.  The motor and propfiles must be in the same folder as the
%   script, or at the path given by the motorfile and propfile variables. 
%
%Note: 
%   This code has the geometry for the MIT POC vehicle hard-coded. This
%   has eight motors and four flaps which can be independently controlled. The
%   method is easily generalizable to an AVL-type input file.  The geometry
%   can be changed by editing the initialize_geometry.m script; the assumption
%   of eight motors and four flaps is hard-coded in some areas. No geometric 
%   or control symmetry is assumed or required, however there is currently an
%   assumpution of a symmetric free-stream flight condition.  This is easy to
%   change, to incorporate quasi-steady rotations and/or sideslip. 
%
%Inputs:
%       V               [m/s] Freestream speed
%       alpha           [deg] Body-axis angle of attack
%       throttles       [-]   1x8 matrix giving the throttle setting of
%                               each motor, on a scale of 0 to 1. 
%       flaps           [rad] 1x4 matrix giving the deflection of
%                               each flap. 
%       vebose          [-]  flag that controls whether the results are
%                               printed. 1 = display, 0 = no output.
%       blow_center     [-]  flag that controls whether the center section
%                               is treated as a continuation of the blown 
%                               wing or an unblown, undeflected airfoil. 1
%                               = blown, 0 = no blowing. 
%       blowing_model   [-]  flag that controls what model is used to
%                           predict the 2D performance of blown wing
%                           sections. 
%                           1 = wind tunnel data. 
%                           2 = TAT model, per Maskell and Spence
%
%Outputs: 
%       CL              [-] Wing lift coefficient
%       CX              [-] Net streamwise force coefficient, including
%                           simple estimates for fuselage and empennage drag
%       CM              [-] Wing pitching moment
%       ai_t_avg        [rad] Average downwash at tail
%       CL_t            [-] Approximate tail lift coefficient, w/o elevator
%                           deflection
    if nargin == 0
        V = 5.5;         %m/s
        V = 10;         %m/s
%                      M1  M2  M3  M4  M5  M6  M7  M8                     
        throttles   =  [0   1   1   1   1   1   1   0]*.75;
%                      F1  F2  F3  F4 
        flaps       =  [40 40 40 40].*pi/180;

        alpha = 25;
        
        verbose = 1;
        blow_center = 0;
        blowing_model = 2;
    end

    useTAT = 0;
    
    
    v_oper = 22;    %Operating voltage
    hd_c = .5;
    Nmotor = 8;     %Assumptions that this is 8 are hard_coded into initialize_geometry.m 
    Nflaps = 4;     %Assumptions that this is 4 are hard_coded into initialize_geometry.m    
    
    motorfile = 'sc3014.txt';
    propfile  = 'cam9x6.txt';
    qprop_path = '~/tools/Qprop/bin/qprop';
    [diff_t, id, it] = unique(throttles);    %Figure out how many different throttle settings there are. 
    T = length(diff_t);
    
    CJ_diff     = zeros(1,T);                    %Get the dCJ of each unique motor
    dCJ_diff    = zeros(1,T);
    Ps_diff     = zeros(1,T);
    
    for i = 1:T
        if diff_t(i) > 0
            v_motor = diff_t(i)*v_oper;         %Operating voltage of each motor
            s = [qprop_path,' ',propfile, ' ', motorfile, ' ', num2str(V),' 0 ',num2str(v_motor) ' > result.out'];
            system(s);
            [~, RPM, T_N, Q_Nm, Pshaft_W, ...
            Volts, Amps, eta_mot, eta_prop, DV] = read_qprop_out('result.out');
            VJ = V+DV;
            CJ_diff(i) = hd_c*(VJ^2/V^2 + VJ/V);
            dCJ_diff(i) = hd_c*(VJ^2/V^2-1)*(V/VJ+1);
            Ps_diff(i) = Pshaft_W;
        else
            CJ_diff(i) = hd_c*2;
            dCJ_diff = 0;
        end
    end
    
    CJs     = zeros(1,Nmotor);
    dCJs    = zeros(1,Nmotor);
    Ps      = zeros(1,Nmotor);
    
    for i = 1:Nmotor
        CJs(i) = CJ_diff(it(i));
        dCJs(i) = dCJ_diff(it(i));
        Ps(i) = Ps_diff(it(i));
    end
    
    %Solve for the induced angles on the wing and resulting net forces
    [CL,~,CXw, CM, ai, y] = run_NWVL3(alpha, dCJs, flaps, verbose, useTAT, blow_center, blowing_model);
    
    %Add simple additional drag estimate for other vehicle components (skin
    %friction only)
    y = y.*.3048;
    rho = 1.225;
    airplane = POC_Input();
    Sref = airplane.geometry.Wing.Sref;
    W = airplane.weight.W;
    bh = airplane.geometry.Htail.b;
    
    [CDp, Drag_Decomp] = getProfileDrag(airplane, 0, 0.05, 'clean');
    fexcr = 1.1;
    CDp = CDp*fexcr;
    CL_stall = 2*W/(rho*V^2*Sref);
    CXt = CXw+CDp;
    disp(['Total shaft power       = ',num2str(sum(Ps)), ' W'])
    disp(['CDp                     = ',num2str(CDp)])
    disp(['CXtot                   = ',num2str(CXt)])
    disp(['CL_stall                = ',num2str(CL_stall)])
    
    %Look at downwash at the horizontal tail
    C = 0;
    for i = 1:length(ai)        %See how many data points are within the htail span
        if abs(y(i)) <= bh/2
            C = C+1;
        end
    end
    ai_tail = zeros(1,C);       %Calculate relevant values at those points
    y_tail  = zeros(1,C);
    c_tail  = zeros(1,C);
    k       = 1;
    for i = 1:length(ai)
        if abs(y(i)) <= bh/2
            y_tail(k) = y(i);
            ai_tail(k) = 2*ai(i);
            k = k+1;
        end
    end
    
    a_eff_t = ai_tail + alpha*pi/180;       %Assume htail uses the same unblown airfoil as the wing
    cl_t = zeros(1,C);
    a_t = alpha;
    load bw02b_polar.mat cl cd cm alpha;
    unblown.cl = cl;
    unblown.cd = cd;
    unblown.alpha = alpha;
    unblown.cm = cm;
    alpha = a_t;
    c_root = 1.36*.3048; %m
    c_tip = .62*.3048; %m
    L = length(y_tail);
    for i = 1:C                             %Get effective chord and cl distribution
        c_tail(i) = interp1([y_tail(1), y_tail(L/2),y_tail(L/2+1), y_tail(end)], [c_tip,c_root,c_root, c_tip], y_tail(i));
        [cl_t(i),~, ~] = get_unblown_coeffs(a_eff_t(i)*180/pi, unblown.cl, unblown.cd, unblown.alpha, unblown.cm);
    end
    CL_tb = trapz(y_tail,cl_t.*c_tail)/airplane.geometry.Htail.S;   %integrate for effective tail CL
    ai_t_avg = mean(ai_tail);
    CLt_max = 1.5;
    Sh = 639;
    S = 2500;
    lw = .16*12;
    lt = 6.5*12;
    cref = 1.41*12;
    de_max = 25;
    de_min = -25;
    cosa = cos(alpha*pi/180);
    sina = sin(alpha*pi/180);
    i_tail = 0*pi/180;
    M_net_w = CM + (CL*cosa+CXw*sina)*lw/cref;          %Net moment withough elevator deflection
    a_tail = alpha*pi/180 + ai_t_avg+i_tail;
    CL_tail_r = M_net_w/(Sh/S*lt/cref*cos(a_tail));
    CL_de = .871*S/Sh*.75;   %per radian
    ARh = 4.6;
    %CLat = 2*pi*(ARh/(2+ARh));
        
    de_req = (CL_tail_r - CL_tb)/CL_de;                 %Required elevator to trim
    if verbose
        figure()
        plot(y_tail/.3048, 2*ai_tail*180/pi);
        disp(['Max dCJ                 = ',num2str(max(dCJs))])
        disp(['ai_tail (avg)           = ',num2str(ai_t_avg*180/pi), ' deg'])
        disp(['CL_t (no elev. defl.)   = ',num2str(CL_tb)])
        disp(['Trim elev. deflection   = ',num2str(de_req*180/pi)])
        xlabel('Tail span (ft)')
        ylabel('Induced AoA at tail (deg)')
        title('Tail Downwash')
    end
end

