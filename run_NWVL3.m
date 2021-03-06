function [CL,CL1,CX,CM, ai, y] = run_NWVL3(alpha_in, dCJs, dFs, verbose, useTAT, blow_center, blowing_model)
%Solves the nonlinear VL system using Newton method.


%Set up geometry
%b = 13;
[vortex, cs, dCJs_v, dFs_v, N, b, Sref, vortex_cs, ra, rb] = initialize_geometry(verbose,dFs, dCJs );
y_fuse = 1.5;
%Load unblown airfoil polar
load bw02b_polar.mat cl cd cm alpha;
unblown.cl = cl;
unblown.cd = cd;
unblown.alpha = alpha;
unblown.cm = cm;

cl_unblown = cl;
cd_unblown = cd;
cm_unblown = cm;
alpha_unblown = alpha;

%Flight condition
alpha = alpha_in*pi/180;
beta = 0;
rho = 1;
Vinf = 1;
Ubar = [-(cos(alpha)*cos(beta));sin(beta);-sin(alpha)*cos(beta)];
dCJ = dCJs(4);
dF = dFs(2);


 if useTAT
     %TAT_model
     cla_b = 2*pi;
     cla_ub = 2*pi;
     a0_b = 0;
     a0_ub = 0;
     cl_nom_b = 2*pi*alpha;
 else
     [~,~,~,~,cl_nom_b] = POC_lin(alpha,dCJ, dF,verbose, blowing_model); %POC aero model
 end


%Locate the control points
for i = 1:N
    %Calculate the chord distribution (hard-coded for wing with a single
    %break (POC vehicle).
    y = vortex(2,i);
    if abs(y) <= y_fuse/2 && ~blow_center
        %c = c_root;
        dCJs_v(i) = 0; %No blowing across fuselage
        dFs_v(i) = 0;
    end
end
 dy = rb(2,:)-ra(2,:);

%Initialize the gamma distribution
y = vortex(2,:);
Gam0 = .5*vortex_cs(N/2)*cl_nom_b*Vinf;
Gam = Gam0.*(1-(2.*y./b).^2).^(.5);
Gam = Gam';
Gami = Gam;

%Solve Gamma using Newton method
err = 1e6;
max_iter = 100;
iter = 0;
thresh = 1e-3;
f = @(t,Gam) get_Gam_R(Gam, vortex, ra, rb, alpha,dCJs_v,dFs_v,N,vortex_cs,unblown, useTAT, blowing_model); 
fac = [];
while err > thresh && iter < max_iter
    r0 = f(1,Gam);
    [dfdy,fac]= numjac(f, 1, Gam, r0, 1e-3, fac,0);
    dGam = -inv(dfdy)*r0;
    Gam = Gam + dGam;
    err = sum(abs(r0));
    iter = iter+1;
end
if iter == max_iter && verbose
    fprintf(1, "Warning: Newton solver terminating on iteration limit\n")
    fprintf(1, "Last residual:\n")
    r0
    
end
%Solve Gamma using matlab build-in solver
% f = @(Gam) get_Gam_R(Gam, vortex, ra, rb, alpha,dCJs_v,dFs_v,N,vortex_cs,unblown, useTAT); 
% Gam = fsolve(f, Gam, optimoptions('fsolve','FunctionTolerance', .01));
% [R, Gam_new] = get_Gam_R(Gam, vortex, ra, rb, alpha,dCJs_v,dFs_v,N,vortex_cs,unblown, useTAT); 

Gamb = (Gam/Vinf)';

%Sanity check lift coefficient
vy = Ubar(1).*-(rb(2,:)-ra(2,:));
CL1 = 2*Gamb*(vy')/Sref;
%Calculate velocity distribution
Vbar = zeros(3, N);
Fi = zeros(3, N);
ls = zeros(3, N);
ai = zeros(1,N);
a_eff = zeros(1,N);

for i = 1:N
    for j = 1:N
        LBOUND = ~(i==j);% || i == j+1 || i == j-1);
        %LBOUND = 0;
        ri = .5*(ra(:,i)+rb(:,i));
        Vh = vorvel2(ri, ra(:,j), rb(:,j),LBOUND);
        Vbar(:,i) = Vbar(:,i) + Vh*Gamb(j);
    end
    Vbar(:,i) = Vbar(:,i)-Ubar;
    a_eff(i) = Vbar(3,i)/Vbar(1,i);
    ai(i) = -(alpha - a_eff(i));
end

T = [cos(alpha),0,sin(alpha);0,1,0;-sin(alpha),0,cos(alpha)];

%Solve for force and downwash distribution
V = Vbar*Vinf;
for i = 1:N
    l = (rb(:,i)-ra(:,i));
    %l(1) = 0;
    %l(3) = 0;
    ls(:,i) = l;
    Fi(:,i) = rho*Vinf*Gamb(i)*cross(V(:,i), l);
end

Fbari = Fi./(.5*rho*Vinf^2*Sref);

Fbar = zeros(3,1);

Fbar(1) = sum(Fbari(1,:));
Fbar(2) = sum(Fbari(2,:));
Fbar(3) = sum(Fbari(3,:));

ai2 = Fbari(1,:)./Fbar(3,:);
%Caluculate coefficients

Coeffs = T*Fbar;
Cdi = Coeffs(1);
CY = Coeffs(2);
CL = Coeffs(3);

%a_eff = alpha+ai;
cbar = Sref/b;
cbar = sum(vortex_cs.*dy)/b;
Scalc = sum(vortex_cs.*dy);


% figure()
% plot(vortex(2,:), ai2)
% hold on
% %plot(vortex(2,:), a_eff)
% title("Downwash Check")
% ylabel("\alpha_i")
% xlabel("Spanwise location (ft)")
cl_t = zeros(N,1);
cx_t = zeros(N,1);
cm_t = zeros(N,1);
cls = zeros(N,1);
cxs = zeros(N,1);
cms = zeros(N,1);

[Di, ai_T] = trefftz(Gamb, vortex, ra, rb, 0);
a_eff = alpha + ai_T;
if useTAT
    %TaT mode
    cl_t = 2*pi.*a_eff';
    cx_t = .0;
    cm_t = -.1;
else
    for i = 1:N
        if dCJs_v(i) == 0
            %Use lookup for bw02b airfoil
            [cl_t(i),cx_t(i)] = get_unblown_coeffs(a_eff(i)*180/pi, cl_unblown, cd_unblown, alpha_unblown, cm_unblown);
            cm_t(i) = -.1;%Placeholder cm
        else
            if blowing_model == 1
            %Use wind tunnel data
            [cl_t(i),cx_t(i),cm_t(i)] = get_coeffs_wing(a_eff(i)*180/pi,dCJs_v(i),dFs_v(i)*180/pi,1);
            elseif blowing_model == 2
            [cl_t(i),cx_t(i),cm_t(i)] = get_TATcl(a_eff(i),dCJs_v(i),dFs_v(i));
            else
                fprintf(1, "Invalid blown lift model")
            end
        end
    end
end

cls = cl_t.*cos(-ai_T') - cx_t.*sin(-ai_T');
t1 = cl_t.*sin(-ai_T');
t2 = cx_t.*cos(-ai_T');
cxs = cl_t.*sin(-ai_T') + cx_t.*cos(-ai_T');
cms = cm_t;
%CM = mean(cms.*vortex_cs'./cbar^2);

if verbose
    %Plot results
    fc = figure();
    plot(vortex(2,:), Gamb*2./(Vinf*cbar), 'x-')
    hold on 
    plot(vortex(2,:), Gamb*2./(Vinf)./vortex_cs)
    plot(vortex(2,:), cls)
    plot(vortex(2,:), Gami, '.-')
    plot(vortex(2,:), Gamb, '.-')
    title("Circulation")
    xlabel("Spanwise location (ft)")
    legend("2\Gamma / (c_{ref} V_\infty)", "c_l (from \Gamma)", "c_l (from \alpha_{eff})", "Initial circulation", "Final circulation")
    saveas(fc, "circ.pdf")
    
    fd = figure();
    %plot(vortex(2,:), ai, 'x-')
    hold on
    plot(vortex(2,:), ai_T.*180/pi)
    title("Downwash")
    ylabel("\alpha_i")
    xlabel("Spanwise location (ft)")
    %legend("Near-field", "Trefftz Plane")
    legend("Trefftz Plane")
    saveas(fd,"Downwash.pdf")
end
CL = trapz(y,cls.*vortex_cs')/Sref;
CXp  = trapz(y,cxs.*vortex_cs')/Sref;
CM  = trapz(y,cms.*vortex_cs')/(Sref*cbar);
Cl  = trapz(y,cls.*vortex_cs'.*y')/(Sref*b);
Cni  = trapz(y,cxs.*vortex_cs'.*y')/(Sref*b);
CX = CXp;
AR = b^2/Sref;
CDi2 = CL^2/(pi*AR);
if verbose
    disp('  ')
    disp('Initial Conditions')
    disp('==================')
    disp(['AoA      = ',num2str(alpha*180/pi),' deg'])
    disp(['AR       = ',num2str(AR)])
    disp(['Sref     = ',num2str(Sref/.3048^2),' m^2 (',num2str(Sref),' ft^2)'])
    disp(['dCJ      = ',num2str(dCJ)])
    disp(['dF       = ',num2str(dF*180/pi)])


    disp('  ')
    disp('Resulting coefficients')
    disp('==================')
    %disp(['CL (from F matrix)      = ',num2str(CL)])
    disp(['CL (from a_eff)         = ',num2str(CL)])
    disp(['CL (check, from gamma)  = ',num2str(CL1)])
    
    disp(['cl (nominal)            = ',num2str(cl_nom_b)])
    %disp(['CDi (from F matrix)     = ',num2str(Cdi)])
    disp(['CDi (trefftz))          = ',num2str(Di/(.5*rho*Vinf^2*Sref))])
    disp(['CDi (check,CL^2/(AR*pi))= ',num2str(CDi2)])
    disp(['CXp                     = ',num2str(CXp)])
    disp(['CXnetw                  = ',num2str(CX)])
    disp(['CM                      = ',num2str(CM)])
    disp(['CY                      = ',num2str(CY)])
    disp(['Cl (roll moment coeff.) = ',num2str(Cl)])
    disp(['Cni                     = ',num2str(Cni)])
    %disp(['cl_section (nominal)    = ',num2str(cl_section)])
    %disp(['e                        = ',num2str(e)])
end
ai = ai_T;
end
