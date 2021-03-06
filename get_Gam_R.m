function [R, Gam_new] = get_Gam_R(Gam, vortex, ra, rb, alpha,dCJs,dF,N,cs,unblown, useTAT, blowing_model)
%GET_GAM_R returns the residual, which is the difference between the
%assumed Gamma distribution Gam, and the new Gamma computed from that
%distribution. 
    [~, ai] = trefftz(Gam, vortex, ra, rb,0);
    a_eff = alpha + ai;

    cl_t = zeros(N,1);
    cx_t = zeros(N,1);
    cls = zeros(N,1);
    cxs = zeros(N,1);

if useTAT
    %TaT mode
    cl_t = 2*pi.*a_eff';
    cx_t = .0;
    cm_t = -.1;
else
    for i = 1:N
        if dCJs(i) == 0 && dF(i) == 0
            %Use lookup for bw02b airfoil
            [cl_t(i),cx_t(i), cm_t(i)] = get_unblown_coeffs(a_eff(i)*180/pi, unblown.cl, unblown.cd, unblown.alpha, unblown.cm);
            %cm_t(i) = -.1;%Placeholder cm
        else
            if blowing_model == 1
            %Use wind tunnel data
            [cl_t(i),cx_t(i),cm_t(i)] = get_coeffs_wing(a_eff(i)*180/pi,dCJs(i),dF(i)*180/pi,1);
            elseif blowing_model == 2
            [cl_t(i),cx_t(i),cm_t(i)] = get_TATcl(a_eff(i),dCJs(i),dF(i));
            else
                fprintf(1, "Invalid blown lift model")
            end
        end
    end
end

cls = cl_t.*cos(-ai') - cx_t.*sin(-ai');
cxs = cl_t.*sin(-ai') + cx_t.*cos(-ai');

Gam_new = .5*cs.*cls';
%Gam_new(1) = 0;
%Gam_new(end) = 0;
R = Gam - Gam_new';
end

