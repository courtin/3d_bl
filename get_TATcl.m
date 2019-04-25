function [cl, cx, cm] = get_TATcl(alfa,dCJ,dF)
%GET_TATCL Get section cl according to ideal TAT. 
%alfa and dF in radians
%Assumes dCJ = CJ (Vj >> Vinf)
    cLdf    = 2*sqrt(pi*dCJ)*sqrt(1+.151*sqrt(dCJ) + .139*dCJ);
    cLa     = 2*pi.*(1+.151*sqrt(dCJ) + .219*dCJ);
    
    cl = cLa*alfa + cLdf*dF;
    cx = -dCJ * cos(alfa+ dF); %Assumes VJ >> V_inf (need info on h/c to compute exact)
    cm = 0; %Not actual zero, but currently no TAT estimate. 
end

