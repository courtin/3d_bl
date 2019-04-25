function [cla, cl0, a0, cl1] = get_blown_cla(alfa, dCJ, dF, cl_u, cd_u, alpha_u, cm_u, blowing_model)
    %Calculate the cla of a blown wing, based on local linearization of WT
    %fit. If there is no blowing revert to unblown model.
    dalfa = -1 * pi/180;
    if dCJ == 0
        %Use lookup for bw02b airfoil
        [cl1,~] = get_unblown_coeffs(alfa*180/pi, cl_u, cd_u, alpha_u, cm_u);
        [cl2,~] = get_unblown_coeffs((alfa+dalfa)*180/pi, cl_u, cd_u, alpha_u, cm_u);
        cla = (cl2-cl1)/(dalfa);
        %cla = 2*pi;
        cl0 = cl1-cla*alfa;
    else
        if blowing_model == 1
            %Use wind tunnel data
            [cl1,~,~] = get_coeffs_wing(alfa*180/pi,dCJ,dF*180/pi,1);

            [cl2,~,~] = get_coeffs_wing((alfa+dalfa)*180/pi,dCJ,dF*180/pi,1);
        elseif blowing_model == 2
            %Use M&S 2D TAT Model
            [cl1, ~, ~] = get_TATcl(alfa,dCJ,dF);
            [cl2, ~, ~] = get_TATcl(alfa+dalfa,dCJ,dF);
        else
            fprintf(1,"Invalid blown lift model\n")
        end
        cla = (cl2-cl1)/(dalfa);
        cl0 = cl1-cla*alfa;
    end
    
    a0 = -cl0/cla;
end

