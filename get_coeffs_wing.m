function [cl,cx,cm] = get_coeffs_wing(a_w_deg,dCJ_B,flap_deg,airplane)
%Returns 2D cl,cx,cm based on fit to wind tunnel data

%for 10 degree motor mount (get_regression_coefficients(MAT,10,1,1,1), accurate as of 2/11/2019):
cl_coeffs =[0.4807    0.0801   -0.0001    0.0137   -0.0000;
           -0.0459    0.0392   -0.0000    0.0214   -0.0000];


cx_coeffs =[0.3430   -0.0086    0.0005   -0.0034    0.0000; 
           -1.1986    0.0136    0.0004    0.0158   -0.0001];

cm_coeffs = [-0.0346   -0.0009   -0.0000   -0.0031    0.0000;
             0.0419    0.0003    0.0000   -0.0037    0.0000];



[cl,cx,cm]=regression_results(a_w_deg, flap_deg, dCJ_B,cl_coeffs,cx_coeffs,cm_coeffs);
    