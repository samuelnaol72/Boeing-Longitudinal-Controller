function [X_T, Z_T, X_u, X_w, Z_u, Z_w, M_u, M_w, M_wdot, M_q, X_delta_elev, Z_delta_elev, M_delta_elev, M_T] = Gen_dimensional_derivatives(stb_coef, gen_info)
    % Inputs: 
    %   Stability coefficients (stb_coef): 
    %       stb_coef = [CD, CD_u, CD_alpha, CL, CL_u, CL_alpha, CM_u, CM_alpha, CM_alpha_dot, CM_q, U0, CL_deltaelevator, CM_deltaelevator]
    %   General information about the aircraft (gen_info):
    %       gen_info = [M, Q, S, c, I_y]

    % Stability coefficients
    CD = stb_coef(1); CD_u = stb_coef(2);   CD_alpha = stb_coef(3);
    CL = stb_coef(4);CL_u = stb_coef(5);  CL_alpha = stb_coef(6);
    CM_u = stb_coef(7);  CM_alpha = stb_coef(8);CM_q = stb_coef(10);
    CM_alpha_dot = stb_coef(9);U0 = stb_coef(11);CL_delta_elevator = stb_coef(12);
    CM_delta_elevator = stb_coef(13);
    % General information about the aircraft
    M = gen_info(1); Q = gen_info(2); S = gen_info(3); c = gen_info(4); I_y = gen_info(5);

    % Calculate aerodynamic derivatives
    X_u = -(CD_u + 2 * CD) * Q * S / (M * U0);
    X_w = -(CD_alpha - CL) * Q * S / (M * U0);
    Z_u = -(CL_u + 2 * CL) * Q * S / (M * U0);
    Z_w = -(CL_alpha + CD) * Q * S / (M * U0);
    M_u = CM_u * Q * S * c / (I_y * U0);
    M_w = CM_alpha * Q * S * c / (I_y * U0);
    M_wdot = CM_alpha_dot * Q * S * c * c / (2 * I_y * U0 *U0);
    M_q = CM_q * Q * S * c * c / (2 * I_y * U0);

    X_delta_elev = 0; % This is usually assumed to be zero in practice.
    Z_delta_elev = -CL_delta_elevator * Q * S / M;
    M_delta_elev = CM_delta_elevator * Q * S * c / I_y;

    X_T=.3*32.2;Z_T=0;M_T=0;

end
