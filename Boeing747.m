
U0= 870.767717; Iy= .331e8; m= 636636/32.2;
X_u= (-3.954e2)/m; X_w= (3.144e2)/m; Xde= (1.544e4)/m;
Z_u= (-8.383e2)/m; Z_w= (-7.928e3)/m; Z_q= (-1.327e5)/m; Z_wd= (1.214e2)/m; Zde=(-3.677e5)/m;
M_u= (-2.062e3)/Iy; M_w= (-6.289e4)/Iy; M_q= (-1.327e7)/Iy; M_wd=(-5.296e3)/Iy; Mde=(-4.038e7)/Iy;
X_T=.3*32.2;

%A = [X_u, X_w, 0, -32.2;
%    Z_u, Z_w, U0, 0;
%   M_u + M_wd * Z_u, M_w + M_wd * Z_w, M_q + M_wd * U0, 0;
%  0, 0, 1, 0]
%  This is only for boeing at the 40,000ft
%B = [Xde X_T; Zde 0; Mde + M_wd * Zde 0; 0 0] 

% the following is more accurate
A = [X_u, X_w, 0, -32.2;
     (Z_u)/(1-Z_wd), (Z_w)/(1-Z_wd), (U0+Z_q)/(1-Z_wd), 0;
     M_u +  (M_wd/(1-Z_wd) )* Z_u,  M_w + (M_wd/(1-Z_wd)) * Z_w, M_q + (M_wd/(1-Z_wd)) *(U0+Z_q), 0;
     0, 0, 1, 0];
B = [Xde X_T; Zde 0; Mde + M_wd * Zde 0; 0 0];
% Find eigenvalues and eigenvectors
[vecs, vals] = eig(A);

% Extract eigenvalues from the diagonal matrix
eigenvalues = diag(vals);

% Identify short and long period modes
[~, max_real_index] = max(abs(real(eigenvalues)));
short_period_mode_index = [max_real_index, max_real_index + 1];
long_period_mode_index = setdiff(1:length(eigenvalues), short_period_mode_index);

% Display original eigenvalues and eigenvectors
%disp('Original Eigenvalues:');
%disp(eigenvalues);
%disp('Original Eigenvectors:');
%disp(vecs);

% Display short period mode eigenvalues and eigenvectors
disp('Short Period Mode Eigenvalues:');
disp(eigenvalues(short_period_mode_index));
disp('Short Period Mode Eigenvectors:');
disp(vecs(:, short_period_mode_index));

% Display long period mode eigenvalues and eigenvectors
disp('Long Period Mode Eigenvalues:');
disp(eigenvalues(long_period_mode_index));
disp('Long Period Mode Eigenvectors:');
disp(vecs(:, long_period_mode_index));

% Construct the state-space system
sys_exact = ss(A, B, eye(4), 0);
% Define transfer functions for each state variable using elevator input
[num_u, den_u] = tfdata(minreal(tf(sys_exact(1, 1))), 'v');
[num_w, den_w] = tfdata(minreal(tf(sys_exact(2, 1))), 'v');
[num_q, den_q] = tfdata(minreal(tf(sys_exact(3, 1))), 'v');
[num_theta, den_theta] = tfdata(minreal(tf(sys_exact(4, 1))), 'v');
% Display the exact fourth-order transfer functions for each state variable
TF_u = tf(num_u, den_u);
TF_alpha = tf(num_w/U0, den_w);
TF_q = tf(num_q, den_q);
TF_theta = tf(num_theta, den_theta);
TF_h= tf(U0*num_theta-num_w, conv([1 0], den_theta));

% Define transfer functions for each state variable using throttle input
[num_ut, den_ut] = tfdata(minreal(tf(sys_exact(1, 2))), 'v');
[num_wt, den_wt] = tfdata(minreal(tf(sys_exact(2, 2))), 'v');
[num_qt, den_qt] = tfdata(minreal(tf(sys_exact(3, 2))), 'v');
[num_thetat, den_thetat] = tfdata(minreal(tf(sys_exact(4, 2))), 'v');
TF_ut = tf(num_ut, den_ut);  % transfer function from speed to thrust input

% Short period approximation
A_sp = [Z_w, 1; M_w * U0 + Z_w * U0 * M_wd,  M_q + M_wd * U0];
B_sp = [Zde / U0; Mde + Zde * M_wd];

% Construct the state-space system
sys_sp = ss(A_sp, B_sp, eye(2), 0);

% Convert the state-space system to transfer functions
[num_alpha_sp, den_alpha_sp] = tfdata(minreal(tf(sys_sp(1, 1))), 'v');
[num_q_sp, den_q_sp] = tfdata(minreal(tf(sys_sp(2, 1))), 'v');
TF_alpha_sp = tf(num_alpha_sp, den_alpha_sp);
TF_q_sp = tf(num_q_sp, den_q_sp);


% Phugoid approximation
%A_p = [X_u, -32.2; -Z_u / U0,  0];
%B_p = [Xde - Mde * (X_w/M_w); (-Zde + Mde * (Z_w/M_w)) / U0];

A_p = [X_u + (X_w)*(m*U0*M_u*Iy-Z_u*m*M_q*Iy)/ (Z_w*m*M_q*Iy-m*U0*M_w*Iy), -32.2; (Z_u*m*M_w*Iy-Z_w*m*M_u*Iy)/(Z_w*m*M_q*Iy-m*U0*M_w*Iy),  0]; % to be more accurate i dont want to use course approximation
B_p = [Xde + (X_w)*(m*U0*Mde*Iy-Zde*m*M_q*Iy)/(Z_w*m*M_q*Iy-m*U0*M_w*Iy) X_T; (Zde*m*M_w*Iy-Z_w*m*Mde*Iy)/(Z_w*m*M_q*Iy-m*U0*M_w*Iy) 0];
% Construct the state-space system
sys_p = ss(A_p, B_p, eye(2), 0);

% Convert the state-space system to transfer functions for elevator input
[num_u_p, den_u_p] = tfdata(minreal(tf(sys_p(1, 1))), 'v');
[num_theta_p, den_theta_p] = tfdata(minreal(tf(sys_p(2, 1))), 'v');
TF_u_p = tf(num_u_p, den_u_p);
TF_theta_p = tf(num_theta_p, den_theta_p);

% Convert the state-space system to transfer functions for Thrust input
[num_ut_p, den_ut_p] = tfdata(minreal(tf(sys_p(1, 2))), 'v');
TF_ut_p = tf(num_ut_p, den_ut_p);

% Compare step responses for short period modes
t_small = linspace(0, 10, 1000);
figure;
subplot(2, 1, 1);
step(TF_alpha, TF_alpha_sp, t_small);
title('Short Period Mode - Angle of Attack Approximation');
legend('Δα/Δδ_e', '≈ Δα/Δδ_e');
subplot(2, 1, 2);
step(TF_q, TF_q_sp, t_small);
title('Short Period Mode - Pitch Rate Approximation');
legend('Δq/Δδ_e', '≈ Δq/Δδ_e');

% Compare step responses for phugoid modes
t_long = linspace(0, 1000, 1000);
figure;
subplot(2, 1, 1);
step(TF_u, TF_u_p, t_long);
title('Phugoid Mode - Speed Approximation');
legend('Δu/Δδ_e', '≈ Δu/Δδ_e');
subplot(2, 1, 2);
step(TF_theta, TF_theta_p, t_long);
title('Phugoid Mode - Pitch Angle Approximation');
legend('Δθ/Δδ_e', '≈ Δθ/Δδ_e');

% Design Section

% Plot the Root Locus for pitch angle controller
% Let's use short period approximation: theta-dot= q_sp

%P controller
figure;
subplot(2, 1, 1);
rlocus(num_q_sp * 10, conv([1 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Negative Feedback');

subplot(2, 1, 2);
rlocus(num_q_sp * -10, conv([1 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Positive Feedback/ Negative Gain');

%PI controller: K(-1- 0.8/s)  
figure;
subplot(2, 1, 1);
rlocus(conv(num_q_sp * 10, [-1 -0.8]), conv([1 0 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Negative Feedback');

subplot(2, 1, 2);
rlocus(conv(num_q_sp * 10, [1 0.8]), conv([1 0 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Positive Feedback/ Negative Gain');

%PD controller: K(-1- 0.8s)  
figure;
subplot(2, 1, 1);
rlocus(conv(num_q_sp * 10, [-1 -0.8]), conv([1 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Negative Feedback');

subplot(2, 1, 2);
rlocus(conv(num_q_sp * 10, [1 0.8]), conv([1 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Positive Feedback/ Negative Gain');

%PID controller: K(-1- 0.8s-0.8/s) 
figure;
subplot(2, 1, 1);
rlocus(conv(num_q_sp * 10,[-0.8 -1 -0.8]), conv([1 0 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Negative Feedback');

subplot(2, 1, 2);
rlocus(conv(num_q_sp * 10,[0.8 1 0.8]), conv([1 0 0], conv(den_q_sp, [1, 10]))); % servo model included
title('Root Locus Plot for TF_{\theta} - Positive Feedback/ Negative Gain');



% Plot the Root Locus for speed hold-using thrust

%Openloop transfer function fn
tfu= tf(conv(num_ut_p * 10, [10 1]), conv(conv(den_ut_p, [1, 0.2]), [1, 10]));
% P controller: K
figure;
subplot(2, 1, 1);
rlocus(tfu.Numerator{1}, tfu.Denominator{1} )
title('Root Locus Plot for TF_U - Negative Feedback');
subplot(2, 1, 2);
rlocus(-tfu.Numerator{1}, tfu.Denominator{1} )
title('Root Locus Plot for TF_U - Positive Feedback/ Negative Gain');


% PI controller: K(1+0.8/s)
figure;
subplot(2, 1, 1);
rlocus(conv(tfu.Numerator{1}, [1 0.8]),conv(tfu.Denominator{1}, [1 0]) );
title('Root Locus Plot for TF_U - Negative Feedback');
subplot(2, 1, 2);
rlocus(- conv(tfu.Numerator{1}, [1 0.8]),conv(tfu.Denominator{1}, [1 0]));
title('Root Locus Plot for TF_U - Positive Feedback/ Negative Gain');

% PD controller: k(1+0.8s)
figure;
subplot(2, 1, 1);
tfuu=tf(conv(tfu.Numerator{1}, [0.8 1]), tfu.Denominator{1});
rlocus(tfuu);
title('Root Locus Plot for TF_U - Negative Feedback');
subplot(2, 1, 2);
rlocus(-tfuu);
title('Root Locus Plot for TF_U - Positive Feedback/ Negative Gain');
% Plot the Root Locus for altitude hold
figure;
subplot(2, 1, 1);
rlocus(TF_h);
title('Root Locus Plot for TF_h - Negative Feedback');

subplot(2, 1, 2);
rlocus(-TF_h);
title('Root Locus Plot for TF_h - Positive Feedback/ Negative Gain');

%Inner loop parameters design
AA=[Z_w ,U0 ,-32.2;
    (M_w+ Z_w*M_wd), (M_q+ U0*M_wd) ,0;
    0 ,1 ,0];
BB= [Zde; Mde; 0];

co= ctrb(AA, BB);
pdesired= [-0.7 -5 + 5i -5-5i];
[k]= place(AA, BB, pdesired);
Kw= k(1);
Kq=k(2);
Ktheta= k(3);
