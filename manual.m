% Verify manually for known optimal values
kp_opt = 0;
kr_opt = 0.95;
Gpr_manual = tf([kp_opt kr_opt kp_opt*((2*pi*50)^2)], [1 0 (2*pi*50)^2]);
sys_manual = series(Gpr_manual, G);
figure;
bode(sys_manual);
grid on;
margin(sys_manual);
grid on;
