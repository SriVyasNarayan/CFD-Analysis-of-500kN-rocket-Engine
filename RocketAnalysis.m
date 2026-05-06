clc; clear; close all; 
%% ── USER INPUTS ────────────────────────────────────────── 
F = input('Enter required thrust (N): '); 
P_c = input('Enter chamber pressure (Pa): '); 
OF = input('Enter oxygen to fuel ratio: '); 
% Fixed conditions 
P_exit = 101325; % Exit (ambient) pressure [Pa] — sea level 
g0 = 9.80665; % Standard gravity [m/s²] 
R_univ = 8314; % Universal gas constant [J/(kmol·K)] 
cea_table = [ 
1.0, 2972.4, 1.167, 7.24; 
1.5, 3116.9, 1.175, 8.18; 
2.0, 3243.8, 1.182, 9.11; 
2.5, 3397.1, 1.196, 10.12; 
3.0, 3479.6, 1.207, 11.14; 
3.2, 3517.0, 1.213, 11.82; 
3.5, 3545.2, 1.219, 12.20; 
4.0, 3563.6, 1.228, 13.27; 
4.5, 3531.2, 1.237, 14.24; 
5.0, 3474.8, 1.247, 15.17; 
5.5, 3399.5, 1.257, 16.01; 
6.0, 3310.7, 1.265, 16.75; 
6.5, 3211.5, 1.273, 17.40; 
7.0, 3104.3, 1.280, 17.96; 
7.5, 2991.8, 1.286, 18.44; 
8.0, 2875.9, 1.292, 18.86; 
9.0, 2636.8, 1.301, 19.54; 
10.0, 2397.2, 1.309, 20.08; 
]; 
% Clamp O/F to table range 
OF_clamped = max(cea_table(1,1), min(cea_table(end,1), OF)); 
if OF < cea_table(1,1) || OF > cea_table(end,1) 
warning('O/F = %.2f is outside table range [%.1f, %.1f]. Clamped.', ... 
OF, cea_table(1,1), cea_table(end,1)); 
end 
% Interpolate all three properties 
T_c = interp1(cea_table(:,1), cea_table(:,2), OF_clamped, 'pchip'); 
gamma = interp1(cea_table(:,1), cea_table(:,3), OF_clamped, 'pchip'); 
M_mol = interp1(cea_table(:,1), cea_table(:,4), OF_clamped, 'pchip'); 
% Derived gas properties 
R_gas = R_univ / M_mol; % Specific gas constant [J/(kg·K)] 
Cp = gamma * R_gas / (gamma - 1); % Specific heat [J/(kg·K)] 
fprintf('\n--- INPUT PARAMETERS ---\n'); 
fprintf(' Thrust F = %.0f N\n', F); 
fprintf(' Chamber pressure P_c = %.3f MPa\n', P_c/1e6); 
fprintf(' O/F ratio = %.2f\n', OF); 
fprintf(' Exit pressure P_e = %.0f Pa\n', P_exit); 
fprintf('\n--- CEA COMBUSTION GAS PROPERTIES (interpolated) ---\n'); 
fprintf(' Flame temperature T_c = %.2f K\n', T_c); 
fprintf(' Gamma γ = %.4f\n', gamma); 
fprintf(' Mol. weight M_mol = %.3f kg/kmol\n',M_mol); 
fprintf(' Gas constant R = %.3f J/(kg·K)\n',R_gas); 
fprintf(' Cp = %.3f J/(kg·K)\n',Cp); 
%% ── THROAT CONDITIONS (M = 1) ───────────────────────────── 
M_t = 1.0; 
T_t = T_c * (2 / (gamma + 1)); 
P_t = P_c * (2 / (gamma + 1))^(gamma / (gamma - 1)); 
rho_t = P_t / (R_gas * T_t); 
a_t = sqrt(gamma * R_gas * T_t); 
V_t = a_t; % velocity at throat = speed of sound 
fprintf('\n--- THROAT CONDITIONS ---\n'); 
fprintf(' Mach number M_t = %.4f\n', M_t); 
fprintf(' Temperature T_t = %.2f K\n', T_t); 
fprintf(' Pressure P_t = %.4f MPa\n',P_t/1e6); 
fprintf(' Speed of sound a_t = %.2f m/s\n',a_t); 
fprintf(' Density rho_t = %.4f kg/m³\n',rho_t); 
%% ── EXIT MACH NUMBER ───────────────────────────────────── 
% Solve isentropic pressure ratio for M_e (supersonic root) 
pressure_ratio = P_c / P_exit; 
f_mach = @(M) (1 + (gamma-1)/2 * M^2)^(gamma/(gamma-1)) - pressure_ratio; 
options = optimset('Display','off','TolX',1e-12,'TolFun',1e-12); 
M_e = fzero(f_mach, 4.0, options); 
T_e = T_c / (1 + (gamma-1)/2 * M_e^2); 
P_e = P_c / (1 + (gamma-1)/2 * M_e^2)^(gamma/(gamma-1)); 
rho_e = P_e / (R_gas * T_e); 
a_e = sqrt(gamma * R_gas * T_e); 
V_e = M_e * a_e; 
fprintf('\n--- EXIT CONDITIONS ---\n'); 
fprintf(' Mach number M_e = %.4f\n', M_e); 
fprintf(' Temperature T_e = %.2f K\n', T_e); 
fprintf(' Pressure P_e = %.2f Pa\n', P_e); 
fprintf(' Exit velocity V_e = %.2f m/s\n', V_e); 
fprintf(' Density rho_e = %.6f kg/m³\n',rho_e); 
C_star = (1/gamma) * sqrt(gamma * R_gas * T_c) * ... 
((gamma+1)/2)^((gamma+1)/(2*(gamma-1))); 
% Ideal thrust coefficient (perfectly expanded, P_e = P_amb) 
term1 = (2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1)); 
term2 = 1 - (P_exit/P_c)^((gamma-1)/gamma); 
C_F = sqrt(term1 * term2); 
% Note: pressure thrust = 0 when perfectly expanded 
% Effective exhaust velocity 
Ve_eff = C_F * C_star; 
% Specific impulse 
Isp = Ve_eff / g0; 
fprintf('\n--- EXHAUST VELOCITY & PERFORMANCE ---\n'); 
fprintf(' Char. velocity C* = %.2f m/s\n', C_star); 
fprintf(' Thrust coefficient C_F = %.4f\n', C_F); 
fprintf(' Eff. exhaust vel. Ve = %.2f m/s\n', Ve_eff); 
fprintf(' Specific impulse Isp = %.2f s\n', Isp); 
m_dot = F / Ve_eff; 
m_dot_ox = m_dot * OF / (1 + OF); 
m_dot_fuel = m_dot / (1 + OF); 
fprintf('\n--- MASS FLOW RATES ---\n'); 
fprintf(' Total mdot = %.4f kg/s\n', m_dot); 
fprintf(' LOX (O2) mdot_ox = %.4f kg/s\n', m_dot_ox); 
fprintf(' LH2 (H2) mdot_fuel = %.4f kg/s\n', m_dot_fuel); 
% Verify thrust 
F_check = m_dot * Ve_eff; 
fprintf(' Thrust check F = mdot x Ve = %.2f N\n', F_check); 
%% ── THROAT AND EXIT AREAS ──────────────────────────────── 
A_t = m_dot / (rho_t * V_t); 
D_t = 2 * sqrt(A_t / pi); 
% Area-Mach relation at exit 
A_ratio_exit = (1/M_e) * ((2/(gamma+1)) * ... 
(1 + (gamma-1)/2 * M_e^2))^((gamma+1)/(2*(gamma-1))); 
A_e = A_ratio_exit * A_t; 
D_e = 2 * sqrt(A_e / pi); 
fprintf('\n--- NOZZLE AREA RATIOS ---\n'); 
fprintf(' Throat area A_t = %.6f m²\n', A_t); 
fprintf(' Throat diameter D_t = %.4f m (%.1f mm)\n', D_t, D_t*1000); 
fprintf(' Exit/throat ratio eps = %.4f\n', A_ratio_exit); 
fprintf(' Exit area A_e = %.6f m²\n', A_e); 
fprintf(' Exit diameter D_e = %.4f m (%.1f mm)\n', D_e, D_e*1000); 
%% ── CHAMBER GEOMETRY ───────────────────────────────────── 
contraction_ratio = 4.0; 
A_c = contraction_ratio * A_t; 
D_c = 2 * sqrt(A_c / pi); 
L_star = 1.07; % characteristic length H2/O2 [m] 
L_c = L_star * A_t / A_c; 
L_conv = (sqrt(A_c/pi) - sqrt(A_t/pi)) / tan(30*pi/180); 
L_n = 0.8 * (sqrt(A_e/pi) - sqrt(A_t/pi)) / tan(15*pi/180); 
L_total = L_c + L_conv + L_n; 
fprintf('\n--- CHAMBER GEOMETRY ---\n'); 
fprintf(' Contraction ratio Ac/At = %.2f\n', contraction_ratio); 
fprintf(' Chamber diameter D_c = %.4f m (%.1f mm)\n', D_c, D_c*1000); 
fprintf(' Chamber area A_c = %.6f m²\n', A_c); 
fprintf(' Chamber length L_c = %.4f m (%.1f mm)\n', L_c, L_c*1000); 
fprintf(' Convergence length L_cv = %.4f m (%.1f mm)\n', L_conv, L_conv*1000); 
fprintf(' Nozzle length L_n = %.4f m (%.1f mm)\n', L_n, L_n*1000); 
fprintf(' Total engine length L_tot = %.4f m (%.1f mm)\n', L_total, L_total*1000); 
%% ── INJECTOR INLET SIZING (5 inlets, 11 MPa) ───────────── 
P_inj = 11e6; 
T_inj = 300; 
n_inlets = 5; 
M_O2 = 0.032; 
M_H2 = 0.002; 
rho_O2 = P_inj * M_O2 / (R_univ/1000 * T_inj); 
rho_H2 = P_inj * M_H2 / (R_univ/1000 * T_inj); 
v_O2 = 20; 
v_H2 = 30; 
A_O2_per = (m_dot_ox / n_inlets) / (rho_O2 * v_O2); 
A_H2_per = (m_dot_fuel / n_inlets) / (rho_H2 * v_H2); 
D_O2_per = 2 * sqrt(A_O2_per / pi); 
D_H2_per = 2 * sqrt(A_H2_per / pi); 
fprintf('\n--- INJECTOR INLET SIZING (5 inlets each, 11 MPa, 300 K) ---\n'); 
fprintf(' LOX inlet diameter (each) = %.4f m (%.1f mm)\n', D_O2_per, D_O2_per*1000); 
fprintf(' LOX slot height (each) = %.1f mm | velocity = %.1f m/s\n', D_O2_per*1000, 
v_O2); 
fprintf(' LH2 inlet diameter (each) = %.4f m (%.1f mm)\n', D_H2_per, D_H2_per*1000); 
fprintf(' LH2 slot height (each) = %.1f mm | velocity = %.1f m/s\n', D_H2_per*1000, 
v_H2); 
%% ── MACH / T / P DISTRIBUTIONS FOR PLOTTING ────────────── 
M_sub = linspace(0.01, 1.0, 300); 
M_sup = linspace(1.0, M_e, 400); 
AR_fn = @(M) (1./M) .* ((2/(gamma+1)) .* ... 
(1 + (gamma-1)/2 .* M.^2)).^((gamma+1)/(2*(gamma-1))); 
T_fn = @(M) T_c ./ (1 + (gamma-1)/2 .* M.^2); 
P_fn = @(M) P_c ./ (1 + (gamma-1)/2 .* M.^2).^(gamma/(gamma-1)); 
AR_sub = AR_fn(M_sub); T_sub = T_fn(M_sub); P_sub = P_fn(M_sub); 
AR_sup = AR_fn(M_sup); T_sup = T_fn(M_sup); P_sup = P_fn(M_sup); 
%% ── PLOTS ──────────────────────────────────────────────── 
figure('Name','Rocket Engine Design v2', ... 
'Position',[80 80 1150 820], 'Color','white'); 
% 1 — Mach vs Area ratio 
subplot(2,3,1); 
plot(AR_sub, M_sub, 'b-', 'LineWidth',2); hold on; 
plot(AR_sup, M_sup, 'r-', 'LineWidth',2); 
plot(1.0, 1.0, 'ko', 'MarkerFaceColor','k', 'MarkerSize',8); 
plot(A_ratio_exit, M_e, 'rs', 'MarkerFaceColor','r', 'MarkerSize',8); 
xlabel('Area ratio A/A_t [-]'); ylabel('Mach number M [-]'); 
title('Mach number vs Area ratio'); 
legend('Subsonic','Supersonic','Throat M=1', ... 
sprintf('Exit M=%.2f',M_e),'Location','best','FontSize',8); 
grid on; xlim([0 A_ratio_exit+2]); 
% 2 — Temperature vs Mach 
subplot(2,3,2); 
plot(M_sub, T_sub, 'b-', 'LineWidth',2); hold on; 
plot(M_sup, T_sup, 'r-', 'LineWidth',2); 
plot(M_t, T_t, 'ko','MarkerFaceColor','k','MarkerSize',8); 
plot(M_e, T_e, 'rs','MarkerFaceColor','r','MarkerSize',8); 
xlabel('Mach number M [-]'); ylabel('Temperature T [K]'); 
title('Temperature vs Mach number'); 
legend('Subsonic','Supersonic', ... 
sprintf('Throat %.0fK',T_t), ... 
sprintf('Exit %.0fK',T_e),'Location','best','FontSize',8); 
grid on; 
% 3 — Pressure vs Mach 
subplot(2,3,3); 
plot(M_sub, P_sub/1e6, 'b-','LineWidth',2); hold on; 
plot(M_sup, P_sup/1e6, 'r-','LineWidth',2); 
plot(M_t, P_t/1e6,'ko','MarkerFaceColor','k','MarkerSize',8); 
plot(M_e, P_exit/1e6,'rs','MarkerFaceColor','r','MarkerSize',8); 
xlabel('Mach number M [-]'); ylabel('Pressure P [MPa]'); 
title('Pressure vs Mach number'); 
legend('Subsonic','Supersonic', ... 
sprintf('Throat %.3fMPa',P_t/1e6), ... 
sprintf('Exit %.4fMPa',P_exit/1e6),'Location','best','FontSize',8); 
grid on; 
% 4 — Engine contour 
subplot(2,3,4); 
r_sub_flip = fliplr(sqrt(AR_sub * A_t / pi)); 
r_sup_line = sqrt(AR_sup * A_t / pi); 
x_sub = linspace(0, L_c+L_conv, length(M_sub)); 
x_sup = linspace(L_c+L_conv, L_c+L_conv+L_n, length(M_sup)); 
x_all = [x_sub, x_sup(2:end)]; 
r_all = [r_sub_flip, r_sup_line(2:end)]; 
fill([x_all, fliplr(x_all)],[r_all, -fliplr(r_all)], ... 
[0.85 0.92 1.0],'EdgeColor',[0.2 0.4 0.8],'LineWidth',1.5); 
hold on; yline(0,'k--','LineWidth',0.8); 
xline(L_c+L_conv,'k:','Throat','LabelHorizontalAlignment','left','FontSize',8); 
xlabel('Axial position x [m]'); ylabel('Radius r [m]'); 
title('Engine contour (schematic)'); grid on; axis equal; 
% 5 — Isp vs O/F sensitivity 
subplot(2,3,5); 
OF_range = linspace(cea_table(1,1), cea_table(end,1), 120); 
Isp_range = zeros(size(OF_range)); 
Ve_range = zeros(size(OF_range)); 
for k = 1:length(OF_range) 
Tc_k = interp1(cea_table(:,1), cea_table(:,2), OF_range(k), 'pchip'); 
gam_k = interp1(cea_table(:,1), cea_table(:,3), OF_range(k), 'pchip'); 
Mm_k = interp1(cea_table(:,1), cea_table(:,4), OF_range(k), 'pchip'); 
Rg_k = R_univ / Mm_k; 
Cs_k = (1/gam_k)*sqrt(gam_k*Rg_k*Tc_k)*((gam_k+1)/2)^((gam_k+1)/(2*(gam_k-1))); 
t1 = (2*gam_k^2/(gam_k-1))*(2/(gam_k+1))^((gam_k+1)/(gam_k-1)); 
t2 = 1-(P_exit/P_c)^((gam_k-1)/gam_k); 
CF_k = sqrt(t1*t2); 
Ve_range(k) = CF_k * Cs_k; 
Isp_range(k) = Ve_range(k) / g0; 
end 
yyaxis left; 
plot(OF_range, Ve_range,'b-','LineWidth',2); 
ylabel('Ve [m/s]'); 
yyaxis right; 
plot(OF_range, Isp_range,'r--','LineWidth',2); 
ylabel('Isp [s]'); 
xline(OF,'k:',sprintf('O/F=%.1f',OF),'FontSize',9); 
xlabel('O/F [-]'); 
title('Ve and Isp vs O/F (CEA data)'); 
grid on; 
% 6 — Summary 
subplot(2,3,6); axis off; 
txt = { 
'── SUMMARY ─────────────────────────────' 
sprintf(' Thrust F = %.0f N', F) 
sprintf(' Chamber P P_c = %.2f MPa', P_c/1e6) 
sprintf(' O/F ratio = %.2f', OF) 
'' 
'── COMBUSTION ───────────────────────────' 
sprintf(' T_chamber = %.1f K', T_c) 
sprintf(' Gamma γ = %.4f', gamma) 
sprintf(' M_mol = %.3f kg/kmol', M_mol) 
sprintf(' R_gas = %.2f J/kg/K', R_gas) 
'' 
'── PERFORMANCE ──────────────────────────' 
sprintf(' C* = %.2f m/s', C_star) 
sprintf(' C_F = %.4f', C_F) 
sprintf(' Ve (eff.) = %.2f m/s', Ve_eff) 
sprintf(' Isp = %.2f s', Isp) 
sprintf(' Mass flow mdot = %.4f kg/s', m_dot) 
sprintf(' LOX flow mdot_ox = %.4f kg/s', m_dot_ox) 
sprintf(' LH2 flow mdot_f = %.4f kg/s', m_dot_fuel) 
'' 
'── THROAT ───────────────────────────────' 
sprintf(' M_throat = %.4f', M_t) 
sprintf(' T_throat = %.2f K', T_t) 
sprintf(' P_throat = %.4f MPa', P_t/1e6) 
sprintf(' D_throat = %.1f mm', D_t*1000) 
'' 
'── EXIT ─────────────────────────────────' 
sprintf(' M_exit = %.4f', M_e) 
sprintf(' T_exit = %.2f K', T_e) 
sprintf(' D_exit = %.1f mm', D_e*1000) 
sprintf(' Area ratio eps = %.4f', A_ratio_exit) 
'' 
'── CHAMBER ──────────────────────────────' 
sprintf(' D_chamber = %.1f mm', D_c*1000) 
sprintf(' Contraction ratio = %.2f', contraction_ratio) 
sprintf(' Total length = %.1f mm', L_total*1000) 
}; 
text(0.02, 0.99, txt, 'Units','normalized','VerticalAlignment','top', ... 
'FontSize',8,'FontName','Courier New','Color',[0.08 0.08 0.08]); 
sgtitle(sprintf('H_2/O_2 Rocket Engine | F=%.0fN | P_c=%.1fMPa | O/F=%.2f | 
Isp=%.1fs', ... 
F, P_c/1e6, OF, Isp), 'FontSize',12,'FontWeight','bold');
