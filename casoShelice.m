%% ============================================================
%  MECÁNICA DEL VUELO — Práctica: Maniobra en S Vertical
%  Aeronave: CASA C-101EB Aviojet
%  Nivel: Tercero/Cuarto de Ingeniería Aeroespacial
% ============================================================
%  OBJETIVO: Analizar la viabilidad aerodinámica, estructural
%  y propulsiva de una maniobra en S vertical a radio constante.
%
%  FÍSICA DEL PROBLEMA
%  -------------------
%  El avión describe dos semicírculos consecutivos de radio R.
%  El ángulo phi recorre [0, 2pi]:
%
%      Tramo AB  (0 <= phi <= pi):  concavidad HACIA ARRIBA (+n)
%      Tramo BC  (pi < phi <= 2pi): concavidad HACIA ABAJO  (-n)
%
%  Variable de curvatura:
%      curv = +1  en AB   (fuerza centrípeta apunta arriba)
%      curv = -1  en BC   (fuerza centrípeta apunta abajo)
%
%  Equilibrio en el eje NORMAL (centrípeto):
%      curv * L - W*cos(phi) = curv * (W/g) * (V²/R)
%
%  => CL(phi) = W/(qS) * [ curv*(V²/gR) + cos(phi) ]
%  => n(phi)  = curv*(V²/gR) + cos(phi)
%
%  PUNTO CRÍTICO: En B (phi=pi) hay una inversión brusca de n.
%      n_B- = V²/(gR) - 1     (justo antes de B, curv=+1)
%      n_B+ = -V²/(gR) - 1    (justo después de B, curv=-1)
%
%  HIPÓTESIS: V = cte, R = cte, vuelo en plano vertical.
% ============================================================

clear; clc; close all;

%% ============================================================
%  BLOQUE 1 — DATOS DEL AVIÓN
% ============================================================

% -- Pesos y geometría --
W   = 32000;    % [N]      Peso en vuelo
S   = 20.0;     % [m²]     Superficie alar

% -- Atmósfera (ISA, nivel del mar) --
rho = 1.225;    % [kg/m³]  Densidad del aire
g   = 9.81;     % [m/s²]   Gravedad

% -- Polar parabólica: CD = CD0 + k·CL² --
CD0 = 0.030;    % [-]      Resistencia parásita
k   = 0.060;    % [-]      Factor de inducción

% -- Límites aerodinámicos --
CLmax =  1.4;   % [-]      Máximo CL (pérdida positiva)
CLmin = -0.8;   % [-]      Mínimo CL (pérdida negativa)

% -- Límites estructurales (C-101EB, Cat. Acrobática) --
n_max =  7.5;   % [g]      Factor de carga límite positivo
n_min = -3.9;   % [g]      Factor de carga límite negativo

% -- Velocidad máxima (Vne) --
Vne   = 232;    % [m/s]    Velocidad nunca exceder

% -- Empuje disponible (Garrett TFE731-3-1J, ISA SL) --
T_max = 15800;  % [N]      Empuje máximo del motor

%% ============================================================
%  BLOQUE 2 — PARÁMETROS DE LA MANIOBRA
% ============================================================
%  *** MODIFICA AQUÍ V y R para explorar distintos casos ***

V = 140;        % [m/s]    Velocidad de vuelo (constante)
R = 500;        % [m]      Radio de cada semicírculo

% Radios adicionales para la comparativa del Bloque 8
radios_comp = [400, 500, 600, 800];   % [m]

fprintf('=== DATOS DE ENTRADA ===\n');
fprintf('  W = %g N  |  S = %g m²  |  rho = %g kg/m³\n', W, S, rho);
fprintf('  V = %g m/s  |  R = %g m\n', V, R);
fprintf('  CLmax = %g  |  n_max = %g g  |  n_min = %g g\n\n', ...
        CLmax, n_max, n_min);

%% ============================================================
%  BLOQUE 3 — CINEMÁTICA: DISCRETIZACIÓN DE LA TRAYECTORIA
% ============================================================
%
%  phi va de 0 a 2*pi (dos semicírculos completos).
%  Con V constante, phi = omega·t donde omega = V/R.

omega   = V / R;                % [rad/s]  Velocidad angular
t_final = 2 * pi * R / V;      % [s]      Duración total de la S

N   = 1000;
t   = linspace(0, t_final, N);
phi = omega * t;                % [rad]    phi in [0, 2pi]

% -- Variable de curvatura --
% Este es el truco clave: una sola variable permite escribir
% CL y n con una fórmula única válida en ambos tramos.
curv = ones(size(phi));         % +1 por defecto (tramo AB)
curv(phi > pi) = -1;            % -1 en tramo BC

fprintf('=== CINEMÁTICA ===\n');
fprintf('  Velocidad angular:  omega   = %.4f rad/s\n', omega);
fprintf('  Duración total S:   t_final = %.2f s\n', t_final);
fprintf('  Puntos:  A (phi=0)  |  B (phi=pi=%.2f rad)  |  C (phi=2pi)\n\n', pi);

%% ============================================================
%  BLOQUE 4 — DINÁMICA: CL Y FACTOR DE CARGA
% ============================================================
%
%  Ecuación normal (centrípeto, con variable curv):
%
%      CL(phi) = W/(qS) * [ curv*(V²/gR) + cos(phi) ]
%      n(phi)  = curv*(V²/gR) + cos(phi)
%
%  El cambio de signo de curv en phi=pi produce el salto de n en B.

q      = 0.5 * rho * V^2;                               % [Pa]
CL_req = (W/(q*S)) * (curv * (V^2/(g*R)) + cos(phi));   % [-]
n_req  = curv * (V^2/(g*R)) + cos(phi);                 % [g]

% -- Valores en los puntos clave A, B-, B+, C --
n_A    =  V^2/(g*R) + 1;    % phi=0,  curv=+1, cos(0)=+1
n_Bm   =  V^2/(g*R) - 1;    % phi=pi, curv=+1, cos(pi)=-1  (justo antes B)
n_Bp   = -V^2/(g*R) - 1;    % phi=pi, curv=-1, cos(pi)=-1  (justo después B)
n_C    = -V^2/(g*R) + 1;    % phi=2pi,curv=-1, cos(2pi)=+1

fprintf('=== DINÁMICA ===\n');
fprintf('  Presión dinámica q = %.1f Pa\n', q);
fprintf('  n en A   (base):           %+.3f g\n', n_A);
fprintf('  n en B-  (antes del salto):%+.3f g\n', n_Bm);
fprintf('  n en B+  (tras el salto):  %+.3f g  <-- CRITICO\n', n_Bp);
fprintf('  n en C   (final):          %+.3f g\n', n_C);
fprintf('  Salto en B: Delta_n = %.3f g\n\n', n_Bp - n_Bm);

%% ============================================================
%  BLOQUE 5 — DIAGNÓSTICO AERODINÁMICO Y ESTRUCTURAL
% ============================================================

fprintf('=== DIAGNÓSTICO DE VIABILIDAD ===\n');

% 5.1 Límite estructural positivo (ocurre en A)
if n_A > n_max
    fprintf('  [ESTRUCTURA +]: CRITICO  — n en A = %+.2f g > +%.1f g\n', n_A, n_max);
else
    fprintf('  [ESTRUCTURA +]: SEGURO   — n en A = %+.2f g <= +%.1f g\n', n_A, n_max);
end

% 5.2 Límite estructural negativo (ocurre en B+, justo después del salto)
%     ¡Este es el límite dominante en la maniobra en S!
if n_Bp < n_min
    fprintf('  [ESTRUCTURA -]: CRITICO  — n en B+ = %+.2f g < %.1f g\n', n_Bp, n_min);
    fprintf('                  Margen negativo: %.2f g\n', n_Bp - n_min);
else
    fprintf('  [ESTRUCTURA -]: SEGURO   — n en B+ = %+.2f g >= %.1f g\n', n_Bp, n_min);
end

% 5.3 Límite aerodinámico positivo
if any(CL_req > CLmax)
    fprintf('  [AERODINAMICA+]: PERDIDA — CL_req supera CLmax (%.1f) en %d puntos.\n', ...
            CLmax, sum(CL_req > CLmax));
else
    fprintf('  [AERODINAMICA+]: SEGURO  — CL_req max = %.3f <= CLmax = %.1f\n', ...
            max(CL_req), CLmax);
end

% 5.4 Límite aerodinámico negativo
if any(CL_req < CLmin)
    fprintf('  [AERODINAMICA-]: PERDIDA — CL_req cae bajo CLmin (%.1f) en %d puntos.\n', ...
            CLmin, sum(CL_req < CLmin));
else
    fprintf('  [AERODINAMICA-]: SEGURO  — CL_req min = %.3f >= CLmin = %.1f\n', ...
            min(CL_req), CLmin);
end
fprintf('\n');

%% ============================================================
%  BLOQUE 6 — EMPUJE REQUERIDO Y POTENCIA DE EXCESO
% ============================================================
%
%  Eje tangencial (V = cte):
%      T(phi) = D(phi) + W*sin(phi)
%
%  NOTA: sin(phi) cambia de signo de forma natural al pasar
%  de phi=pi a phi=2pi, cubriendo correctamente el descenso.

CD_req = CD0 + k * CL_req.^2;
D_req  = q * S * CD_req;
sin_gamma = abs(sin(phi)); 
T_req = D_req + W * sin_gamma;
P_eje = 1500000;    % Potencia del motor en Watios 
eta_p = 0.8;        % Eficiencia de la hélice (80%)
T_prop = (P_eje * eta_p) / V; % Empuje hélice

% --- Comparación de Déficit ---
% T_req es el vector que ya calculaste (D + W*sin(phi))
deficit_reactor = T_req - 15800;      % 15800 era tu T_max del C-101
deficit_helice  = T_req - T_prop;

Ps     = V * (T_max - T_req) / W;   % [m/s]  Potencia de exceso específica

[Tmax_val, idx_Tmax] = max(T_req);
phi_Tmax = rad2deg(phi(idx_Tmax));

fprintf('=== EMPUJE Y POTENCIA (Reactor vs Hélice) ===\n');
fprintf('  T_req máximo:     %.0f N  en phi = %.1f deg\n', Tmax_val, phi_Tmax);
fprintf('  T_req mínimo:     %.0f N\n', min(T_req));
fprintf('  T_max Reactor:    %.0f N  (TFE731)\n', T_max);
fprintf('  T_max Hélice:     %.0f N  (Motor de %.0f kW)\n', T_prop, P_eje/1000);

% Análisis del Reactor
if Tmax_val > T_max
    fprintf('  [!] DÉFICIT REACTOR: %.0f N\n', Tmax_val - T_max);
else
    fprintf('  Reactor: SUFICIENTE\n');
end

% Análisis de la Hélice
if Tmax_val > T_prop
    fprintf('  [!] DÉFICIT HÉLICE:  %.0f N\n', Tmax_val - T_prop);
else
    fprintf('  Hélice: SUFICIENTE\n');
end

% Diagnóstico final de la hipótesis V=cte
if Tmax_val > T_max || Tmax_val > T_prop
    fraccion = 100 * mean(T_req > T_max);
    fprintf('\n  [!] DIAGNÓSTICO: La hipotesis V=cte no es sostenible.\n');
    fprintf('      El déficit ocurre en el %.0f%% de la maniobra.\n', fraccion);
    fprintf('      En el ascenso el avion desacelerará por falta de empuje.\n\n');
else
    fprintf('\n  DIAGNÓSTICO: V=cte es viable con ambos motores.\n\n');
end
%% ============================================================
%  BLOQUE 7 — ENVOLVENTE V-n
% ============================================================
%
%  La envolvente define la zona volable del avión.
%  Sus fronteras son:
%
%  PARÁBOLAS DE PÉRDIDA (frontera izquierda):
%      n_stall_pos(V) = rho*V²*S*CLmax / (2W)  [rama positiva]
%      n_stall_neg(V) = rho*V²*S*CLmin / (2W)  [rama negativa, <0]
%
%  VELOCIDADES CARACTERÍSTICAS:
%      VS    = sqrt(2W / (rho*S*CLmax))           pérdida 1g positiva
%      VSneg = sqrt(2W / (rho*S*|CLmin|))         pérdida 1g negativa
%      VA    = VS * sqrt(n_max)                   maniobra positiva
%      VAneg: donde n_stall_neg = n_min           maniobra negativa
%
%  Sobre esta envolvente se superponen las curvas de la maniobra S:
%      n_man_A(V) = +V²/(gR) + 1    máximo (en A)
%      n_man_B(V) = -V²/(gR) - 1    mínimo (en B+, tras el salto)

V_vec = linspace(0, Vne*1.15, 500);

% Velocidades características
VS    = sqrt(2*W / (rho*S*CLmax));
VSneg = sqrt(2*W / (rho*S*abs(CLmin)));
VA    = VS * sqrt(n_max);
VAneg = sqrt(2*W*n_min / (rho*S*CLmin));   % CLmin<0, n_min<0 → positivo

% Fronteras de pérdida
n_stall_pos = (rho * V_vec.^2 * S * CLmax) / (2*W);
n_stall_neg = (rho * V_vec.^2 * S * CLmin) / (2*W);

% Curvas de la maniobra S
n_man_A =  1 + V_vec.^2/(g*R);    % máximo (en A)
n_man_B = -1 - V_vec.^2/(g*R);    % mínimo (en B+)

fprintf('=== VELOCIDADES CARACTERÍSTICAS (ENVOLVENTE) ===\n');
fprintf('  VS    = %.1f m/s  (perdida 1g positiva)\n', VS);
fprintf('  VSneg = %.1f m/s  (perdida 1g negativa)\n', VSneg);
fprintf('  VA    = %.1f m/s  (maniobra positiva)\n', VA);
fprintf('  VAneg = %.1f m/s  (maniobra negativa)\n', VAneg);
fprintf('  Vne   = %.1f m/s  (nunca exceder)\n\n', Vne);

% Verificación de viabilidad a V_ensayo
n_man_A_ensayo = 1 + V^2/(g*R);
n_man_B_ensayo = -1 - V^2/(g*R);
fprintf('  A V=%.0f m/s: n_man,A = %+.2f g | n_man,B = %+.2f g\n', ...
        V, n_man_A_ensayo, n_man_B_ensayo);
if n_man_A_ensayo > n_max || n_man_B_ensayo < n_min
    fprintf('  [!] La maniobra SALE de la envolvente a esta velocidad.\n\n');
else
    fprintf('  La maniobra permanece dentro de la envolvente.\n\n');
end

%% ============================================================
%  BLOQUE 8 — COMPARATIVA DE RADIOS
% ============================================================

fprintf('=== COMPARATIVA DE RADIOS (V = %.0f m/s) ===\n', V);
fprintf('  %-10s %-12s %-12s %-12s %-10s\n', ...
        'R [m]','n en A [g]','n en B+ [g]','T_max [N]','Viable?');
fprintf('  %s\n', repmat('-',1,60));

for Ri = radios_comp
    nA_i    =  1 + V^2/(g*Ri);
    nBp_i   = -1 - V^2/(g*Ri);
    CLi_max = (W/(q*S)) * (V^2/(g*Ri) + 1);
    CDi_max = CD0 + k*CLi_max^2;
    Di_max  = q*S*CDi_max;
    Ti_max  = Di_max + W;   % empuje máximo (phi=pi/2, sin=1)

    ok_A  = nA_i  <= n_max;
    ok_B  = nBp_i >= n_min;
    ok_CL = CLi_max <= CLmax;
    viable = ok_A && ok_B && ok_CL;
    str_v = 'SI';
    if ~viable, str_v = 'NO'; end

    fprintf('  %-10d %+11.2f  %+11.2f  %11.0f  %s\n', ...
            Ri, nA_i, nBp_i, Ti_max, str_v);
end
fprintf('\n  Limites: n_max=+%.1f g  |  n_min=%.1f g  |  CLmax=%.1f\n\n', ...
        n_max, n_min, CLmax);

%% ============================================================
%  BLOQUE 9 — GRÁFICAS
% ============================================================

phi_deg = rad2deg(phi);
phi_360 = linspace(0, 360, N);

% ── Figura 1: Trayectoria coloreada por n ────────────────────
%
%  Tramo AB: x = R*sin(phi),       z = R*(1-cos(phi))
%  Tramo BC: la curvatura se invierte. Centro en (0, 2R).
%            x = R*sin(phi),       z = 2R - R*(1+cos(phi-pi))
%                                    = 2R + R*cos(phi)      ... simplificado

idx_AB = phi <= pi;
idx_BC = phi >  pi;

x_tray = zeros(size(phi));
z_tray = zeros(size(phi));

% Tramo AB: centro en (0, 0) → z sube desde -R a +R
x_tray(idx_AB) = R * sin(phi(idx_AB));
z_tray(idx_AB) = -R * cos(phi(idx_AB));          % de -R (base A) a +R (cima B)

% Tramo BC: centro en (0, 2R) → continúa desde B hasta C
x_tray(idx_BC) = R * sin(phi(idx_BC));
z_tray(idx_BC) =  R * cos(phi(idx_BC)) + 2*R;    % de +R (cima B) a +3R (C)

figure('Name','Trayectoria de la Maniobra en S','Color','w', ...
       'Position',[50 50 520 600]);

scatter(x_tray, z_tray, 18, n_req, 'filled');
colormap(jet); cb = colorbar;
cb.Label.String = 'Factor de carga n [g]';

hold on;
% Puntos A, B, C
plot(x_tray(1),       z_tray(1),       'ro','MarkerFaceColor','r','MarkerSize',10);
plot(x_tray(sum(idx_AB)), z_tray(sum(idx_AB)), ...
                                        'bs','MarkerFaceColor','b','MarkerSize',10);
plot(x_tray(end),     z_tray(end),     'g^','MarkerFaceColor','g','MarkerSize',10);
text(x_tray(1)+20,           z_tray(1)-60,         'A','Color','r','FontWeight','bold','FontSize',11);
text(x_tray(sum(idx_AB))+20, z_tray(sum(idx_AB))+40,'B (salto n)','Color','b','FontWeight','bold','FontSize',10);
text(x_tray(end)+20,         z_tray(end)+40,        'C','Color',[0 0.6 0],'FontWeight','bold','FontSize',11);

% Separador entre AB y BC
xline(0,'k--','LineWidth',0.8,'HandleVisibility','off');

annotation('textbox',[0.13 0.01 0.75 0.05], ...
    'String','x(\phi)=R\cdot\sin\phi   |   Tramo AB: z=-R\cos\phi   |   Tramo BC: z=R\cos\phi+2R', ...
    'FitBoxToText','off','EdgeColor','none','FontSize',8,'Interpreter','tex');

axis equal; grid on;
xlabel('x [m]'); ylabel('z [m]');
title('Trayectoria de la maniobra en S (color = factor de carga n)','FontWeight','bold');

% ── Figura 2: CL(phi) ────────────────────────────────────────
figure('Name','CL requerido — Maniobra S','Color','w', ...
       'Position',[590 50 700 350]);

plot(phi_deg, CL_req, 'b-', 'LineWidth', 2, 'DisplayName','C_L requerido');
hold on;
yline(CLmax, 'r--', 'LineWidth', 1.5, 'DisplayName',['C_{L,max} = ' num2str(CLmax)]);
yline(CLmin, 'r:',  'LineWidth', 1.5, 'DisplayName',['C_{L,min} = ' num2str(CLmin)]);
yline(0,     'k-',  'LineWidth', 0.8, 'HandleVisibility','off');

% Zona de pérdida positiva (si la hay)
if any(CL_req > CLmax)
    m = CL_req > CLmax;
    fill([phi_deg(m), fliplr(phi_deg(m))], [CL_req(m), CLmax*ones(1,sum(m))], ...
         'r','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Zona perdida +');
end
% Zona de pérdida negativa (si la hay)
if any(CL_req < CLmin)
    m = CL_req < CLmin;
    fill([phi_deg(m), fliplr(phi_deg(m))], [CL_req(m), CLmin*ones(1,sum(m))], ...
         'm','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Zona perdida -');
end

% Marcadores A, B, C y línea divisoria AB/BC
xline(180,'k--','LineWidth',1.2,'DisplayName','B (phi=180°)');
plot(0,   CL_req(1),   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(360, CL_req(end), 'g^','MarkerFaceColor','g','MarkerSize',8,'HandleVisibility','off');
text(5,   CL_req(1)+0.03,   'A','Color','r','FontWeight','bold');
text(355, CL_req(end)+0.03, 'C','Color',[0 0.6 0],'FontWeight','bold');

% Flecha indicando el salto en B
plot([180 180], [CL_req(sum(idx_AB)), CL_req(sum(idx_AB)+1)], ...
     'k-','LineWidth',2,'HandleVisibility','off');

grid on; legend('Location','best');
xlabel('\phi [°]'); ylabel('C_L [-]');
title('Coeficiente de sustentación requerido en la maniobra S','FontWeight','bold');
xlim([0 360]); xticks(0:45:360);

% ── Figura 3: n(phi) ─────────────────────────────────────────
figure('Name','Factor de carga — Maniobra S','Color','w', ...
       'Position',[590 440 700 350]);

plot(phi_deg, n_req, 'b-', 'LineWidth', 2, 'DisplayName','n requerido');
hold on;
yline(n_max, 'r--', 'LineWidth', 1.5, 'DisplayName',['n_{max} = +' num2str(n_max) ' g']);
yline(n_min, 'r:',  'LineWidth', 1.5, 'DisplayName',['n_{min} = ' num2str(n_min) ' g']);
yline(1,     'k--', 'LineWidth', 0.8, 'DisplayName','n = 1 g');
yline(0,     'k-',  'LineWidth', 0.8, 'HandleVisibility','off');

% Zona fuera de límites
if any(n_req > n_max)
    m = n_req > n_max;
    fill([phi_deg(m),fliplr(phi_deg(m))],[n_req(m),n_max*ones(1,sum(m))], ...
         'r','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Exceso n_max');
end
if any(n_req < n_min)
    m = n_req < n_min;
    fill([phi_deg(m),fliplr(phi_deg(m))],[n_req(m),n_min*ones(1,sum(m))], ...
         'm','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Exceso n_min');
end

% Marcadores y salto en B
xline(180,'k--','LineWidth',1.2,'DisplayName','B (phi=180°)');
plot(0,   n_req(1),   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(360, n_req(end), 'g^','MarkerFaceColor','g','MarkerSize',8,'HandleVisibility','off');
text(5,   n_req(1)+0.2,   'A','Color','r','FontWeight','bold');
text(355, n_req(end)+0.2, 'C','Color',[0 0.6 0],'FontWeight','bold');

% Anotación del salto en B
salto = n_Bp - n_Bm;
annotation('textbox',[0.49 0.55 0.18 0.12], ...
    'String',sprintf('Salto en B:\nΔn = %.2f g\n(n_{B-}=%+.2f → n_{B+}=%+.2f)',salto,n_Bm,n_Bp), ...
    'FitBoxToText','on','BackgroundColor','w','EdgeColor','k','FontSize',8);

grid on; legend('Location','best');
xlabel('\phi [°]'); ylabel('n [g]');
title('Factor de carga en la maniobra S','FontWeight','bold');
xlim([0 360]); xticks(0:45:360);

% ── Figura 4: T_req(phi) ─────────────────────────────────────
figure('Name','Empuje requerido — Maniobra S','Color','w', ...
       'Position',[50 510 700 360]);

plot(phi_deg, T_req/1000, 'k-', 'LineWidth', 2, 'DisplayName','T_{req}');
hold on;
yline(T_max/1000, 'r--', 'LineWidth', 1.5, ...
      'DisplayName',['T_{max} = ' num2str(T_max/1000,'%.1f') ' kN (TFE731)']);
yline(T_prop/1000, '--r', 'Empuje de la Hélice Disponible (8.571 kN)', 'LabelVerticalAlignment', 'top', 'LineWidth', 1.5);

% Zona de déficit
if any(T_req > T_prop)
    m = T_req > T_prop;
    fill([phi_deg(m),fliplr(phi_deg(m))], ...
         [T_req(m)/1000, (T_prop/1000)*ones(1,sum(m))], ...
         'r','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Deficit de empuje');
end

% Zona de Ps positiva (exceso de empuje, en el descenso)
if any(T_req < T_prop)
    m = T_req < T_prop;
    fill([phi_deg(m),fliplr(phi_deg(m))], ...
         [T_req(m)/1000, (T_prop/1000)*ones(1,sum(m))], ...
         'g','FaceAlpha',0.12,'EdgeColor','none','DisplayName','Exceso de empuje (P_s>0)');
end

% Divisor AB/BC y marcadores
xline(180,'k--','LineWidth',1.2,'HandleVisibility','off');
plot(0,   T_req(1)/1000,   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(360, T_req(end)/1000, 'g^','MarkerFaceColor','g','MarkerSize',8,'HandleVisibility','off');
text(5,   T_req(1)/1000+0.5,   'A','Color','r','FontWeight','bold');
text(355, T_req(end)/1000+0.5, 'C','Color',[0 0.6 0],'FontWeight','bold');
text(185, T_max/1000+0.8, 'B','Color','b','FontWeight','bold');

yline(0,'k-','LineWidth',0.8,'HandleVisibility','off');

grid on; legend('Location','best');
xlabel('\phi [°]'); ylabel('T [kN]');
title('Empuje requerido en la maniobra S','FontWeight','bold');
xlim([0 360]); xticks(0:45:360);

% ── Figura 5: Envolvente V-n ─────────────────
figure('Name','Envolvente V-n — Maniobra S','Color','w', ...
       'Position',[790 440 680 430]);

% Polígono de la envolvente (zona volable)
mask_p = V_vec >= VS    & V_vec <= VA;
mask_n = V_vec >= VSneg & V_vec <= VAneg;
Vp = V_vec(mask_p);  np = n_stall_pos(mask_p);
Vn = fliplr(V_vec(mask_n));  nn = fliplr(n_stall_neg(mask_n));

V_fill = [Vp,   VA,   Vne,  Vne,   VAneg, Vn,    VSneg, VS];
n_fill = [np, n_max, n_max, n_min, n_min, nn,      -1,    1];
fill(V_fill, n_fill, [0.85 0.95 0.85], ...
     'EdgeColor','none', 'DisplayName','Zona volable'); hold on;

% Fronteras de pérdida
plot(Vp, np, 'b-', 'LineWidth',2.5, 'DisplayName','Perdida (+)');
plot(Vn, nn, 'b--','LineWidth',2.5, 'DisplayName','Perdida (-)');

% Límites estructurales
plot([VA,    Vne], [n_max, n_max], 'k-','LineWidth',2, ...
     'DisplayName',['n_{max} = +' num2str(n_max) ' g']);
plot([VAneg, Vne], [n_min, n_min], 'k-','LineWidth',2, ...
     'DisplayName',['n_{min} = ' num2str(n_min) ' g']);
plot([Vne,   Vne], [n_min, n_max], 'r-','LineWidth',2, ...
     'DisplayName','V_{ne}');

% Curvas de la maniobra S superpuestas
idx_v = V_vec >= VS & V_vec <= Vne*1.05;
plot(V_vec(idx_v), n_man_A(idx_v), 'm-', 'LineWidth',2.5, ...
     'DisplayName','n_{man,A}(V)  [A, maximo]');
plot(V_vec(idx_v), n_man_B(idx_v), 'm--','LineWidth',2.5, ...
     'DisplayName','n_{man,B}(V)  [B+, minimo]');

% Partes que salen de la envolvente (en rojo)
fuera_A = idx_v & n_man_A > n_max;
fuera_B = idx_v & n_man_B < n_min;
if any(fuera_A)
    plot(V_vec(fuera_A), n_man_A(fuera_A), 'r-','LineWidth',3, ...
         'DisplayName','Fuera de envolvente');
end
if any(fuera_B)
    plot(V_vec(fuera_B), n_man_B(fuera_B), 'r-','LineWidth',3, ...
         'HandleVisibility','off');
end

% Punto de ensayo
plot(V, n_man_A_ensayo, 'ko','MarkerFaceColor','y','MarkerSize',12, ...
     'DisplayName',sprintf('Ensayo A: V=%.0f, n=%+.2f g', V, n_man_A_ensayo));
plot(V, n_man_B_ensayo, 'kd','MarkerFaceColor','c','MarkerSize',12, ...
     'DisplayName',sprintf('Ensayo B+: V=%.0f, n=%+.2f g', V, n_man_B_ensayo));

% Etiquetas velocidades características
text(VS+2,     0.2,    'V_S',      'FontSize',9,'Color','b','FontWeight','bold');
text(VSneg+2, -0.25,   'V_{Sneg}', 'FontSize',9,'Color','b','FontWeight','bold');
text(VA+2,     n_max-0.8,'V_A',    'FontSize',9,'Color','k','FontWeight','bold');
text(VAneg+2,  n_min+0.4,'V_{Aneg}','FontSize',9,'Color','k','FontWeight','bold');
text(Vne+2,    n_max-0.8,'V_{ne}', 'FontSize',9,'Color','r','FontWeight','bold');

yline(0,'k-','LineWidth',0.8,'HandleVisibility','off');
yline(1,'k:','LineWidth',0.8,'HandleVisibility','off');

grid on; legend('Location','northeast','FontSize',8);
xlabel('V [m/s]','FontWeight','bold'); ylabel('n [g]','FontWeight','bold');
title(sprintf('Envolvente V-n + maniobra S  (R = %g m)', R), ...
      'FontWeight','bold');
xlim([0, Vne*1.12]); ylim([n_min-1.5, n_max+1.5]);

% ── Figura 6: Comparativa de radios ──────────────────────────
figure('Name','Comparativa de Radios — Maniobra S','Color','w', ...
       'Position',[790 50 680 360]);

colores = {'b','r','m',[0.8 0.5 0]};
phi_comp = linspace(0, 2*pi, 1000);
curv_comp = ones(size(phi_comp));
curv_comp(phi_comp > pi) = -1;

h_leg = [];
for i = 1:length(radios_comp)
    Ri   = radios_comp(i);
    ni   = curv_comp * (V^2/(g*Ri)) + cos(phi_comp);
    col  = colores{mod(i-1,length(colores))+1};
    h    = plot(rad2deg(phi_comp), ni, '-', 'LineWidth', 2, ...
                'Color', col, 'DisplayName', sprintf('R = %d m', Ri));
    hold on;
    h_leg(end+1) = h; %#ok
end

yline(n_max, 'k--', 'LineWidth', 1.8, 'DisplayName', ...
      ['n_{max} = +' num2str(n_max) ' g']);
yline(n_min, 'k:',  'LineWidth', 1.8, 'DisplayName', ...
      ['n_{min} = '  num2str(n_min) ' g']);
yline(0, 'k-','LineWidth',0.8,'HandleVisibility','off');
xline(180,'k--','LineWidth',1.0,'HandleVisibility','off');

% Marcar qué radios superan n_min en B+
for i = 1:length(radios_comp)
    Ri   = radios_comp(i);
    nBp_i = -1 - V^2/(g*Ri);
    if nBp_i < n_min
        plot(180, nBp_i, 'rx','MarkerSize',14,'LineWidth',3,'HandleVisibility','off');
        text(183, nBp_i-0.2, sprintf('n_{B+}=%.1fg',nBp_i), ...
             'Color','r','FontSize',8);
    end
end

grid on; 
% legend([h_leg, findobj(gca,'Type','ConstantLine')],'Location','best','FontSize',9);
xlabel('\phi [°]'); ylabel('n [g]');
title(sprintf('Sensibilidad del factor de carga al radio  (V = %.0f m/s)', V), ...
      'FontWeight','bold');
xlim([0 360]); xticks(0:45:360);
text(185, n_max*0.92, '\leftarrow B (phi=180°)', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);