%% ============================================================
%  MECÁNICA DEL VUELO — Práctica: Rizo Vertical (Loop)
%  Aeronave: CASA C-101EB Aviojet
%  Nivel: Tercero/Cuarto de Ingeniería Aeroespacial
% ============================================================
%  OBJETIVO: Analizar la viabilidad aerodinámica, estructural
%  y propulsiva de un rizo vertical a radio constante.
%
%  FÍSICA DEL PROBLEMA
%  -------------------
%  El avión describe una circunferencia vertical de radio R.
%  El ángulo phi mide el recorrido desde A (base, phi=0)
%  hasta B (cima, phi=pi). La posición es:
%
%      x(phi) = R·sin(phi)
%      z(phi) = -R·cos(phi)        [origen en el centro]
%
%  El equilibrio en el eje NORMAL (centrípeto) da:
%
%      L - W·cos(phi) = (W/g)·(V²/R)
%
%  de donde se obtiene el CL y el factor de carga requeridos.
%
%  HIPÓTESIS: V = cte, R = cte, vuelo en plano vertical.
% ============================================================

clear; clc; close all;

%% ============================================================
%  BLOQUE 1 — DATOS DEL AVIÓN Y DEL ENSAYO
% ============================================================

% -- Pesos y geometría --
W   = 30000;    % [N]      Peso en vuelo
S   = 20.0;     % [m²]    Superficie alar

% -- Atmósfera (ISA, nivel del mar) --
rho = 1.225;    % [kg/m³]  Densidad del aire
g   = 9.81;     % [m/s²]   Gravedad

% -- Parámetros de la maniobra --
V   = 120;      % [m/s]    Velocidad de vuelo (constante)
R   = 450;      % [m]      Radio del rizo

% -- Polar parabólica: CD = CD0 + k·CL² --
CD0 = 0.025;    % [-]      Resistencia parásita
k   = 0.050;    % [-]      Factor de inducción

% -- Límites aerodinámicos --
CLmax =  1.4;   % [-]      Máximo CL (borde de pérdida positiva)
CLmin = -0.8;   % [-]      Mínimo CL (borde de pérdida negativa)

% -- Límites estructurales (C-101EB, Cat. Acrobática) --
n_max =  7.5;   % [g]      Factor de carga límite positivo
n_min = -3.9;   % [g]      Factor de carga límite negativo

% -- Velocidad máxima (Vne) --
Vne   = 232;    % [m/s]    Velocidad nunca exceder

% -- Empuje disponible (Garrett TFE731-3-1J, ISA SL) --
T_max = 15800;  % [N]      Empuje máximo del motor

fprintf('=== DATOS DE ENTRADA ===\n');
fprintf('  W = %g N  |  S = %g m²  |  rho = %g kg/m³\n', W, S, rho);
fprintf('  V = %g m/s  |  R = %g m\n', V, R);
fprintf('  CLmax = %g  |  n_max = %g g  |  n_min = %g g\n\n', CLmax, n_max, n_min);

%% ============================================================
%  BLOQUE 2 — CINEMÁTICA: DISCRETIZACIÓN DE LA TRAYECTORIA
% ============================================================
%
%  phi va de 0 (base, punto A) a pi (cima, punto B).
%  Con V constante, phi = omega·t donde omega = V/R.

omega = V / R;              % [rad/s]  Velocidad angular
t_man = pi * R / V;         % [s]      Duración del semirizo

N   = 500;                  % Número de puntos de la discretización
t   = linspace(0, t_man, N);
phi = omega * t;            % [rad]   phi ∈ [0, pi]

% Trayectoria en el plano vertical (origen en el centro del rizo)
x = R * sin(phi);           % [m]  Posición horizontal
z = -R * cos(phi);          % [m]  Posición vertical  (z=−R en A, z=+R en B)

fprintf('=== CINEMÁTICA ===\n');
fprintf('  Velocidad angular:  omega = %.4f rad/s\n', omega);
fprintf('  Duración del rizo:  t_man = %.2f s\n\n', t_man);

%% ============================================================
%  BLOQUE 3 — DINÁMICA: CL Y FACTOR DE CARGA
% ============================================================
%
%  Eje NORMAL (perpendicular a la trayectoria, centrípeto):
%
%      L - W·cos(phi) = (W/g)·(V²/R)
%
%      => CL(phi) = (W/qS)·[V²/(gR) + cos(phi)]
%      => n(phi)  = V²/(gR) + cos(phi)
%
%  donde  q = 0.5·rho·V²  es la presión dinámica.

q      = 0.5 * rho * V^2;                      % [Pa]  Presión dinámica

CL_req = (W / (q * S)) * (V^2/(g*R) + cos(phi));  % [-]  CL requerido
n_req  = V^2/(g*R) + cos(phi);                 % [g]   Factor de carga

fprintf('=== DINÁMICA ===\n');
fprintf('  Presión dinámica:   q      = %.1f Pa\n', q);
fprintf('  CL máximo requerido (en A): %.4f\n', max(CL_req));
fprintf('  CL mínimo requerido (en B): %.4f\n', min(CL_req));
fprintf('  n máximo (en A):           %.3f g\n', max(n_req));
fprintf('  n mínimo (en B):           %.3f g\n\n', min(n_req));

%% ============================================================
%  BLOQUE 4 — VELOCIDADES MÍNIMAS DE PÉRDIDA
% ============================================================
%
%  En A (phi=0, base): CL_req es MÁXIMO. Igualando a CLmax:
%
%      Vmin_A = sqrt( 2gR / [rho·g·S·R/W · CLmax − 2] )
%
%  En B (phi=pi, cima): CL_req es MÍNIMO. Igualando a CLmax
%  con signo negativo (vuelo invertido, pérdida por abajo):
%
%      Vmin_B = sqrt( 2gR / [rho·g·S·R/W · CLmax + 2] )
%
%  CONDICIÓN DE EXISTENCIA de Vmin_A:  rho·g·S·R/(2W)·CLmax > 1
%  Si no se cumple, NO existe velocidad que sostenga el rizo.

denom_A = rho*S*CLmax - 2*W/(g*R);
denom_B = rho*S*CLmax + 2*W/(g*R);

fprintf('=== VELOCIDADES MÍNIMAS DE PÉRDIDA ===\n');

if denom_A > 0
    Vmin_A = sqrt(2*W / denom_A);
    fprintf('  Vmin en A (base): %.2f m/s\n', Vmin_A);
else
    Vmin_A = NaN;
    fprintf('  [!] No existe Vmin_A físicamente (denominador <= 0).\n');
    fprintf('      El rizo no puede sostenerse sin pérdida en la base.\n');
end

Vmin_B = sqrt(2*W / denom_B);
fprintf('  Vmin en B (cima): %.2f m/s\n', Vmin_B);
fprintf('  Siempre Vmin_B < Vmin_A: la cima es menos exigente.\n\n');

%% ============================================================
%  BLOQUE 5 — EMPUJE REQUERIDO
% ============================================================
%
%  Eje TANGENCIAL (a lo largo de la trayectoria):
%
%      T - D - W·sin(phi) = 0       [V = cte => aceleración tangencial = 0]
%
%      => T_req(phi) = D(phi) + W·sin(phi)
%
%  donde  D = q·S·CD  y  CD = CD0 + k·CL²

CD_req = CD0 + k * CL_req.^2;      % [-]  CD por polar parabólica
D_req  = q * S * CD_req;            % [N]  Resistencia aerodinámica
T_req  = D_req + W * sin(phi);      % [N]  Empuje requerido

fprintf('=== EMPUJE ===\n');
fprintf('  T_req máximo:     %.0f N  (phi = %.1f°)\n', ...
    max(T_req), rad2deg(phi(T_req == max(T_req))));
fprintf('  T_req mínimo:     %.0f N\n', min(T_req));
fprintf('  T_max disponible: %.0f N  (TFE731)\n', T_max);

if max(T_req) > T_max
    deficit = max(T_req) - T_max;
    fprintf('  [!] DÉFICIT DE EMPUJE: %.0f N. ', deficit);
    fprintf('La hipótesis V=cte no es sostenible.\n\n');
else
    fprintf('  Motor suficiente. V = cte es viable.\n\n');
end

%% ============================================================
%  BLOQUE 6 — DIAGNÓSTICO DE VIABILIDAD
% ============================================================

fprintf('=== DIAGNÓSTICO DE VIABILIDAD ===\n');

% 6.1 Comprobación estructural
if max(n_req) > n_max
    fprintf('  [ESTRUCTURA +]: CRÍTICO  — n_max = %.2f g > %.1f g\n', max(n_req), n_max);
else
    fprintf('  [ESTRUCTURA +]: SEGURO   — n_max = %.2f g <= %.1f g\n', max(n_req), n_max);
end

if min(n_req) < n_min
    fprintf('  [ESTRUCTURA -]: CRÍTICO  — n_min = %.2f g < %.1f g\n', min(n_req), n_min);
else
    fprintf('  [ESTRUCTURA -]: SEGURO   — n_min = %.2f g >= %.1f g\n', min(n_req), n_min);
end

% 6.2 Comprobación aerodinámica
if any(CL_req > CLmax)
    fprintf('  [AERODINÁMICA]: PÉRDIDA  — CL_req supera CLmax en %d puntos.\n', sum(CL_req > CLmax));
elseif any(CL_req < CLmin)
    fprintf('  [AERODINÁMICA]: PÉRDIDA  — CL_req cae bajo CLmin en %d puntos.\n', sum(CL_req < CLmin));
else
    fprintf('  [AERODINÁMICA]: SEGURO   — CL_req ∈ [%.2f, %.2f] en toda la maniobra.\n', min(CL_req), max(CL_req));
end

% 6.3 Comprobación de velocidad
if ~isnan(Vmin_A) && V >= Vmin_A
    fprintf('  [VELOCIDAD]:    SEGURO   — V = %g m/s >= Vmin_A = %.2f m/s\n', V, Vmin_A);
elseif ~isnan(Vmin_A)
    fprintf('  [VELOCIDAD]:    CRÍTICO  — V = %g m/s < Vmin_A = %.2f m/s\n', V, Vmin_A);
end

% 6.4 Comprobación propulsiva
if max(T_req) <= T_max
    fprintf('  [PROPULSIÓN]:   SEGURO   — T_req_max = %.0f N <= T_max = %.0f N\n', max(T_req), T_max);
else
    fprintf('  [PROPULSIÓN]:   DÉFICIT  — T_req_max = %.0f N > T_max = %.0f N\n', max(T_req), T_max);
end
fprintf('\n');

%% ============================================================
%  BLOQUE 7 — ENVOLVENTE V-n
% ============================================================
%
%  
%  Es la frontera estructural/aerodinámica del avión, definida por:
%
%  FRONTERA IZQUIERDA — Pérdida (stall):
%      Rama positiva: n_stall_pos(V) = rho·V²·S·CLmax / (2W)
%      Rama negativa: n_stall_neg(V) = rho·V²·S·CLmin / (2W)   [CLmin < 0]
%
%  VELOCIDADES CARACTERÍSTICAS:
%      VS    = sqrt(2W / (rho·S·CLmax))          1-g stall speed
%      VSneg = sqrt(2W / (rho·S·|CLmin|))        1-g negative stall speed
%      VA    = VS · sqrt(n_max)                   design maneuvering speed
%      VAneg: donde n_stall_neg(V) = n_min        maniobrabilidad negativa
%
%  ZONA VOLABLE (interior de la envolvente):
%      n_stall_neg(V) <= n <= n_stall_pos(V)     [zona de pérdida]
%      n_min          <= n <= n_max               [zona estructural]
%      V <= Vne                                   [velocidad máxima]
%
%  La curva del rizo n_man(V) = 1 + V²/(gR) se superpone
%  para comprobar si la maniobra CABE dentro de la envolvente.
%  ¡Puede haber puntos fuera de n_man sin tocar n_max!

V_vec = linspace(0, Vne*1.15, 600);

% -- Velocidades características --
VS    = sqrt(2*W / (rho*S*CLmax));          % [m/s] Velocidad de pérdida 1g
VSneg = sqrt(2*W / (rho*S*abs(CLmin)));     % [m/s] Velocidad de pérdida 1g negativa
VA    = VS * sqrt(n_max);                   % [m/s] Velocidad de maniobra

% VAneg: velocidad donde la rama negativa toca n_min
%   n_stall_neg = rho*V²*S*CLmin/(2W) = n_min  => V² = 2W*n_min/(rho*S*CLmin)
%   Como CLmin<0 y n_min<0, el cociente es positivo.
VAneg = sqrt(2*W*n_min / (rho*S*CLmin));    % [m/s] Maniobra negativa

fprintf('=== VELOCIDADES CARACTERÍSTICAS (ENVOLVENTE) ===\n');
fprintf('  VS    (pérdida 1g positiva): %.1f m/s\n', VS);
fprintf('  VSneg (pérdida 1g negativa): %.1f m/s\n', VSneg);
fprintf('  VA    (maniobra positiva):   %.1f m/s\n', VA);
fprintf('  VAneg (maniobra negativa):   %.1f m/s\n', VAneg);
fprintf('  Vne   (nunca exceder):       %.1f m/s\n\n', Vne);

% -- Fronteras de pérdida (parábolas) --
n_stall_pos = (rho * V_vec.^2 * S * CLmax)  / (2*W);   % rama positiva
n_stall_neg = (rho * V_vec.^2 * S * CLmin)  / (2*W);   % rama negativa (valores <0)

% -- Curva del rizo superpuesta --
n_man = 1 + V_vec.^2 / (g*R);      % Factor de carga máximo del rizo (en A)

% -- Construcción del perímetro de la envolvente (para el fill) --
%  Sentido antihorario: rama pos (VS->VA) + techo (VA->Vne) +
%  pared der. (n_max->n_min) + suelo (Vne->VAneg) + rama neg (VAneg->VSneg)

% Rama positiva de pérdida: VS a VA
mask_pos = V_vec >= VS & V_vec <= VA;
Vp = V_vec(mask_pos);  np = n_stall_pos(mask_pos);

% Techo estructural: VA a Vne
Vtecho = [VA, Vne];  ntecho = [n_max, n_max];

% Pared derecha: n_max a n_min (en Vne)
Vder = [Vne, Vne];  nder = [n_max, n_min];

% Suelo estructural: Vne a VAneg
Vsuelo = [Vne, VAneg];  nsuelo = [n_min, n_min];

% Rama negativa de pérdida: VAneg a VSneg
mask_neg = V_vec >= VSneg & V_vec <= VAneg;
Vn = fliplr(V_vec(mask_neg));   nn = fliplr(n_stall_neg(mask_neg));

% Cierre (VSneg de vuelta a VS con n=0 aproximadamente)
Vcierre = [VSneg, VS];  ncierre = [-1, 1];

% Coordenadas del polígono completo
V_fill = [Vp,    Vtecho, Vder, Vsuelo, Vn,    Vcierre];
n_fill = [np,    ntecho, nder, nsuelo, nn,    ncierre];

%% ============================================================
%  BLOQUE 8 — GRÁFICAS
% ============================================================

phi_deg = rad2deg(phi);     % Convertir a grados para los ejes

% ── Figura 1: Trayectoria coloreada por factor de carga ──────
figure('Name','Trayectoria del Rizo','Color','w','Position',[50 50 500 520]);

scatter(x, z, 18, n_req, 'filled');
colormap(jet); cb = colorbar;
cb.Label.String = 'Factor de carga n [g]';
clim([min(n_req)-0.1, max(n_req)+0.1]);

hold on;
% Marcar A y B
plot(x(1),   z(1),   'ro', 'MarkerFaceColor','r', 'MarkerSize',10);
plot(x(end), z(end), 'bs', 'MarkerFaceColor','b', 'MarkerSize',10);
text(x(1)+15,   z(1)-50,  'A  (n_{max})', 'Color','r', 'FontWeight','bold');
text(x(end)+15, z(end)+30,'B  (n_{min})', 'Color','b', 'FontWeight','bold');

% Ecuaciones paramétricas
annotation('textbox',[0.13 0.02 0.5 0.07], ...
    'String','x(\phi)=R\cdot\sin\phi \quad z(\phi)=-R\cdot\cos\phi', ...
    'FitBoxToText','on','EdgeColor','none','FontSize',9, ...
    'Interpreter','tex');

axis equal; grid on;
xlabel('x [m]'); ylabel('z [m]');
title('Trayectoria del rizo vertical (color = factor de carga)','FontWeight','bold');

% ── Figura 2: CL(phi) ────────────────────────────────────────
figure('Name','CL requerido','Color','w','Position',[570 50 650 340]);

plot(phi_deg, CL_req, 'b-', 'LineWidth', 2, 'DisplayName','C_L requerido');
hold on;
yline(CLmax, 'r--', 'LineWidth', 1.5, 'DisplayName', ...
    ['C_{L,max} = ' num2str(CLmax)]);
yline(CLmin, 'r:', 'LineWidth', 1.5, 'DisplayName', ...
    ['C_{L,min} = ' num2str(CLmin)]);
yline(0, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');

% Zona de pérdida (si la hay)
if any(CL_req > CLmax)
    area_x = phi_deg(CL_req > CLmax);
    area_y = CL_req(CL_req > CLmax);
    fill([area_x, fliplr(area_x)], [area_y, CLmax*ones(size(area_y))], ...
        'r', 'FaceAlpha', 0.25, 'EdgeColor','none','DisplayName','Zona de pérdida');
end

% Puntos A y B
plot(0,   CL_req(1),   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(180, CL_req(end), 'bs','MarkerFaceColor','b','MarkerSize',8,'HandleVisibility','off');
text(5,   CL_req(1)+0.03,  'A','Color','r','FontWeight','bold');
text(173, CL_req(end)+0.03,'B','Color','b','FontWeight','bold');

grid on;
xlabel('\phi [°]'); ylabel('C_L [-]');
title('Coeficiente de sustentación requerido a lo largo del rizo','FontWeight','bold');
legend('Location','northeast');
xlim([0 180]); xticks(0:30:180);

% ── Figura 3: n(phi) ─────────────────────────────────────────
figure('Name','Factor de carga','Color','w','Position',[570 430 650 340]);

plot(phi_deg, n_req, 'b-', 'LineWidth', 2, 'DisplayName','n requerido');
hold on;
yline(n_max, 'r--', 'LineWidth', 1.5, 'DisplayName', ...
    ['n_{max} = ' num2str(n_max) ' g']);
yline(n_min, 'r:',  'LineWidth', 1.5, 'DisplayName', ...
    ['n_{min} = ' num2str(n_min) ' g']);
yline(1, 'k--', 'LineWidth', 0.8, 'DisplayName','n = 1 g (nivel)');
yline(0, 'k-',  'LineWidth', 0.8, 'HandleVisibility','off');

plot(0,   n_req(1),   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(180, n_req(end), 'bs','MarkerFaceColor','b','MarkerSize',8,'HandleVisibility','off');
text(5,   n_req(1)+0.15,  'A','Color','r','FontWeight','bold');
text(173, n_req(end)+0.15,'B','Color','b','FontWeight','bold');

grid on;
xlabel('\phi [°]'); ylabel('n [g]');
title('Factor de carga a lo largo del rizo','FontWeight','bold');
legend('Location','northeast');
xlim([0 180]); xticks(0:30:180);

% ── Figura 4: T_req(phi) ─────────────────────────────────────
figure('Name','Empuje requerido','Color','w','Position',[50 580 650 340]);

plot(phi_deg, T_req/1000, 'k-',  'LineWidth', 2, 'DisplayName','T_{req}');
hold on;
yline(T_max/1000, 'r--', 'LineWidth', 1.5, 'DisplayName', ...
    ['T_{max} = ' num2str(T_max/1000,'%.1f') ' kN (TFE731)']);

% Zona de déficit (si la hay)
if any(T_req > T_max)
    area_x = phi_deg(T_req > T_max);
    area_y = T_req(T_req > T_max)/1000;
    fill([area_x, fliplr(area_x)], [area_y, (T_max/1000)*ones(size(area_y))], ...
        'r','FaceAlpha',0.25,'EdgeColor','none','DisplayName','Déficit de empuje');
end

plot(0,   T_req(1)/1000,   'ro','MarkerFaceColor','r','MarkerSize',8,'HandleVisibility','off');
plot(180, T_req(end)/1000, 'bs','MarkerFaceColor','b','MarkerSize',8,'HandleVisibility','off');
text(5,   T_req(1)/1000+0.2,  'A','Color','r','FontWeight','bold');
text(173, T_req(end)/1000+0.2,'B','Color','b','FontWeight','bold');

grid on;
xlabel('\phi [°]'); ylabel('T [kN]');
title('Empuje requerido a lo largo del rizo','FontWeight','bold');
legend('Location','northeast');
xlim([0 180]); xticks(0:30:180);

% ── Figura 5: Envolvente V-n ─────
%
%  La zona VERDE es el interior de la envolvente.
%  La zona de pérdida queda FUERA (a la izquierda de las parábolas).
%  La curva del rizo n_man(V) se superpone: si sale de la zona verde,
%  la maniobra no es viable a esa velocidad.

figure('Name','Envolvente V-n','Color','w','Position',[720 580 700 420]);

% -- Zona volable (fill del polígono de la envolvente) --
fill(V_fill, n_fill, [0.85 0.95 0.85], ...
    'EdgeColor','none', 'DisplayName','Zona volable (envolvente)');
hold on;

% -- Fronteras de pérdida (parábolas) --
%    Solo se dibujan hasta VA (positiva) y VAneg (negativa)
V_pos = V_vec(V_vec >= VS  & V_vec <= VA);
V_neg = V_vec(V_vec >= VSneg & V_vec <= VAneg);

plot(V_pos, n_stall_pos(V_vec >= VS  & V_vec <= VA), ...
    'b-', 'LineWidth', 2.5, 'DisplayName','Frontera pérdida (+)');
plot(V_neg, n_stall_neg(V_vec >= VSneg & V_vec <= VAneg), ...
    'b--','LineWidth', 2.5, 'DisplayName','Frontera pérdida (-)');

% -- Límites estructurales horizontales --
plot([VA,  Vne], [n_max, n_max], 'k-', 'LineWidth', 2.5, ...
    'DisplayName',['n_{max} = +' num2str(n_max) ' g']);
plot([VAneg, Vne], [n_min, n_min], 'k-', 'LineWidth', 2.5, ...
    'DisplayName',['n_{min} = ' num2str(n_min) ' g']);

% -- Pared derecha (Vne) --
plot([Vne, Vne], [n_min, n_max], 'r-', 'LineWidth', 2.5, ...
    'DisplayName','V_{ne}');

% -- Zona de pérdida (a la izquierda, sombreada en rojo claro) --
%    Positiva: entre V=0 y la parábola positiva (hasta VS)
V_shade = linspace(0, VS, 50);
fill([V_shade, fliplr(V_shade)], ...
     [n_stall_pos(1:50), zeros(1,50)], ...
     [1 0.88 0.88], 'EdgeColor','none', 'HandleVisibility','off');
%    Negativa: zona bajo la parábola negativa
fill([V_shade, fliplr(V_shade)], ...
     [n_stall_neg(1:50), zeros(1,50)], ...
     [1 0.88 0.88], 'EdgeColor','none', 'HandleVisibility','off');

% -- Curva del rizo superpuesta (la que usamos en la práctica) --
%    Esta curva debe quedar DENTRO de la zona verde para que la maniobra sea viable
idx_visible = V_vec >= VS & V_vec <= Vne*1.1;
plot(V_vec(idx_visible), n_man(idx_visible), ...
    'm-', 'LineWidth', 2.5, 'DisplayName', ...
    sprintf('n_{man}(V) rizo  R=%gm', R));

% Trazar en rojo la parte que sale de la envolvente (n_man > n_max)
idx_fuera = idx_visible & n_man > n_max;
if any(idx_fuera)
    plot(V_vec(idx_fuera), n_man(idx_fuera), ...
        'r-', 'LineWidth', 3, 'DisplayName','Rizo fuera de envolvente');
end

% -- Punto de ensayo --
plot(V, max(n_req), 'ko', 'MarkerFaceColor','y', 'MarkerSize',12, ...
    'DisplayName',['Ensayo: V=' num2str(V) ' m/s  (n=' ...
    num2str(max(n_req),'%.2f') 'g)']);

% -- Etiquetas de velocidades características --
text(VS+1,    0.15, 'V_S',    'FontSize',9, 'Color','b', 'FontWeight','bold');
text(VSneg+1, -0.2, 'V_{Sneg}','FontSize',9,'Color','b', 'FontWeight','bold');
text(VA+1,    n_max-0.6, 'V_A',  'FontSize',9,'Color','k','FontWeight','bold');
text(VAneg+1, n_min+0.3, 'V_{Aneg}','FontSize',9,'Color','k','FontWeight','bold');
text(Vne+2,   n_max-0.6, 'V_{ne}','FontSize',9,'Color','r','FontWeight','bold');

% -- Línea n=0 y n=1g de referencia --
yline(0, 'k-',  'LineWidth', 0.8, 'HandleVisibility','off');
yline(1, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
text(Vne*1.05, 1.1, 'n=1g', 'FontSize', 8, 'Color', [0.4 0.4 0.4]);

% -- Anotación explicativa --
annotation('textbox', [0.13 0.01 0.75 0.06], ...
    'String', ...
    ['Zona verde = envolvente (zona volable). ' ...
     'La curva magenta es n_{man}(V) = 1 + V^2/(gR) del rizo. ' ...
     'Si sale de la zona verde, la maniobra es INVIABLE a esa V.'], ...
    'FitBoxToText','off', 'EdgeColor','none', ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3]);

grid on;
xlabel('V [m/s]', 'FontWeight','bold');
ylabel('n [g]',   'FontWeight','bold');
title(sprintf('Envolvente V-n con curva del rizo (R = %g m)', R), ...
    'FontWeight','bold');
legend('Location','northeast', 'FontSize', 8);
xlim([0, Vne*1.12]);
ylim([n_min - 1.5, n_max + 1.5]);