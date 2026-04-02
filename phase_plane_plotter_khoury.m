%% =========================================================
%  PHASE PLANE & VECTOR FIELD PLOTTER
%  for a user-defined 2D system of ODEs
%  =========================================================
%
%  REQUIRES: Symbolic Math Toolbox
%
%  STEPS:
%    1. Declare symbolic variables.
%    2. Define your equations f1, f2.
%    3. Set parameters, plot window, and options.
%    4. Run — the phase plane, vector field, nullclines,
%       equilibria, and sample trajectories are all plotted.
% =========================================================

clear; clc; close all;

%% ═══════════════════════════════════════════════════════
%  1 — SYMBOLIC VARIABLES
% ═══════════════════════════════════════════════════════
syms x y real

% ---- EXTRA PARAMETERS -----------------
syms mu alpha sigma L w real
% --------------------------------------------------------

%% ═══════════════════════════════════════════════════════
%  2 — USER-DEFINED EQUATIONS
%  Define the right-hand side of:
%    dx/dt = f1(x, y)
%    dy/dt = f2(x, y)
%  Here, I use the model by Khoury et al, 2011, used in the article
% ═══════════════════════════════════════════════════════
f1 =  L*(x+y)/(w+x+y) - x*(alpha - sigma*y/(x+y));
f2 =  x*(alpha - sigma*y/(x+y)) - mu*y;

f_sym = [f1; f2];

%% ═══════════════════════════════════════════════════════
%  3 — PARAMETERS & PLOT SETTINGS
% ═══════════════════════════════════════════════════════
a = 0.25; s = 0.75; LL = 2000; ww = 27000;
mu_values = [0.2 0.4 (LL/ww)*(1/(2*(a-LL/ww)))*(a+s+sqrt(a^2+4*s*LL/ww-2*s*a+s^2))];

% ---- Numeric values for any extra parameters -----------
param_syms = [mu alpha sigma L w ];          % symbolic names 
param_vals = [mu_values(3) a s LL ww];         % numeric values  (same order)
% Leave both as [] if there are no extra parameters.

% ---- Plot window ----------------------------------------
x_range = [0,  15000];    % [x_min, x_max]
y_range = [0,  5000];    % [y_min, y_max]

% ---- Vector field grid resolution -----------------------
n_grid = 22;                % grid points per axis (≈20 is clean)

% ---- Trajectory initial conditions ----------------------
% Each row is one initial condition [x0, y0]
IC = [
     9000.0,  2000.0;
     9000.0,  0.0;
     10,      10.0;
     5000,    5000;
     2000,    4000;
     8000,    100;
];

t_span = [0, 500];           % integration time for each trajectory

% ---- Display toggles ------------------------------------
show_nullclines  = true;    % x-nullcline and y-nullcline
show_equilibria  = true;    % locate and mark equilibria
show_trajectories = true;   % integrate and draw trajectories

%% ═══════════════════════════════════════════════════════
%  4 — BUILD NUMERIC FUNCTIONS  
% ═══════════════════════════════════════════════════════

% Substitute extra parameters into symbolic expressions
if isempty(param_syms)
    f1_sub = f1;
    f2_sub = f2;
else
    f1_sub = subs(f1, param_syms, param_vals);
    f2_sub = subs(f2, param_syms, param_vals);
end

% Convert to fast numeric handles
F1 = matlabFunction(f1_sub, 'Vars', [x, y]);
F2 = matlabFunction(f2_sub, 'Vars', [x, y]);
rhs = @(t, X) [F1(X(1), X(2)); F2(X(1), X(2))];

% Print system info
fprintf('==============================================\n');
fprintf('  PHASE PLANE PLOTTER\n');
fprintf('==============================================\n');
fprintf('System:\n');
fprintf('  dx/dt = %s\n', char(f1_sub));
fprintf('  dy/dt = %s\n\n', char(f2_sub));

%% ═══════════════════════════════════════════════════════
%  5 — FIND EQUILIBRIA  (symbolic solve)
% ═══════════════════════════════════════════════════════
eq_points = [];

if show_equilibria
    fprintf('Searching for equilibria...\n');
    try
        sol_sym = solve([f1_sub == 0, f2_sub == 0], [x, y], ...
                        'Real', true, 'ReturnConditions', false);
        if isstruct(sol_sym)
            xs = double(sol_sym.x);
            ys = double(sol_sym.y);
        else
            xs = double(sol_sym{1});
            ys = double(sol_sym{2});
        end
        % Keep only real, finite solutions inside the window
        mask = isreal(xs) & isreal(ys) & isfinite(xs) & isfinite(ys) & ...
               xs >= x_range(1) & xs <= x_range(2) & ...
               ys >= y_range(1) & ys <= y_range(2);
        eq_points = [xs(mask), ys(mask)];
    catch
        fprintf('  Symbolic solve failed — trying numerical search.\n');
    end

    % Numerical fallback: multi-start fsolve on the window
    opts_fs = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);
    F_root  = @(X) [F1(X(1),X(2)); F2(X(1),X(2))];
    rng(0);
    n_try = 80;
    x_starts = x_range(1) + rand(n_try,1)*diff(x_range);
    y_starts = y_range(1) + rand(n_try,1)*diff(y_range);
    for k = 1:n_try
        try
            p = fsolve(F_root, [x_starts(k); y_starts(k)], opts_fs);
            if norm(F_root(p)) < 1e-8 && ...
               p(1) >= x_range(1) && p(1) <= x_range(2) && ...
               p(2) >= y_range(1) && p(2) <= y_range(2)
                % Deduplicate
                if isempty(eq_points) || min(vecnorm(eq_points - p', 2, 2)) > 1e-4
                    eq_points = [eq_points; p']; %#ok<AGROW>
                end
            end
        catch
        end
    end

    if isempty(eq_points)
        fprintf('  No equilibria found in the plot window.\n');
    else
        fprintf('  Found %d equilibrium/a:\n', size(eq_points,1));
        J_sym = jacobian(f_sym, [x, y]);
        for k = 1:size(eq_points,1)
            px = eq_points(k,1);  py = eq_points(k,2);
            J_num = double(subs(J_sym, [x,y,param_syms], [px,py,param_vals]));
            ev    = eig(J_num);
            re    = real(ev);
            im    = imag(ev);
            if max(re) < -1e-8
                if max(abs(im)) > 1e-8, stab = 'Stable spiral';
                else,                    stab = 'Stable node'; end
            elseif max(re) > 1e-8
                if min(re) < -1e-8,     stab = 'Saddle';
                elseif max(abs(im))>1e-8,stab = 'Unstable spiral';
                else,                    stab = 'Unstable node'; 
                end
            else
                stab = 'Non-hyperbolic / center';
            end
            fprintf('    E%d: (%.4f, %.4f)  →  %s\n', k, px, py, stab);
        end
    end
    fprintf('\n');
end

%% ═══════════════════════════════════════════════════════
%  6 — GRID FOR VECTOR FIELD
% ═══════════════════════════════════════════════════════
xv = linspace(x_range(1), x_range(2), n_grid);
yv = linspace(y_range(1), y_range(2), n_grid);
[Xg, Yg] = meshgrid(xv, yv);

U = zeros(size(Xg));
V = zeros(size(Xg));
for r = 1:size(Xg,1)
    for c = 1:size(Xg,2)
        U(r,c) = F1(Xg(r,c), Yg(r,c));
        V(r,c) = F2(Xg(r,c), Yg(r,c));
    end
end

% Normalize arrows to uniform length for readability
mag  = sqrt(U.^2 + V.^2) + 1e-12;
Un   = U ./ mag;
Vn   = V ./ mag;
% Color arrows by log-speed for perceptual clarity
logmag = log1p(mag);

%% ═══════════════════════════════════════════════════════
%  7 — NULLCLINES  (symbolic)
% ═══════════════════════════════════════════════════════
xv_fine = linspace(x_range(1), x_range(2), 600);
yv_fine = linspace(y_range(1), y_range(2), 600);
[Xnc, Ync] = meshgrid(xv_fine, yv_fine);

NC1 = zeros(size(Xnc));   % f1 = 0  (x-nullcline)
NC2 = zeros(size(Xnc));   % f2 = 0  (y-nullcline)
for r = 1:size(Xnc,1)
    for c = 1:size(Xnc,2)
        NC1(r,c) = F1(Xnc(r,c), Ync(r,c));
        NC2(r,c) = F2(Xnc(r,c), Ync(r,c));
    end
end

%% ═══════════════════════════════════════════════════════
%  8 — INTEGRATE TRAJECTORIES
% ═══════════════════════════════════════════════════════
trajs = cell(size(IC,1), 1);
if show_trajectories
    ode_opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, ...
                      'Events', @(t,X) exit_event(t, X, x_range, y_range));
    for k = 1:size(IC,1)
        try
            [~, S] = ode45(rhs, t_span, IC(k,:)', ode_opts);
            trajs{k} = S;
        catch
            trajs{k} = [];
        end
    end
end

%% ═══════════════════════════════════════════════════════
%  9 — FIGURE
% ═══════════════════════════════════════════════════════
figure('Color','w','Position',[60 60 880 680], ...
       'Name','Phase Plane & Vector Field');

ax = axes('Position',[0.10 0.10 0.82 0.83]);
hold(ax,'on');

% ── 9a. Vector field (colored by speed) ─────────────────
scale = 0.85 * min(diff(x_range), diff(y_range)) / n_grid;
q = quiver(ax, Xg, Yg, Un*scale, Vn*scale, 0, ...
           'LineWidth', 1.1, 'MaxHeadSize', 0.6);

% Color each arrow by log-speed
q.XData = Xg;  q.YData = Yg;
% Build per-arrow RGB from colormap
cmap    = colormap(ax, parula(256));
logmag_flat = logmag(:);
lmin = min(logmag_flat);  lmax = max(logmag_flat);
cidx = round(1 + 255*(logmag_flat - lmin)/(lmax - lmin + eps));
cidx(isnan(cidx)) = 1;
arrow_colors = cmap(cidx, :);

% Manually draw colored arrows (quiver doesn't support per-arrow color)
delete(q);
arrow_scale = scale;
for r = 1:size(Xg,1)
    for c = 1:size(Xg,2)
        col = arrow_colors((r-1)*size(Xg,2)+c, :);
        x0 = Xg(r,c); y0 = Yg(r,c);
        dx = Un(r,c)*arrow_scale; dy = Vn(r,c)*arrow_scale;
        annotation_arrow = plot(ax, [x0, x0+dx], [y0, y0+dy], ...
            '-', 'Color', [col, 0.75], 'LineWidth', 0.9);
        % Arrowhead dot
        plot(ax, x0+dx, y0+dy, '.', 'Color', col, 'MarkerSize', 5);
    end
end

% ── 9b. Colorbar for speed ───────────────────────────────
cb = colorbar(ax);
cb.Label.String = 'log(1 + |f|)  (flow speed)';
cb.Label.FontSize = 15;
clim(ax, [lmin, lmax]);

% ── 9c. Nullclines ───────────────────────────────────────
if show_nullclines
    contour(ax, Xnc, Ync, NC1, [0 0], 'r-',  'LineWidth', 2.0);
    contour(ax, Xnc, Ync, NC2, [0 0], 'b--', 'LineWidth', 2.0);
    % Legend patches (invisible plots for legend)
    plot(ax, NaN, NaN, 'r-',  'LineWidth', 2.0, 'DisplayName', 'x-nullcline  (dx/dt=0)');
    plot(ax, NaN, NaN, 'b--', 'LineWidth', 2.0, 'DisplayName', 'y-nullcline  (dy/dt=0)');
end

% ── 9d. Trajectories ─────────────────────────────────────
if show_trajectories
    traj_cmap = lines(size(IC,1));
    for k = 1:size(IC,1)
        S = trajs{k};
        if isempty(S), continue; end
        % Clip to window
        in_win = S(:,1) >= x_range(1) & S(:,1) <= x_range(2) & ...
                 S(:,2) >= y_range(1) & S(:,2) <= y_range(2);
        S = S(in_win, :);
        if size(S,1) < 2, continue; end
        plot(ax, S(:,1), S(:,2), '-', 'Color', [traj_cmap(k,:), 0.85], ...
             'LineWidth', 1.8, 'HandleVisibility','off');
        % Start marker
        plot(ax, S(1,1), S(1,2), 'o', ...
             'Color', traj_cmap(k,:), 'MarkerFaceColor', traj_cmap(k,:), ...
             'MarkerSize', 7, 'HandleVisibility','off');
        % Direction arrow mid-trajectory
        mid = max(1, round(size(S,1)*0.45));
        if mid < size(S,1)
            dx = S(mid+1,1)-S(mid,1);  dy = S(mid+1,2)-S(mid,2);
            quiver(ax, S(mid,1), S(mid,2), dx, dy, 0, ...
                   'Color', traj_cmap(k,:), 'LineWidth', 1.5, ...
                   'MaxHeadSize', 3, 'HandleVisibility','off');
        end
    end
end

% ── 9e. Equilibria ───────────────────────────────────────
if show_equilibria && ~isempty(eq_points)
    J_sym = jacobian(f_sym, [x, y]);
    for k = 1:size(eq_points,1)
        px = eq_points(k,1);  py = eq_points(k,2);
        J_num = double(subs(J_sym, [x,y,param_syms], [px,py,param_vals]));
        ev    = eig(J_num);
        re    = real(ev);
        im    = imag(ev);
        if max(re) < -1e-8
            if max(abs(im))>1e-8, marker='o'; fc=[0.2 0.6 1.0]; label='stable spiral';
            else,                  marker='o'; fc=[0.0 0.4 0.8]; label='stable node'; end
        elseif min(re) < -1e-8
            marker = 's'; fc = [0.9 0.5 0.1]; label = 'saddle';
        elseif max(re) > 1e-8
            if max(abs(im))>1e-8, marker='o'; fc=[1.0 0.3 0.3]; label='unstable spiral';
            else,                  marker='o'; fc=[0.8 0.1 0.1]; label='unstable node'; end
        else
            marker = 'd'; fc = [0.5 0.5 0.5]; label = 'non-hyperbolic';
        end
        plot(ax, px, py, marker, 'MarkerFaceColor', fc, ...
             'MarkerEdgeColor', 'k', 'MarkerSize', 11, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('E%d (%s)', k, label));
        text(ax, px + 0.05*diff(x_range), py + 0.03*diff(y_range), ...
             sprintf('E_%d', k), 'FontSize', 12, 'FontWeight', 'bold', ...
             'Color', 'k');
    end
end

% ── 9f. Formatting ───────────────────────────────────────
xlim(ax, x_range);  ylim(ax, y_range);
xlabel(ax, 'H', 'FontSize', 20);
ylabel(ax, 'F', 'FontSize', 20);

% % Build title string
% param_str = '';
% for k = 1:length(param_syms)
%     param_str = [param_str, sprintf('%s=%g', char(param_syms(k)), param_vals(k))]; %#ok<AGROW>
%     if k < length(param_syms), param_str = [param_str, ',  ']; end %#ok<AGROW>
% end
% if isempty(param_str)
%     title(ax, sprintf('Phase Plane:  dx/dt = %s,   dy/dt = %s', ...
%           char(f1_sub), char(f2_sub)), 'FontSize', 11, 'FontWeight','bold');
% else
%     title(ax, sprintf('Phase Plane  [%s]\ndx/dt = %s,   dy/dt = %s', ...
%           param_str, char(f1_sub), char(f2_sub)), ...
%           'FontSize', 10, 'FontWeight','bold');
% end

% if show_nullclines || show_equilibria
%     legend(ax, 'Location','best','FontSize',9,'Box','on');
% end

grid(ax,'on'); box(ax,'on'); ax.GridAlpha = 0.18;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

fprintf('Done. Phase plane displayed.\n');

%% ═══════════════════════════════════════════════════════
%  LOCAL FUNCTION — stop integration at plot boundary
% ═══════════════════════════════════════════════════════
function [value, isterminal, direction] = exit_event(~, X, xr, yr)
    % Triggers when trajectory leaves the plot window
    dist = min([X(1)-xr(1), xr(2)-X(1), X(2)-yr(1), yr(2)-X(2)]);
    value      = dist;
    isterminal = 1;
    direction  = -1;
end
