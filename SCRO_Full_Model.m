function SCRO_Full_Model()
% =========================================================================
% SEMI-CLOSED REVERSE OSMOSIS (SCRO) OPTIMIZATION MODEL
% =========================================================================

    clear variables; clc; close all;

    % --- 1. SYSTEM PARAMETERS ---
    geo.L = 6.0;                     % Length [m]
    geo.W = 35.0;                    % Width [m]
    geo.H = 0.7e-3;                  % Channel height [m]
    geo.epsilon = 0.85;              % Porosity
    geo.dh = 4 * geo.epsilon * geo.H / (2 + 2*geo.epsilon);      % Hydraulic diameter
    geo.Area = geo.L * geo.W;        % Membrane area

    mem.A = 3.0e-12;                 % Water Permeability (~1.1 LMH/bar)
    mem.B = 1.0e-7;                  % Salt Permeability [m/s]

    % Efficiencies & Losses
    sys.eff_pump = 0.80;             % High-pressure pump efficiency
    sys.eff_erd  = 0.98;             % ERD efficiency
    sys.dP_loss_circ = 0.5;          % Circulation pressure drop penalty [bar]

    % Feed Conditions
    feed_init.C_gL = 35.0;           % Feed salinity [g/L]
    feed_init.Temp_K = 298.15;       % Temperature [K]
    feed_init.Vol = 10.0;            % Batch volume [m^3]

    % Osmotic Pressure
    props_f = get_fluid_properties(feed_init.C_gL, feed_init.Temp_K);
    Pi_feed_bar = props_f.osmotic_pressure_bar;

    % =========================================================================
    % PHASE 1: EFFICIENCY ANALYSIS (Variable Cycles)
    % =========================================================================
    fprintf('\n--- Phase 1: Efficiency Analysis (Fixed RR=0.60) ---\n');
    RR_fixed = 0.60;

    N_range = 1:6;

    % Thermodynamic Minimum SEC
    SEC_ideal_fixed = (Pi_feed_bar * 1e5 / RR_fixed * log(1/(1-RR_fixed))) / 3.6e6;

    results_p1.N = N_range;
    results_p1.SEC = nan(length(N_range), 1);
    results_p1.Efficiency = nan(length(N_range), 1);

    for i = 1:length(N_range)
        N = N_range(i);
        SEC = run_scro_batch(N, RR_fixed, feed_init, geo, mem, sys);
        results_p1.SEC(i) = SEC;

        if ~isnan(SEC)
            results_p1.Efficiency(i) = (SEC_ideal_fixed / SEC) * 100;
            fprintf('N=%d | SEC=%.3f | Eff=%.1f%%\n', N, SEC, results_p1.Efficiency(i));
        end
    end

    % =========================================================================
    % PHASE 2: RECOVERY SWEEP
    % =========================================================================
    fprintf('\n--- Phase 2: Recovery Sweep ---\n');

    RR_sweep = 0.30:0.05:0.85;
    N_selected = 1:6;

    results_p2.RR = RR_sweep;
    results_p2.Ideal_Curve = zeros(length(RR_sweep), 1);
    results_p2.Model_Curves = nan(length(RR_sweep), length(N_selected));

    % 1. Ideal Curve Calculation
    for j = 1:length(RR_sweep)
        rr = RR_sweep(j);
        val_ideal_J = (Pi_feed_bar * 1e5 / rr) * log(1 / (1 - rr));
        results_p2.Ideal_Curve(j) = val_ideal_J / 3.6e6;
    end

    % 2. Model Simulation
    for k = 1:length(N_selected)
        N_curr = N_selected(k);
        fprintf('Simulating N=%d... ', N_curr);
        for j = 1:length(RR_sweep)
            rr = RR_sweep(j);
            results_p2.Model_Curves(j, k) = run_scro_batch(N_curr, rr, feed_init, geo, mem, sys);
        end
        fprintf('Done.\n');
    end

    % =========================================================================
    % PLOTTING RESULTS
    % =========================================================================

    % --- GRAPH 1: SEC vs Recovery Ratio ---
    figure(1); clf; hold on; box on;

    legend_handles = [];
    legend_strings = {};

    % Plot Ideal Limit
    h_ideal = plot(results_p2.RR, results_p2.Ideal_Curve, '--k', 'LineWidth', 1);
    legend_handles(end+1) = h_ideal;
    legend_strings{end+1} = 'Thermodynamic Limit';

    % Plot Simulation Results
    colors = jet(length(N_selected));
    for k = 1:length(N_selected)
        valid = ~isnan(results_p2.Model_Curves(:, k));
        if any(valid)
            % Determine line color
            if N_selected(k) == 1
                line_color = [1, 0, 0];         % Red for SSRO
            elseif N_selected(k) == 5
                line_color = [0.85, 0.65, 0.1]; % Gold
            else
                line_color = colors(k,:);
            end

            h = plot(results_p2.RR(valid), results_p2.Model_Curves(valid, k), '-o', ...
                'Color', line_color, 'LineWidth', 1.2, 'MarkerSize', 4);

            legend_handles(end+1) = h;

            if N_selected(k) == 1
                legend_strings{end+1} = 'SSRO (N=1)';
            else
                legend_strings{end+1} = sprintf('SCRO N=%d', N_selected(k));
            end
        end
    end

    xlabel('Water Recovery Ratio (RR)');
    ylabel('Specific Energy Consumption (kWh/m^3)');
    title('SEC vs. Recovery Ratio');
    legend(legend_handles, legend_strings, 'Location', 'NorthWest', 'NumColumns', 2);
    grid on;

    % Axes Formatting
    ylim([0.0, 6.0]);
    xlim([0.3, 0.8]);
    xticks(0.3:0.05:0.8);
    yticks(0:1.0:6.0);
    set(gca, 'FontSize', 12);

    % --- GRAPH 2: SEC vs Number of Cycles ---
    figure(2); clf; hold on; box on;

    valid = ~isnan(results_p1.SEC);
    h_sec = plot(results_p1.N(valid), results_p1.SEC(valid), '-s', ...
        'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'b');

    [min_sec, min_idx] = min(results_p1.SEC);
    h_min = plot(results_p1.N(min_idx), min_sec, 'rp', 'MarkerSize', 12, ...
        'MarkerFaceColor', 'r');

    xlabel('Number of Cycles (N)');
    ylabel('Specific Energy Consumption (kWh/m^3)');
    title(sprintf('SEC vs Number of Cycles (RR=%.2f)', RR_fixed));
    legend([h_sec, h_min], {'SCRO SEC', 'Optimal Point'}, 'Location', 'NorthEast');
    grid on;
    xticks(results_p1.N);
    set(gca, 'FontSize', 12);

    % --- GRAPH 3: Second-Law Efficiency ---
    figure(3); clf; hold on; box on;

    valid = ~isnan(results_p1.Efficiency);
    plot(results_p1.N(valid), results_p1.Efficiency(valid), '-s', ...
        'LineWidth', 2, 'Color', [0, 0.5, 0], 'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.8, 0.2]);

    line([0.5, 6.5], [100, 100], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
    text(1, 102, 'Max Thermodynamic Efficiency (100%)', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 10);

    xlabel('Number of Cycles (N)');
    ylabel('Second-Law Efficiency (%)');
    title('Second-Law Efficiency vs Cycle Count');
    grid on;
    xticks(results_p1.N);
    ylim([0, 110]);
    xlim([0.5, 6.5]);
    set(gca, 'FontSize', 12);
end

% =========================================================================
% Helper Functions
% =========================================================================

function SEC = run_scro_batch(N, RR_target, feed_initial, geo, mem, sys)
    r_cycle = 1 - (1 - RR_target)^(1/N);
    current_feed = feed_initial;
    total_net_energy_J = 0;
    total_permeate_vol = 0;

    for cyc = 1:N
        target_ret_vol = current_feed.Vol * (1 - r_cycle);
        props = get_fluid_properties(current_feed.C_gL, current_feed.Temp_K);
        P_guess = (props.osmotic_pressure_bar / (1 - r_cycle)) * 1.1 + 5;

        [P_applied, ret_state] = solve_pressure(current_feed, target_ret_vol, P_guess, geo, mem);

        if isnan(P_applied)
            SEC = NaN; return;
        end

        Hyd_Work_In = (P_applied + sys.dP_loss_circ) * 1e5 * current_feed.Vol;
        Hyd_Work_Rec = (ret_state.P_out * 1e5) * ret_state.Vol * sys.eff_erd;
        Net_Energy_cycle = (Hyd_Work_In - Hyd_Work_Rec) / sys.eff_pump;

        total_net_energy_J = total_net_energy_J + Net_Energy_cycle;
        vol_perm = current_feed.Vol - ret_state.Vol;
        total_permeate_vol = total_permeate_vol + vol_perm;

        current_feed.C_gL = ret_state.C_gL;
        current_feed.Vol  = ret_state.Vol;

        if current_feed.Vol < 0.01; break; end
    end

    if total_permeate_vol > 0
        SEC = (total_net_energy_J / 3.6e6) / total_permeate_vol;
    else
        SEC = NaN;
    end
end

function [P_res, retentate] = solve_pressure(feed, target_vol, P_guess_init, geo, mem)
    P_min = 1.0;
    P_max = 300.0;
    tol = 1e-3;

    P_res = NaN; retentate.Vol = NaN; retentate.C_gL = NaN; retentate.P_out = NaN;

    if P_guess_init > P_max; P_guess_init = P_max; end
    [v_max, ~, ~] = run_module_1D(feed, P_max, geo, mem);

    if v_max > target_vol; return; end

    P_low = P_min; P_high = P_max;

    for k = 1:40
        P_mid = (P_low + P_high)/2;
        [ret_vol, ret_C, ret_P] = run_module_1D(feed, P_mid, geo, mem);

        err = ret_vol - target_vol;
        if abs(err/target_vol) < tol
            P_res = P_mid;
            retentate.Vol = ret_vol;
            retentate.C_gL = ret_C;
            retentate.P_out = ret_P;
            return;
        end

        if err > 0; P_low = P_mid; else; P_high = P_mid; end
    end

    P_res = P_mid;
    retentate.Vol = ret_vol;
    retentate.C_gL = ret_C;
    retentate.P_out = ret_P;
end

function [vol_out, C_out, P_out, Cp_bulk_avg] = run_module_1D(feed, P_in, geo, mem)
    N_nodes = 20;
    dz = geo.L / N_nodes;

    chan_area = geo.W * geo.H * geo.epsilon;
    v_in = 0.15;
    Q = v_in * chan_area;

    C = feed.C_gL;
    M_salt_flow = Q * C;
    P = P_in;

    total_salt_perm_kg_s = 0;
    total_water_perm_m3_s = 0;

    for i = 1:N_nodes
        props = get_fluid_properties(C, feed.Temp_K);
        v = Q / chan_area;
        Re = (props.rho * v * geo.dh) / props.mu;
        Sc = props.mu / (props.rho * props.D);
        Sh = 0.46 * (Re * Sc)^0.36;
        k_mt = (props.D * Sh) / geo.dh;

        pi_bulk = props.osmotic_pressure_bar;
        Jw = 0;
        for iter = 1:10
             cp_factor = exp(min(5, Jw / k_mt));
             Cm = C * cp_factor;
             pi_m = pi_bulk * (Cm / C);
             NDP = P - pi_m;
             if NDP < 0; NDP = 0; end
             Jw_new = mem.A * NDP * 1e5;
             if abs(Jw_new - Jw) < 1e-7
                 Jw = Jw_new;
                 break;
             end
             Jw = 0.5*Jw + 0.5*Jw_new;
        end

        if Cm > 260
            Cm = 260;
            pi_m = pi_bulk * (260 / C);
            NDP = P - pi_m;
            if NDP < 0; NDP = 0; end
            Jw = mem.A * NDP * 1e5;
        end

        Js = Cm * (mem.B * Jw) / (Jw + mem.B);
        f = 0.42 + 189.3/(Re+0.1);
        dP = (f * props.rho * v^2 * dz) / (2 * geo.dh);

        Q_perm = Jw * geo.W * dz;
        M_salt_perm_node = Js * geo.W * dz;

        Q = Q - Q_perm;
        if Q < 1e-6; Q = 1e-6; end
        M_salt_flow = M_salt_flow - M_salt_perm_node;

        C = M_salt_flow / Q;
        P = P - dP/1e5;

        total_salt_perm_kg_s = total_salt_perm_kg_s + M_salt_perm_node;
        total_water_perm_m3_s = total_water_perm_m3_s + Q_perm;
    end

    vol_out = feed.Vol * (Q / (v_in * chan_area));
    C_out = C;
    P_out = P;
    if total_water_perm_m3_s > 0
        Cp_bulk_avg = total_salt_perm_kg_s / total_water_perm_m3_s;
    else
        Cp_bulk_avg = 0;
    end
end

function props = get_fluid_properties(C_gL, ~)
    rho = 1000 + 0.7*C_gL;
    X = C_gL / rho;
    mu = (2.15e-3 * X) + 9.80e-4;
    pi_bar = 0.848 * (3.14e-6 * C_gL^2 + 2.13e-4 * C_gL + 0.917) * C_gL;
    props.rho = rho; props.mu = mu; props.osmotic_pressure_bar = pi_bar; props.D = 1.5e-9;
end
