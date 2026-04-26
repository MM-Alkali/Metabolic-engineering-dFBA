%% ------------
%% MAIN SCRIPT 
%% ------------

clc; clear; close all;

% Set paths
dataPath = 'C:\Users\muawi\Documents\MatLab\Metabolic-engineering-dFBA';
outDir = fullfile(dataPath, 'Results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Load model
modelFile = fullfile(dataPath, 'yeast-GEM.mat');
if exist(modelFile, 'file')
    tmp = load(modelFile);
    flds = fieldnames(tmp);
    model = tmp.(flds{1});
else
    error('Model file not found');
end

% Set solver
changeCobraSolver('glpk', 'LP', 0);

% Biomass reaction
biomassRxnID = 'r_2111';

% Conditions
conditions = {
    'Anaerobic_S_P', ...
    'Anaerobic_R_P', ...
    'Aerobic_S_P', ...
    'Aerobic_R_P'};

% Reaction mapping
fluxMap.qs = 'r_1714';     % Glucose
fluxMap.qco2 = 'r_1672';   % CO2
fluxMap.qo2 = 'r_1992';    % O2
fluxMap.qace = 'r_1634';   % Acetate
fluxMap.qcit = 'r_1687';   % Citrate
fluxMap.qlac = 'r_1546';   % Lactate
fluxMap.qpyr = 'r_2033';   % Pyruvate
fluxMap.qeth = 'r_1761';   % Ethanol
fluxMap.qgly = 'r_1808';   % Glycerol

% Pyruvate metabolism reactions
fluxMap.pyk = 'r_0962';    % Pyruvate kinase (PEP → pyruvate) - CDC19
fluxMap.pdc = 'r_0959';    % Pyruvate decarboxylase (main PDC)
fluxMap.pdh = 'r_0961';    % Pyruvate dehydrogenase (mitochondrial)

% Run dFBA for each condition
results = struct();
for c = 1:numel(conditions)
    condName = conditions{c};
    fprintf('\n=== %s ===\n', condName);
    
    xlsFile = fullfile(dataPath, [condName '.xlsx']);
    if ~exist(xlsFile, 'file')
        warning('File %s not found', xlsFile);
        continue;
    end
    
    % Read data
    T = readtable(xlsFile, 'Sheet', 'mmol');
    
    % Run dFBA
    res = run_dFBA(model, T, fluxMap, biomassRxnID, condName);
    results.(condName) = res;
end

%% 1a: MAIN COMPARISON SP VS RP WITH ALL BYPRODUCTS (NO LEGENDS)
fprintf('\n=== Generating Figure 1: Main Comparison SP vs RP ===\n');
createMainComparisonPlot(results, fluxMap, biomassRxnID, outDir);

%% 1b: LEGEND FOR MAIN COMPARISON PLOT
fprintf('\n=== Generating Legend for Main Comparison ===\n');
createMainLegend(outDir);

%% 1c: LEGEND FOR PYRUVATE ROUTING ANALYSIS
fprintf('\n=== Generating Legend for Pyruvate Routing ===\n');
createPyruvateRoutingLegend(outDir);

%% 1d: SUPPLEMENTARY FIGURE: Model Validation (Experimental vs Model)
fprintf('\n=== Generating Supplementary Figures: Model Validation ===\n');
createValidationFigures(results, fluxMap, biomassRxnID, outDir);

%% 2a: SENSITIVITY ANALYSIS
fprintf('\n=== Running PYK Sensitivity Analysis (Generating CSV data) ===\n');
for c = 1:numel(conditions)
    condName = conditions{c};
    xlsFile = fullfile(dataPath, [condName '.xlsx']);
    if exist(xlsFile, 'file')
        T_sens = readtable(xlsFile, 'Sheet', 'mmol');
        analyzePyruvateKinaseSensitivity(model, T_sens, fluxMap, biomassRxnID, condName, outDir);
    end
end

%% 2b: CREATE COMBINED FIGURES FROM CSV DATA 
fprintf('\n=== Creating Combined Sensitivity Figures ===\n');
analyzeCombinedSensitivity(outDir);

%% 3: FLUX VARIABILITY ANALYSIS (Combined on one plot)
fprintf('\n=== Performing Flux Variability Analysis ===\n');
analyzeFluxVariabilityCombined(model, T, fluxMap, biomassRxnID, outDir);

%% --------------------------------
%% 0: dFBA FUNCTION - WITH CONSTRAINTS
%% --------------------------------
function res = run_dFBA(model0, T, fluxMap, biomassRxnID, condName)

    time_s = T.Time;
    nT = height(T);
    nRxn = numel(model0.rxns);
    
    fluxMatrix = zeros(nT, nRxn);
    mu_predicted = zeros(nT, 1);
    status = cell(nT, 1);
    
    % Detect anaerobic/aerobic from condition name
    isAnaerobic = contains(lower(condName), 'anaerobic');
    
    % Feed phase detection (high glucose uptake > 5 mmol/gCDW/h)
    feed_phase = false(nT, 1);
    if ismember('qs', T.Properties.VariableNames)
        feed_phase = abs(T.qs) > 5;
    end
    
    % Loop through time points
    for k = 1:nT
        model = model0;
        
        %% 1. MAXIMISE BIOMASS
        model.c(:) = 0;
        biomassIdx = findRxnIDs(model, biomassRxnID);
        model.c(biomassIdx) = 1;
        
        %% 2. ATP MAINTENANCE
        atpm_idx = findRxnIDs(model, 'r_4046');
        if atpm_idx > 0
            model.lb(atpm_idx) = 1.5;
            model.ub(atpm_idx) = 2.5;
        end
        
        %% 3. GROWTH RATE LIMITS
        if isAnaerobic
            model.ub(biomassIdx) = 0.25;
            model.lb(biomassIdx) = 0;
            
            % Allow anaerobic supplements
            erg_idx = findRxnIDs(model, 'r_1897');
            if erg_idx > 0
                model.lb(erg_idx) = -0.005;
                model.ub(erg_idx) = 0;
            end
            
            oleate_idx = findRxnIDs(model, 'r_2134');
            if oleate_idx > 0
                model.lb(oleate_idx) = -0.002;
                model.ub(oleate_idx) = 0;
            end
        else
            model.ub(biomassIdx) = 0.45;
            model.lb(biomassIdx) = 0;
        end
        
        if feed_phase(k)
            if isAnaerobic
                model.ub(biomassIdx) = 0.25;
            else
                model.ub(biomassIdx) = 0.45;
            end
        else
            model.ub(biomassIdx) = 0.05;
        end
        
        %% 4. CONSTRAIN GLUCOSE UPTAKE
        if ismember('qs', T.Properties.VariableNames)
            qs = T.qs(k);
            idx = findRxnIDs(model, fluxMap.qs);
            if idx > 0
                model.lb(idx) = qs;
                model.ub(idx) = 0;
            end
        end
        
        %% 5. BLOCK O2 FOR ANAEROBIC
        if isAnaerobic && isfield(fluxMap, 'qo2')
            o2_idx = findRxnIDs(model, fluxMap.qo2);
            if o2_idx > 0
                model.lb(o2_idx) = 0;
                model.ub(o2_idx) = 0;
            end
        end
        
        %% 6. PYRUVATE KINASE CONSTRAINT
        if isfield(fluxMap, 'pyk')
            pyk_idx = findRxnIDs(model, fluxMap.pyk);
            if pyk_idx > 0
                q_glc_abs = abs(T.qs(k));
                
                if feed_phase(k)  % During glucose pulse
                    if isAnaerobic
                        % ANAEROBIC: Keep at 60% 
                        normal_pyk = q_glc_abs * 2.0;
                        constraint_level = 0.6;  % 60% 
                        model.ub(pyk_idx) = normal_pyk * constraint_level;
                        fprintf('PYK=%.1f(60%%) ', model.ub(pyk_idx));
                    else
                        % AEROBIC: NO CONSTRAINT 
                        % Let the model use full PYK capacity
                        model.ub(pyk_idx) = 1000;  % Unconstrained
                        fprintf('PYK=unconstrained ');
                    end
                else
                    % POST-FEED: No constraint for both conditions
                    model.ub(pyk_idx) = 1000;  % Unconstrained
                end
                model.lb(pyk_idx) = 0;
            end
        end
        
        %% 7. Solve
        sol = optimizeCbModel(model, 'max');
        
        if sol.stat == 1
            mu_predicted(k) = sol.v(biomassIdx);
            fluxMatrix(k,:) = sol.v(:)';
            status{k} = 'optimal';
            
            % DEBUG: Print actual PYK and PDC fluxes for verification
            pyk_idx = findRxnIDs(model, fluxMap.pyk);
            pdc_idx = findRxnIDs(model, fluxMap.pdc);
            
            if pyk_idx > 0
                fprintf(' Actual PYK=%.1f', sol.v(pyk_idx));
            end
            if pdc_idx > 0
                fprintf(' PDC=%.1f', sol.v(pdc_idx));
                if pyk_idx > 0 && sol.v(pyk_idx) > 0
                    fprintf(' Ratio=%.0f%%', sol.v(pdc_idx)/sol.v(pyk_idx)*100);
                end
            end
            fprintf('\n');
            
        else
            status{k} = 'infeasible';
        end
    end
    
    %% Package results
    res = struct();
    res.condition = condName;
    res.time_s = time_s;
    res.mu_predicted = mu_predicted;
    res.fluxMatrix = fluxMatrix;
    res.rxnIDs = model0.rxns;
    res.rxnNames = model0.rxnNames;
    res.status = status;
    res.feed_phase = feed_phase;
    res.isAnaerobic = isAnaerobic;
    
    % Store experimental data
    res.measured = struct();
    res.measured.time = time_s;
    res.measured.qs = T.qs;
    
    if ismember('mu', T.Properties.VariableNames)
        res.measured.mu = T.mu;
    end
    if ismember('qco2', T.Properties.VariableNames)
        res.measured.qco2 = T.qco2;
    end
    if ismember('qo2', T.Properties.VariableNames)
        res.measured.qo2 = T.qo2;
    end
    
    products = {'qeth','qgly','qace','qcit','qlac','qpyr'};
    for i = 1:length(products)
        if ismember(products{i}, T.Properties.VariableNames)
            res.measured.(products{i}) = T.(products{i});
        end
    end
end

%% -------------------------------------------------------------
%% 1a: MAIN COMPARISON + PYRUVATE ROUTING ANALYSIS
%% -------------------------------------------------------------
function createMainComparisonPlot(results, fluxMap, biomassRxnID, outDir)
    
    % Define colors
    colors = struct();
    colors.Aerobic_S_P = [0.6294, 0.8078, 0.9804]; % Light blue
    colors.Aerobic_R_P = [0, 0, 0.6451];            % Dark blue
    colors.Anaerobic_S_P = [1.0000, 0.6471, 0];     % Light orange
    colors.Anaerobic_R_P = [0.8000, 0.5000, 0];     % Dark orange
    
    % Order of conditions for plotting
    cond_order = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};
    
    % Parameters to plot (name, rxnID, ylabel, invert)
    % First 8 panels: main comparison
    params = {
        'Glucose uptake', fluxMap.qs, '[mmol/gDW/h]', true;
        'Growth rate', biomassRxnID, '[h^{-1}]', false;
        'CO_2 production', fluxMap.qco2, '[mmol/gDW/h]', false;
        'O_2 uptake', fluxMap.qo2, '[mmol/gDW/h]', true;
        'Ethanol', fluxMap.qeth, '[mmol/gDW/h]', false;
        'Pyruvate', fluxMap.qpyr, '[mmol/gDW/h]', false; 
        'Pyruvate kinase (PYK)', fluxMap.pyk, '[mmol/gDW/h]', false;
        'PDC & PDH routing', '', '[mmol/gDW/h]', false;
    };
    
    % Set font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
    
    % Create figure with dimensions (3 rows x 4 columns)
    fig = figure('Position', [50, 50, 1800, 1200], 'Color', 'white');
    
    % Get time points and feed phase from first condition
    first_cond = cond_order{1};
    if ~isfield(results, first_cond)
        first_cond = fieldnames(results);
        first_cond = first_cond{1};
    end
    time_s = results.(first_cond).time_s;
    feed_phase = results.(first_cond).feed_phase;
    
    %% ROWS 1-2: Main Comparison (8 panels)
    for p = 1:8
        subplot(3, 4, p);
        hold on;
        box on;
        
        param_name = params{p, 1};
        rxn_id = params{p, 2};
        ylabel_str = params{p, 3};
        invert_plot = params{p, 4};
        
        % Special handling for last panel of main comparison (PDC & PDH)
        if p == 8
            % Plot PDC and PDH for all conditions
            for c = 1:length(cond_order)
                condName = cond_order{c};
                if ~isfield(results, condName)
                    continue;
                end
                res = results.(condName);
                
                % Determine line style based on condition
                if contains(condName, '_S_')
                    line_style = '-';  % SP - solid
                else
                    line_style = '--'; % RP - broken
                end
                
                % Plot PDC
                pdc_idx = findRxnIDs(res.rxnIDs, fluxMap.pdc);
                if pdc_idx > 0
                    pdc_flux = res.fluxMatrix(:, pdc_idx);
                    plot(time_s, pdc_flux, 'Color', colors.(condName), ...
                        'LineStyle', line_style, 'LineWidth', 1.5);
                end
                
                % Plot PDH (dotted line to distinguish from PDC)
                pdh_idx = findRxnIDs(res.rxnIDs, fluxMap.pdh);
                if pdh_idx > 0
                    pdh_flux = res.fluxMatrix(:, pdh_idx);
                    plot(time_s, pdh_flux, 'Color', colors.(condName), ...
                        'LineStyle', ':', 'LineWidth', 1.5);
                end
            end
            
        else
            % Regular parameter - plot all conditions
            for c = 1:length(cond_order)
                condName = cond_order{c};
                if ~isfield(results, condName)
                    continue;
                end
                res = results.(condName);
                
                idx = findRxnIDs(res.rxnIDs, rxn_id);
                if idx > 0
                    flux = res.fluxMatrix(:, idx);
                    if invert_plot
                        flux = -flux;
                    end
                    
                    % Determine line style
                    if contains(condName, '_S_')
                        line_style = '-';  % SP - solid
                    else
                        line_style = '--'; % RP - broken
                    end
                    
                    plot(time_s, flux, 'Color', colors.(condName), ...
                        'LineStyle', line_style, 'LineWidth', 2);
                end
            end
        end
        
        % Set x-axis limits
        xlim([0 150]);
        
        % Calculate y-limits with padding
        children = get(gca, 'Children');
        y_min = inf;
        y_max = -inf;
        for i = 1:length(children)
            if strcmp(get(children(i), 'Type'), 'line')
                ydata = get(children(i), 'YData');
                y_min = min(y_min, min(ydata));
                y_max = max(y_max, max(ydata));
            end
        end
        
        if isfinite(y_min) && isfinite(y_max)
            y_padding = 0.1 * (y_max - y_min);
            if y_padding < 0.01
                y_padding = 0.1;
            end
            ylim([y_min - y_padding, y_max + y_padding]);
        end
        
        % Highlight feed phase
        if any(feed_phase)
            yl = ylim;
            feed_end = max(time_s(feed_phase));
            fill([0 feed_end feed_end 0], [yl(1) yl(1) yl(2) yl(2)], ...
                [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
            line([feed_end feed_end], yl, 'Color', [0.5 0.5 0.5], ...
                'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
        end
        
        xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel(ylabel_str, 'FontSize', 10, 'FontWeight', 'bold');
        title(param_name, 'FontSize', 14, 'FontWeight', 'bold');
        
        grid on;
    end
    
    %% ROW 3: Pyruvate Routing Analysis (4 panels)
    subplot_order = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};

    for c = 1:4
        condName = subplot_order{c};
        panel_idx = 8 + c;  % Panels 9-12
        
        subplot(3, 4, panel_idx);
        hold on;
        box on;
        
        if ~isfield(results, condName)
            text(0.5, 0.5, sprintf('%s\nData not available', condName), ...
                'HorizontalAlignment', 'center', 'FontSize', 12, ...
                'FontName', 'Helvetica', 'FontWeight', 'bold');
            axis off;
            continue;
        end
        
        res = results.(condName);
        
        % Get indices for pyruvate-related reactions
        pyk_idx = findRxnIDs(res.rxnIDs, fluxMap.pyk);
        pdc_idx = findRxnIDs(res.rxnIDs, fluxMap.pdc);
        pdh_idx = findRxnIDs(res.rxnIDs, fluxMap.pdh);
        
        time_s = res.time_s;
        fluxMatrix = res.fluxMatrix;
        feed_phase = res.feed_phase;
        
        % Get fluxes
        pyk_flux = fluxMatrix(:, pyk_idx);
        pdc_flux = zeros(length(time_s), 1);
        pdh_flux = zeros(length(time_s), 1);
        
        if pdc_idx > 0
            pdc_flux = fluxMatrix(:, pdc_idx);
        end
        
        if pdh_idx > 0
            pdh_flux = fluxMatrix(:, pdh_idx);
        end
        
        % Plot fluxes with condition-specific colors
        plot(time_s, pyk_flux, 'Color', colors.(condName), 'LineStyle', '-', ...
            'LineWidth', 2);
        plot(time_s, pdc_flux, 'Color', colors.(condName), 'LineStyle', '--', ...
            'LineWidth', 2);
        plot(time_s, pdh_flux, 'Color', colors.(condName), 'LineStyle', ':', ...
            'LineWidth', 2);
        
        % Set x-axis limits
        xlim([0 150]);
        
        % Calculate y-limits with padding
        all_flux_data = [pyk_flux(:); pdc_flux(:); pdh_flux(:)];
        y_min = min(all_flux_data);
        y_max = max(all_flux_data);
        
        if y_max > y_min
            y_padding = 0.15 * (y_max - y_min);
            y_min = y_min - y_padding;
            y_max = y_max + y_padding;
        else
            y_min = y_min - 1;
            y_max = y_max + 1;
        end
        
        % Highlight feed phase
        if any(feed_phase)
            feed_end = max(time_s(feed_phase));
            fill([0 feed_end feed_end 0], [y_min y_min y_max y_max], ...
                 [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
            line([feed_end feed_end], [y_min y_max], ...
                 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
                 'LineStyle', '--', 'HandleVisibility', 'off');
        end
        
        ylim([y_min y_max]);
        
        xlabel('Time [s]', 'FontSize', 10, 'FontWeight', 'bold');
        ylabel('Flux [mmol/gDW/h]', 'FontSize', 10, 'FontWeight', 'bold');
        
        % Format condition name for title
        title_cond = strrep(condName, '_', ' ');
        title(title_cond, 'FontSize', 14, 'FontWeight', 'bold');
        
        grid on;
        
        % Calculate and display statistics in top right corner
        if any(feed_phase)
            mean_pyk = mean(pyk_flux(feed_phase));
            mean_pdc = mean(pdc_flux(feed_phase));
            mean_pdh = mean(pdh_flux(feed_phase));
            
            if mean_pyk > 0
                pdc_percent = mean_pdc/mean_pyk*100;
                pdh_percent = mean_pdh/mean_pyk*100;
            else
                pdc_percent = 0;
                pdh_percent = 0;
            end
            
            text_str = {
                sprintf('PYK: %.1f', mean_pyk);
                sprintf('PDC: %.1f (%.0f%%)', mean_pdc, pdc_percent);
                sprintf('PDH: %.1f (%.0f%%)', mean_pdh, pdh_percent)
            };
            
            % MOVED to TOP RIGHT corner
            text(0.9, 0.9, text_str, 'Units', 'normalized', ...
                 'BackgroundColor', 'white', 'EdgeColor', [0.5 0.5 0.5], ...
                 'FontSize', 9, 'FontName', 'Helvetica', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
    end
        
    % Add overall title
    sgtitle('Figure 1: Model Predictions and Pyruvate Routing Analysis under Glucose Pulsing', ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
    
    % Save figure
    pngFile = fullfile(outDir, 'Main_Comparison_SP_RP.png');
    exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Main Comparison + Pyruvate Routing Plot: %s\n', pngFile);
    
    %% Reset default font
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
end

%% ------------------------------------------------------------------------
%% 1b: LEGEND FOR MAIN COMPARISON PLOT
%% ------------------------------------------------------------------------
function createMainLegend(outDir)
    % Creates a separate figure with legend for the main comparison plot
    
    figure('Position', [100, 100, 400, 300], 'Color', 'white');
    axis off;
    
    % Create dummy lines for legend
    hold on;
    
    % Aerobic SP
    plot(nan, nan, 'Color', [0.6294, 0.8078, 0.9804], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic SP');
    
    % Aerobic RP
    plot(nan, nan, 'Color', [0, 0, 0.6451], 'LineStyle', '--', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic RP');
    
    % Anaerobic SP
    plot(nan, nan, 'Color', [1.0000, 0.6471, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic SP');
    
    % Anaerobic RP
    plot(nan, nan, 'Color', [0.8000, 0.5000, 0], 'LineStyle', '--', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic RP');
    
    % Feed phase indicator
    fill([0 1 1 0], [0 0 1 1], [0.8 0.8 0.8], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName', 'Feed Phase');
    
    legend('Location', 'best', 'FontSize', 10, 'FontName', 'Helvetica', ...
        'FontWeight', 'bold', 'Box', 'on');
    
    title('Legend: Main Comparison Plot', 'FontSize', 14, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Save legend
    pngFile = fullfile(outDir, 'Legend_Main_Comparison.png');
    exportgraphics(gcf, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Main Legend: %s\n', pngFile);
end

%% ------------------------------------------------------------------------
%% 1c: LEGEND FOR PYRUVATE ROUTING ANALYSIS
%% ------------------------------------------------------------------------
function createPyruvateRoutingLegend(outDir)
    % Creates a separate figure with legend for pyruvate routing analysis
    
    figure('Position', [100, 100, 400, 350], 'Color', 'white');
    axis off;
    
    % Create dummy lines for legend
    hold on;
    
    % Line styles
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', '-', ...
        'LineWidth', 2.5, 'DisplayName', 'PYK Flux');
    
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', '--', ...
        'LineWidth', 2.5, 'DisplayName', 'PDC Flux');
    
    plot(nan, nan, 'Color', [0, 0, 0], 'LineStyle', ':', ...
        'LineWidth', 2.5, 'DisplayName', 'PDH Flux');
    
    % Colors
    text(0.1, 0.6, 'Color Coding:', 'FontSize', 11, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Color samples
    plot(nan, nan, 'Color', [0.6294, 0.8078, 0.9804], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic SP');
    
    plot(nan, nan, 'Color', [0, 0, 0.6451], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Aerobic RP');
    
    plot(nan, nan, 'Color', [1.0000, 0.6471, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic SP');
    
    plot(nan, nan, 'Color', [0.8000, 0.5000, 0], 'LineStyle', '-', ...
        'LineWidth', 2, 'DisplayName', 'Anaerobic RP');
    
    % Feed phase
    fill([0 1 1 0], [0 0 1 1], [0.8 0.8 0.8], 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'DisplayName', 'Feed Phase');
    
    % Statistics box
    rectangle('Position', [0.1, 0.05, 0.8, 0.15], 'EdgeColor', 'black', ...
        'LineWidth', 1);
    text(0.5, 0.12, 'Feed Phase Averages Box', 'HorizontalAlignment', 'center', ...
        'FontSize', 10, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    
    legend('Location', 'best', 'FontSize', 9, 'FontName', 'Helvetica', ...
        'FontWeight', 'bold', 'Box', 'on');
    
    title('Legend: Pyruvate Routing Analysis', 'FontSize', 14, 'FontWeight', 'bold', ...
        'FontName', 'Helvetica');
    
    % Save legend
    pngFile = fullfile(outDir, 'Legend_Pyruvate_Routing.png');
    exportgraphics(gcf, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved Pyruvate Routing Legend: %s\n', pngFile);
end

%% ------------------------------------------------------------------------
%% 1d: SUPPLEMENTARY FIGURE: Model Validation
%% ------------------------------------------------------------------------
function createValidationFigures(results, fluxMap, biomassRxnID, outDir)

    fprintf('\n=== Creating Model Validation Figures ===\n');
    
    % Define colors for model and experimental
    model_color = [0, 0.4470, 0.7410];  % Blue
    exp_color = [0.8500, 0.3250, 0.0980];  % Orange-red
    
     validation_fluxes = {
        'Glucose uptake', fluxMap.qs, 'qs', true;               % flip to positive
        'Growth rate [h^{-1}]', biomassRxnID, 'mu', false;      % keep as is (always positive)
        'CO_2', fluxMap.qco2, 'qco2', false;                    % production positive, consumption negative
        'O_2 uptake', fluxMap.qo2, 'qo2', true;                 % flip to positive
        'Ethanol', fluxMap.qeth, 'qeth', false;                 % production positive, consumption negative
        'Glycerol', fluxMap.qgly, 'qgly', false;                % production positive, consumption negative
        'Acetate', fluxMap.qace, 'qace', false;                 % production positive, consumption negative
        'Pyruvate', fluxMap.qpyr, 'qpyr', false;                % production positive, consumption negative
        'Citrate', fluxMap.qcit, 'qcit', false;                 % production positive, consumption negative
        'Lactate', fluxMap.qlac, 'qlac', false;                 % production positive, consumption negative
    };
    
    n_fluxes = size(validation_fluxes, 1);
    
    % Layout: 3 rows, 4 columns
    n_rows = 3;
    n_cols = 4;
    total_panels = n_rows * n_cols; 
    
    % Conditions to plot
    cond_order = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};
    
    for c = 1:length(cond_order)
        condName = cond_order{c};
        
        if ~isfield(results, condName)
            fprintf('Warning: %s not found, skipping\n', condName);
            continue;
        end
        
        res = results.(condName);
        time_s = res.time_s;
        feed_phase = res.feed_phase;
        
        % Create figure for this condition
        fig = figure('Position', [100, 100, 1600, 1000], 'Color', 'white');
        set(0, 'DefaultAxesFontName', 'Helvetica');
        set(0, 'DefaultTextFontName', 'Helvetica');
        
        for f = 1:total_panels
            subplot(n_rows, n_cols, f);
            hold on;
            box on;
            
            if f > n_fluxes
                axis off;
                continue;
            end
            
            flux_name = validation_fluxes{f, 1};
            rxn_id = validation_fluxes{f, 2};
            exp_field = validation_fluxes{f, 3};
            flip_to_positive = validation_fluxes{f, 4};
            
            % Get model prediction (raw from FBA)
            idx = findRxnIDs(res.rxnIDs, rxn_id);
            if idx > 0
                flux_model_raw = res.fluxMatrix(:, idx);
                if flip_to_positive
                    flux_model = -flux_model_raw;   % make uptake positive
                else
                    flux_model = flux_model_raw;    % keep original sign
                end
                
                plot(time_s, flux_model, 'Color', model_color, 'LineStyle', '-', ...
                    'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4, ...
                    'MarkerFaceColor', model_color, 'DisplayName', 'Model');
            end
            
            % Get experimental data
            if isfield(res.measured, exp_field)
                exp_raw = res.measured.(exp_field);
                if flip_to_positive
                    exp_data = -exp_raw;
                else
                    exp_data = exp_raw;
                end
                plot(time_s, exp_data, 'Color', exp_color, 'LineStyle', '--', ...
                    'LineWidth', 2, 'Marker', 's', 'MarkerSize', 5, ...
                    'MarkerFaceColor', exp_color, 'DisplayName', 'Experiment');
            end
            
            % Highlight feed phase
            if any(feed_phase)
                yl = ylim;
                feed_end = max(time_s(feed_phase));
                fill([0 feed_end feed_end 0], [yl(1) yl(1) yl(2) yl(2)], ...
                    [0.85 0.85 0.85], 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');
                line([feed_end feed_end], yl, 'Color', [0.5 0.5 0.5], ...
                    'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
            end
            
            xlabel('Time [s]', 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Helvetica');
            ylabel(flux_name, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Helvetica');
            
            % Adjust y-limits with padding using visible data
            y_vals = [];
            h_lines = findobj(gca, 'Type', 'line');
            for i = 1:length(h_lines)
                ydata = get(h_lines(i), 'YData');
                y_vals = [y_vals; ydata(:)];
            end
            y_vals = y_vals(isfinite(y_vals));
            if ~isempty(y_vals)
                y_min = min(y_vals);
                y_max = max(y_vals);
                y_padding = 0.1 * (y_max - y_min);
                if y_padding < 0.01
                    y_padding = 0.1;
                end
                ylim([y_min - y_padding, y_max + y_padding]);
            end
            
            xlim([0 150]);
            grid on;
            
            if f == 1
                legend('Location', 'best', 'FontSize', 9, 'FontName', 'Helvetica');
            end
        end
        
        title_cond = strrep(condName, '_', ' ');
        sgtitle(sprintf('Flux [mmol/gDW/h] - %s', title_cond), ...
            'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Helvetica');
        
        pngFile = fullfile(outDir, sprintf('Supplementary_Model_Validation_%s.png', condName));
        exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
        fprintf('Saved: %s\n', pngFile);
        
        close(fig);
    end
    
    set(0, 'DefaultAxesFontName', 'Helvetica');
    set(0, 'DefaultTextFontName', 'Helvetica');
    fprintf('=== All validation figures saved ===\n');
end

%% ------------------------------------------------------------------------
%% 2a: SENSITIVITY ANALYSIS - GENERATE CSV DATA
%% ------------------------------------------------------------------------
function analyzePyruvateKinaseSensitivity(model0, T, fluxMap, biomassRxnID, condName, outDir)
    % Tests different PYK constraint levels over the ENTIRE FEED PHASE
    % Saves results to CSV for later plotting
    
    fprintf('\n=== PYK SENSITIVITY ANALYSIS: %s ===\n', condName);
    
    % Detect conditions
    isAnaerobic = contains(lower(condName), 'anaerobic');
    
    % Identify feed phase (high glucose uptake > 5 mmol/g/h)
    feed_phase = abs(T.qs) > 5;
    
    if ~any(feed_phase)
        fprintf('No feed phase detected (glucose uptake < 5). Skipping.\n');
        return;
    end
    
    feed_times = T.Time(feed_phase);
    feed_glc = abs(T.qs(feed_phase));
    
    fprintf('Feed phase: t = %.1f to %.1f s (%d points)\n', ...
        min(feed_times), max(feed_times), sum(feed_phase));
    fprintf('Mean glucose uptake during feed: %.2f ± %.2f mmol/g/h\n', ...
        mean(feed_glc), std(feed_glc));
    
    % Define constraint levels to test
    constraint_levels = 0:0.1:1.0;  % 0% to 100%
    n_levels = length(constraint_levels);
    
    % Pre-allocate results (will average over feed phase)
    growth_rates_avg = zeros(n_levels, 1);
    ethanol_yields_avg = zeros(n_levels, 1);
    pdc_fluxes_avg = zeros(n_levels, 1);
    pdh_fluxes_avg = zeros(n_levels, 1);
    pyk_fluxes_avg = zeros(n_levels, 1);
    pdc_pyk_ratios_avg = zeros(n_levels, 1);
    
    % Get reaction indices once
    pyk_idx = findRxnIDs(model0, fluxMap.pyk);
    pdc_idx = findRxnIDs(model0, fluxMap.pdc);
    pdh_idx = findRxnIDs(model0, fluxMap.pdh);
    eth_idx = findRxnIDs(model0, fluxMap.qeth);
    glc_idx = findRxnIDs(model0, fluxMap.qs);
    biomass_idx = findRxnIDs(model0, biomassRxnID);
    
    % Loop through constraint levels
    for i = 1:n_levels
        constraint = constraint_levels(i);
        fprintf('  Testing PYK constraint: %.0f%%... ', constraint * 100);
        
        % Initialize temporary storage for this constraint level
        growth_rates_temp = [];
        ethanol_yields_temp = [];
        pdc_fluxes_temp = [];
        pdh_fluxes_temp = [];
        pyk_fluxes_temp = [];
        pdc_pyk_ratios_temp = [];
        
        % Test at EACH feed phase time point
        feed_indices = find(feed_phase);
        for t_idx = 1:length(feed_indices)
            actual_idx = feed_indices(t_idx);
            model = model0;
            
            % Set objective
            model.c(:) = 0;
            model.c(biomass_idx) = 1;
            
            % Apply ATP maintenance
            atpm_idx = findRxnIDs(model, 'r_4046');
            if atpm_idx > 0
                model.lb(atpm_idx) = 1.5;
                model.ub(atpm_idx) = 2.5;
            end
            
            % Apply growth limits
            if isAnaerobic
                model.ub(biomass_idx) = 0.25;
            else
                model.ub(biomass_idx) = 0.45;
            end
            model.lb(biomass_idx) = 0;
            
            % Apply glucose constraint
            q_glc = T.qs(actual_idx);
            if glc_idx > 0
                model.lb(glc_idx) = q_glc;
                model.ub(glc_idx) = 0;
            end
            
            % Apply O2 constraint
            if isAnaerobic && isfield(fluxMap, 'qo2')
                o2_idx = findRxnIDs(model, fluxMap.qo2);
                if o2_idx > 0
                    model.lb(o2_idx) = 0;
                    model.ub(o2_idx) = 0;
                end
            end
            
            % Apply PYK constraint
            if pyk_idx > 0
                normal_pyk = abs(q_glc) * 2.0;
                if constraint == 0
                    model.ub(pyk_idx) = 0;
                else
                    model.ub(pyk_idx) = normal_pyk * constraint;
                end
                model.lb(pyk_idx) = 0;
            end
            
            % Solve
            sol = optimizeCbModel(model, 'max');
            
            if sol.stat == 1
                growth_rates_temp(end+1) = sol.v(biomass_idx);
                pyk_flux = sol.v(pyk_idx);
                pyk_fluxes_temp(end+1) = pyk_flux;
                
                if pdc_idx > 0
                    pdc_flux = sol.v(pdc_idx);
                    pdc_fluxes_temp(end+1) = pdc_flux;
                    if pyk_flux > 0
                        pdc_pyk_ratios_temp(end+1) = pdc_flux / pyk_flux * 100;
                    else
                        pdc_pyk_ratios_temp(end+1) = 0;
                    end
                end
                
                if pdh_idx > 0
                    pdh_fluxes_temp(end+1) = sol.v(pdh_idx);
                end
                
                if eth_idx > 0 && glc_idx > 0
                    glucose_uptake = abs(sol.v(glc_idx));
                    if glucose_uptake > 0
                        ethanol_yields_temp(end+1) = sol.v(eth_idx) / glucose_uptake;
                    else
                        ethanol_yields_temp(end+1) = 0;
                    end
                end
            end
        end
        
        % Average over feed phase time points
        if ~isempty(growth_rates_temp)
            growth_rates_avg(i) = mean(growth_rates_temp);
            ethanol_yields_avg(i) = mean(ethanol_yields_temp);
            pdc_fluxes_avg(i) = mean(pdc_fluxes_temp);
            pdh_fluxes_avg(i) = mean(pdh_fluxes_temp);
            pyk_fluxes_avg(i) = mean(pyk_fluxes_temp);
            pdc_pyk_ratios_avg(i) = mean(pdc_pyk_ratios_temp);
            
            fprintf('μ=%.3f, EtOH=%.2f, PDC/PYK=%.0f%%\n', ...
                growth_rates_avg(i), ethanol_yields_avg(i), pdc_pyk_ratios_avg(i));
        else
            fprintf('INFEASIBLE\n');
            growth_rates_avg(i) = NaN;
            ethanol_yields_avg(i) = NaN;
            pdc_fluxes_avg(i) = NaN;
            pdh_fluxes_avg(i) = NaN;
            pyk_fluxes_avg(i) = NaN;
            pdc_pyk_ratios_avg(i) = NaN;
        end
    end
    
    % Save data to CSV (no plotting)
    csvFile = fullfile(outDir, [condName '_PYK_Sensitivity_Data.csv']);
    T_sens = table();
    T_sens.PYK_Capacity = constraint_levels' * 100;
    T_sens.Growth_Rate = growth_rates_avg;
    T_sens.Ethanol_Yield = ethanol_yields_avg;
    T_sens.PDC_Flux = pdc_fluxes_avg;
    T_sens.PDH_Flux = pdh_fluxes_avg;
    T_sens.PYK_Flux = pyk_fluxes_avg;
    T_sens.PDC_PYK_Ratio = pdc_pyk_ratios_avg;
    writetable(T_sens, csvFile);
    fprintf('Saved sensitivity data: %s\n', csvFile);
end

%% ------------------------------------------------------------------------
%% 2b: SENSITIVITY ANALYSIS: EFFECT OF PYK CONSTRAINT ON METABOLIC ROUTING
%% ------------------------------------------------------------------------
function analyzeCombinedSensitivity(outDir)
    % Combines SP and RP for same aeration type into one figure (2x4)
    % Reads CSV files directly from outDir
    
    fprintf('\n=== Creating Combined Sensitivity Figures ===\n');
    
    % Define aeration types and their conditions
    aeration_config = {
        'Aerobic', {'Aerobic_S_P', 'Aerobic_R_P'}, [0.6294, 0.8078, 0.9804], [0, 0, 0.6451];
        'Anaerobic', {'Anaerobic_S_P', 'Anaerobic_R_P'}, [1.0000, 0.6471, 0], [0.8000, 0.5000, 0]
    };
    
    for a = 1:size(aeration_config, 1)
        aeration = aeration_config{a, 1};
        cond_names = aeration_config{a, 2};
        sp_color = aeration_config{a, 3};
        rp_color = aeration_config{a, 4};
        
        fprintf('  Processing %s...\n', aeration);
        
        % Load data for SP and RP
        sp_data = [];
        rp_data = [];
        
        for c = 1:2
            condName = cond_names{c};
            csvFile = fullfile(outDir, [condName '_PYK_Sensitivity_Data.csv']);
            if exist(csvFile, 'file')
                T = readtable(csvFile);
                if c == 1
                    sp_data = T;
                else
                    rp_data = T;
                end
            else
                fprintf('    Warning: %s not found\n', condName);
            end
        end
        
        if isempty(sp_data) || isempty(rp_data)
            fprintf('  Missing data for %s, skipping...\n', aeration);
            continue;
        end
        
        % Create figure (2 rows x 4 columns = 8 subplots)
        fig = figure('Position', [100, 100, 1600, 800], 'Color', 'white');
        set(0, 'DefaultAxesFontName', 'Helvetica');
        set(0, 'DefaultTextFontName', 'Helvetica');
        
        % ===== ROW 1: SP (solid line with circles) =====
        
        % Subplot 1: Growth rate
        subplot(2, 4, 1);
        valid_idx = ~isnan(sp_data.Growth_Rate);
        plot(sp_data.PYK_Capacity(valid_idx), sp_data.Growth_Rate(valid_idx), ...
            'o-', 'Color', sp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', sp_color);
        xlabel('PYK capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Growth rate [h^{-1}]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s SP: Growth rate', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        if strcmp(aeration, 'Anaerobic')
            ylim([-0.02 0.3]);
        else
            ylim([-0.02 0.6]);
        end
        
        % Subplot 2: Ethanol yield
        subplot(2, 4, 2);
        valid_eth = ~isnan(sp_data.Ethanol_Yield);
        plot(sp_data.PYK_Capacity(valid_eth), sp_data.Ethanol_Yield(valid_eth), ...
            's-', 'Color', sp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', sp_color);
        xlabel('PYK capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Ethanol yield [mol/mol]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s SP: Ethanol yield', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        ylim([0 2.2]);
        line([0 100], [2 2], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
        
        % Subplot 3: PDC/PYK ratio
        subplot(2, 4, 3);
        valid_ratio = ~isnan(sp_data.PDC_PYK_Ratio);
        plot(sp_data.PYK_Capacity(valid_ratio), sp_data.PDC_PYK_Ratio(valid_ratio), ...
            '^-', 'Color', sp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', sp_color);
        xlabel('PYK capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('PDC/PYK ratio [%]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s SP: Routing efficiency', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        ylim([0 120]);
        line([0 100], [100 100], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
        
        % Subplot 4: PDC & PDH fluxes
        subplot(2, 4, 4);
        hold on;
        valid_pdc = ~isnan(sp_data.PDC_Flux);
        if any(valid_pdc)
            plot(sp_data.PYK_Capacity(valid_pdc), sp_data.PDC_Flux(valid_pdc), ...
                's-', 'Color', [0.2, 0.6, 0.2], 'LineWidth', 2, ...
                'MarkerSize', 5, 'MarkerFaceColor', [0.2, 0.6, 0.2], 'DisplayName', 'PDC');
        end
        valid_pdh = ~isnan(sp_data.PDH_Flux);
        if any(valid_pdh)
            plot(sp_data.PYK_Capacity(valid_pdh), sp_data.PDH_Flux(valid_pdh), ...
                'd-', 'Color', [0.8, 0.3, 0.6], 'LineWidth', 2, ...
                'MarkerSize', 5, 'MarkerFaceColor', [0.8, 0.3, 0.6], 'DisplayName', 'PDH');
        end
        xlabel('PYK capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Flux [mmol/gDW/h]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s SP: PDC & PDH fluxes', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 7);
        grid on;
        xlim([0 100]);
        
        % ===== ROW 2: RP (dashed line with squares) =====
        
        % Subplot 5: Growth rate
        subplot(2, 4, 5);
        valid_idx = ~isnan(rp_data.Growth_Rate);
        plot(rp_data.PYK_Capacity(valid_idx), rp_data.Growth_Rate(valid_idx), ...
            's--', 'Color', rp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', rp_color);
        xlabel('PYK capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Growth rate [h^{-1}]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s RP: Growth rate', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        if strcmp(aeration, 'Anaerobic')
            ylim([-0.02 0.3]);
        else
            ylim([-0.02 0.6]);
        end
        
        % Subplot 6: Ethanol yield
        subplot(2, 4, 6);
        valid_eth = ~isnan(rp_data.Ethanol_Yield);
        plot(rp_data.PYK_Capacity(valid_eth), rp_data.Ethanol_Yield(valid_eth), ...
            's--', 'Color', rp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', rp_color);
        xlabel('PYK Capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Ethanol yield [mol/mol]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s RP: Ethanol Yield', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        ylim([0 2.2]);
        line([0 100], [2 2], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
        
        % Subplot 7: PDC/PYK ratio
        subplot(2, 4, 7);
        valid_ratio = ~isnan(rp_data.PDC_PYK_Ratio);
        plot(rp_data.PYK_Capacity(valid_ratio), rp_data.PDC_PYK_Ratio(valid_ratio), ...
            '^--', 'Color', rp_color, 'LineWidth', 2, ...
            'MarkerSize', 6, 'MarkerFaceColor', rp_color);
        xlabel('PYK Capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('PDC/PYK ratio [%]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s RP: Routing Efficiency', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0 100]);
        ylim([0 120]);
        line([0 100], [100 100], 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.8);
        
        % Subplot 8: PDC & PDH fluxes
        subplot(2, 4, 8);
        hold on;
        valid_pdc = ~isnan(rp_data.PDC_Flux);
        if any(valid_pdc)
            plot(rp_data.PYK_Capacity(valid_pdc), rp_data.PDC_Flux(valid_pdc), ...
                's--', 'Color', [0.2, 0.6, 0.2], 'LineWidth', 2, ...
                'MarkerSize', 5, 'MarkerFaceColor', [0.2, 0.6, 0.2], 'DisplayName', 'PDC');
        end
        valid_pdh = ~isnan(rp_data.PDH_Flux);
        if any(valid_pdh)
            plot(rp_data.PYK_Capacity(valid_pdh), rp_data.PDH_Flux(valid_pdh), ...
                'd--', 'Color', [0.8, 0.3, 0.6], 'LineWidth', 2, ...
                'MarkerSize', 5, 'MarkerFaceColor', [0.8, 0.3, 0.6], 'DisplayName', 'PDH');
        end
        xlabel('PYK Capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Flux [mmol/gDW/h]', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s RP: PDC & PDH Fluxes', aeration), 'FontSize', 11, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 7);
        grid on;
        xlim([0 100]);
        
        sgtitle(sprintf('Pyruvate Kinase Sensitivity Analysis: %s Conditions', aeration), ...
            'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Helvetica');
        
        % Save figure
        pngFile = fullfile(outDir, sprintf('PYK_Sensitivity_%s.png', aeration));
        exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
        fprintf('  Saved: %s\n', pngFile);
        close(fig);
    end
    
    fprintf('=== Combined sensitivity figures saved ===\n');
end

%% ------------------------------------------------------------------------
%% 3a: FLUX VARIABILITY ANALYSIS - SIMPLIFIED (All conditions combined)
%% ------------------------------------------------------------------------
function analyzeFluxVariabilityCombined(model0, T_data, fluxMap, biomassRxnID, outDir)
    % Performs Flux Variability Analysis for all conditions
    % Combines results into one figure (2x4 layout)
    % Saves detailed summary to text file
    
    fprintf('\n=== FLUX VARIABILITY ANALYSIS (Combined) ===\n');
    
    % Define conditions
    conditions = {'Aerobic_S_P', 'Aerobic_R_P', 'Anaerobic_S_P', 'Anaerobic_R_P'};
    aeration = {'Aerobic', 'Aerobic', 'Anaerobic', 'Anaerobic'};
    
    % Colors for conditions
    colors = {[0.6294, 0.8078, 0.9804], [0, 0, 0.6451], [1.0000, 0.6471, 0], [0.8000, 0.5000, 0]};
    
    % Results storage
    all_results = {};
    
    for cond_idx = 1:4
        condName = conditions{cond_idx};
        fprintf('\n--- Analyzing: %s ---\n', condName);
        
        % Load data for this condition
        xlsFile = fullfile(outDir, '..', [condName '.xlsx']);
        if ~exist(xlsFile, 'file')
            fprintf('  File not found: %s\n', xlsFile);
            continue;
        end
        
        T = readtable(xlsFile, 'Sheet', 'mmol');
        
        isAnaerobic = contains(lower(condName), 'anaerobic');
        
        % Use peak glucose time point
        [~, max_glc_idx] = max(abs(T.qs));
        q_glc = abs(T.qs(max_glc_idx));
        fprintf('  Time point: %.1f s, Glucose uptake: %.2f mmol/gDW/h\n', T.Time(max_glc_idx), q_glc);
        
        % Define PYK constraints to test
        if isAnaerobic
            constraint_levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
        else
            constraint_levels = [0, 0.5, 1.0];
        end
        
        % Reactions to track
        reactions = {'PYK', 'PDC', 'Ethanol', 'Biomass', 'CO2'};
        rxn_ids = {fluxMap.pyk, fluxMap.pdc, fluxMap.qeth, biomassRxnID, fluxMap.qco2};
        
        % Storage for this condition
        cond_results = struct();
        cond_results.name = condName;
        cond_results.constraints = constraint_levels;
        cond_results.growth = zeros(length(constraint_levels), 1);
        cond_results.ethanol_range = zeros(length(constraint_levels), 2);
        cond_results.pdc_range = zeros(length(constraint_levels), 2);
        cond_results.pyk_range = zeros(length(constraint_levels), 2);
        
        for c = 1:length(constraint_levels)
            constraint = constraint_levels(c);
            model = model0;
            
            % Get indices
            pyk_idx = findRxnIDs(model, fluxMap.pyk);
            glc_idx = findRxnIDs(model, fluxMap.qs);
            biomass_idx = findRxnIDs(model, biomassRxnID);
            
            % Apply constraints
            if glc_idx > 0
                model.lb(glc_idx) = T.qs(max_glc_idx);
                model.ub(glc_idx) = 0;
            end
            
            % ATP maintenance
            atpm_idx = findRxnIDs(model, 'r_4046');
            if atpm_idx > 0
                model.lb(atpm_idx) = 1.5;
                model.ub(atpm_idx) = 2.5;
            end
            
            % Growth limits
            if isAnaerobic
                model.ub(biomass_idx) = 0.25;
            else
                model.ub(biomass_idx) = 0.45;
            end
            
            % O2 constraint
            if isAnaerobic && isfield(fluxMap, 'qo2')
                o2_idx = findRxnIDs(model, fluxMap.qo2);
                if o2_idx > 0
                    model.lb(o2_idx) = 0;
                    model.ub(o2_idx) = 0;
                end
            end
            
            % PYK constraint
            if pyk_idx > 0
                normal_pyk = q_glc * 2.0;
                model.ub(pyk_idx) = normal_pyk * constraint;
                model.lb(pyk_idx) = 0;
            end
            
            % Find optimal growth
            model.c(:) = 0;
            model.c(biomass_idx) = 1;
            sol_opt = optimizeCbModel(model, 'max');
            
            if sol_opt.stat ~= 1
                cond_results.growth(c) = NaN;
                continue;
            end
            
            optimal_growth = sol_opt.f;
            cond_results.growth(c) = optimal_growth;
            
            % Skip FVA if growth too low
            if optimal_growth < 0.001
                continue;
            end
            
            % Perform FVA at 95% growth
            model.lb(biomass_idx) = 0.95 * optimal_growth;
            model.ub(biomass_idx) = optimal_growth;
            
            % Get reaction names for FVA
            rxn_names = {};
            for r = 1:length(rxn_ids)
                idx = findRxnIDs(model, rxn_ids{r});
                if idx > 0
                    rxn_names{end+1} = rxn_ids{r};
                end
            end
            
            try
                [flux_min, flux_max] = fluxVariability(model, 0, 'max', rxn_names);
                
                for r = 1:length(rxn_names)
                    if strcmp(rxn_names{r}, fluxMap.qeth)
                        cond_results.ethanol_range(c, :) = [flux_min(r), flux_max(r)];
                    elseif strcmp(rxn_names{r}, fluxMap.pdc)
                        cond_results.pdc_range(c, :) = [flux_min(r), flux_max(r)];
                    elseif strcmp(rxn_names{r}, fluxMap.pyk)
                        cond_results.pyk_range(c, :) = [flux_min(r), flux_max(r)];
                    end
                end
            catch ME
                fprintf('    FVA failed: %s\n', ME.message);
            end
        end
        
        all_results{cond_idx} = cond_results;
    end
    
    %% Create combined figure (2 rows x 4 columns)
    fig = figure('Position', [100, 100, 1400, 800], 'Color', 'white');
    set(0, 'DefaultAxesFontName', 'Helvetica');
    
    metrics = {'Growth rate', 'Ethanol yield', 'PDC flux'};
    
    for m = 1:3
        subplot(2, 3, m);
        hold on;
        
        for cond_idx = 1:4
            if isempty(all_results{cond_idx})
                continue;
            end
            res = all_results{cond_idx};
            valid = ~isnan(res.growth);
            
            if m == 1
                plot(res.constraints(valid) * 100, res.growth(valid), ...
                    'o-', 'Color', colors{cond_idx}, 'LineWidth', 2, ...
                    'MarkerSize', 6, 'MarkerFaceColor', colors{cond_idx});
            elseif m == 2
                ethanol_mean = mean(res.ethanol_range, 2);
                ethanol_mean = ethanol_mean / q_glc;
                plot(res.constraints(valid) * 100, ethanol_mean(valid), ...
                    's-', 'Color', colors{cond_idx}, 'LineWidth', 2, ...
                    'MarkerSize', 6, 'MarkerFaceColor', colors{cond_idx});
            else
                pdc_mean = mean(res.pdc_range, 2);
                plot(res.constraints(valid) * 100, pdc_mean(valid), ...
                    '^-', 'Color', colors{cond_idx}, 'LineWidth', 2, ...
                    'MarkerSize', 6, 'MarkerFaceColor', colors{cond_idx});
            end
        end
        
        xlabel('PYK Capacity [%]', 'FontSize', 11, 'FontWeight', 'bold');
        if m == 1
            ylabel('Growth rate [h^{-1}]', 'FontSize', 11, 'FontWeight', 'bold');
            ylim([0, 0.5]);
        elseif m == 2
            ylabel('Ethanol yield [mol/mol]', 'FontSize', 11, 'FontWeight', 'bold');
            ylim([0, 2.2]);
            line([0, 100], [2, 2], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
        else
            ylabel('PDC flux [mmol/gDW/h]', 'FontSize', 11, 'FontWeight', 'bold');
        end
        
        title(metrics{m}, 'FontSize', 13, 'FontWeight', 'bold');
        grid on;
        xlim([0, 100]);
        
        if m == 1
            legend(conditions, 'Location', 'east', 'FontSize', 8);
        end
    end
    
    sgtitle('Flux Variability Analysis: Metabolic Flexibility at 95% Optimal Growth', ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    pngFile = fullfile(outDir, 'FVA_Combined.png');
    exportgraphics(fig, pngFile, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('Saved FVA combined figure: %s\n', pngFile);
    
    %% Save detailed summary to text file
    saveFVASummary(all_results, conditions, outDir);
end

%% ------------------------------------------------------------------------
%% 3b: Helper: Save FVA summary to text file
%% ------------------------------------------------------------------------
function saveFVASummary(all_results, conditions, outDir)
    txtFile = fullfile(outDir, 'FVA_Summary.txt');
    fid = fopen(txtFile, 'w');
    
    fprintf(fid, '========================================\n');
    fprintf(fid, 'FLUX VARIABILITY ANALYSIS (FVA) SUMMARY\n');
    fprintf(fid, '========================================\n\n');
    fprintf(fid, 'FVA answers: For each PYK constraint level, what is the\n');
    fprintf(fid, 'range of possible fluxes while maintaining 95%% of optimal growth?\n\n');
    fprintf(fid, 'INTERPRETATION:\n');
    fprintf(fid, '  • Wide range = flexible pathway (can be up/down regulated)\n');
    fprintf(fid, '  • Narrow range = constrained pathway (essential for growth)\n');
    fprintf(fid, '  • Range = 0 = flux is fixed at that level\n\n');
    fprintf(fid, '========================================\n\n');
    
    for cond_idx = 1:4
        if isempty(all_results{cond_idx})
            continue;
        end
        
        res = all_results{cond_idx};
        condName = conditions{cond_idx};
        
        fprintf(fid, 'CONDITION: %s\n', condName);
        fprintf(fid, '----------------------------------------\n');
        fprintf(fid, 'PYK%%\tGrowth\tEthanol Min\tEthanol Max\tPDC Min\tPDC Max\n');
        
        for c = 1:length(res.constraints)
            if isnan(res.growth(c))
                continue;
            end
            fprintf(fid, '%.0f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
                res.constraints(c)*100, res.growth(c), ...
                res.ethanol_range(c,1), res.ethanol_range(c,2), ...
                res.pdc_range(c,1), res.pdc_range(c,2));
        end
        fprintf(fid, '\n');
        
        % Add interpretation
        if contains(condName, 'Anaerobic')
            fprintf(fid, 'Interpretation: Anaerobic metabolism is constrained.\n');
            fprintf(fid, '  • Ethanol production is tightly coupled to growth\n');
            fprintf(fid, '  • PYK constraint directly affects ethanol yield\n\n');
        else
            fprintf(fid, 'Interpretation: Aerobic metabolism is flexible.\n');
            fprintf(fid, '  • Wide ranges indicate multiple metabolic routes\n');
            fprintf(fid, '  • PYK constraint has minimal impact\n\n');
        end
    end
    
    fclose(fid);
    fprintf('Saved FVA summary: %s\n', txtFile);
end

%% Helper function to find reaction index
function idx = findRxnIDs(rxnList, rxnID)
    if iscell(rxnList)
        idx = find(strcmp(rxnList, rxnID));
    else
        idx = find(strcmp(rxnList.rxns, rxnID));
    end
    if isempty(idx)
        idx = 0;
    else
        idx = idx(1);
    end
end