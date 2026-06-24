function MC_Neurophotonic_Assignment_Plots(cfg, results)
%% GPU PMC Neurophotonic Assignment - Plot Version
% This file saves the paper-style summary figures for the assignment.
%
% Teaching organization:
%   1. Prepare the data that will be plotted.
%   2. Draw the simplified tissue geometry.
%   3. Draw the DCS reference figure.
%   4. Draw the SCOS reference figure.
%   5. Draw the multi-geometry comparison figure.
%
% The plotting workflow is kept in one top-down function. Small helpers that
% only selected fields, styled axes, or skipped empty curves have been inlined
% so first-time users can follow the figure logic without jumping around.

%% 1. Prepare plotting data
% Purpose:
%   Convert the result tree into one geometry array and choose the reference
%   geometry used for the DCS and SCOS detail figures.
%
% Variables to notice:
%   outdir          - folder where figure files are saved
%   geometry_runs   - one entry per extracerebral thickness
%   reference_run   - geometry used for the reference panels

outdir = cfg.output_dir;

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

if isfield(results, 'geometry_runs')
    geometry_runs = results.geometry_runs;
else
    geometry_runs = [];
    geometry_runs(1).extracerebral_thickness_mm = single(sum(cfg.layer_thickness_mm(1:3)));
    geometry_runs(1).cfg_snapshot = cfg;
    geometry_runs(1).transport_file = '';
    geometry_runs(1).analysis_file = '';
    geometry_runs(1).results = results;
end

if isfield(results, 'reference_geometry_index') && ~isempty(results.reference_geometry_index)
    reference_index = results.reference_geometry_index;
else
    thickness_values = arrayfun(@(entry) double(entry.extracerebral_thickness_mm), geometry_runs);
    [~, reference_index] = min(abs(thickness_values - 15));
end

reference_run = geometry_runs(reference_index);

% Colors are grouped near the top so students can quickly customize figures.
layer_colors = [ ...
    0.86 0.66 0.58; ... % Scalp
    0.78 0.78 0.84; ... % Skull
    0.78 0.90 0.97; ... % CSF
    0.93 0.76 0.66];    % Brain

sd_colors = [ ...
    0.16 0.34 0.76; ...
    0.21 0.60 0.90; ...
    0.24 0.82 0.74; ...
    0.43 0.87 0.53; ...
    0.78 0.88 0.24; ...
    0.98 0.74 0.20; ...
    0.92 0.43 0.22; ...
    0.68 0.24 0.18];

dcs850_color = [0.54 0.29 0.56];
dcs1064_color = [0.39 0.79 0.43];
scos850_color = [0.92 0.32 0.22];

%% 2. Draw the simplified layered geometry
% Purpose:
%   Show the tissue layers and detector positions for each simulated
%   extracerebral thickness.
%
% Variables to notice:
%   layer_thickness_mm  - thickness of scalp, skull, CSF, and brain
%   z_edges             - cumulative layer boundaries in depth

f = figure('Color', 'w', 'Position', [100 100 1400 420]);
tiledlayout(1, numel(geometry_runs), 'Padding', 'compact', 'TileSpacing', 'compact');

for geometry_index = 1:numel(geometry_runs)
    nexttile;
    hold on;
    axis equal;
    axis off;

    layer_thickness_mm = double(geometry_runs(geometry_index).cfg_snapshot.layer_thickness_mm(:));
    z_edges = [0; cumsum(layer_thickness_mm)];
    x_limits = [-5, 45];

    for layer_index = 1:numel(layer_thickness_mm)
        patch([x_limits(1), x_limits(2), x_limits(2), x_limits(1)], ...
            [z_edges(layer_index), z_edges(layer_index), z_edges(layer_index + 1), z_edges(layer_index + 1)], ...
            layer_colors(layer_index, :), 'EdgeColor', [0.35 0.35 0.35], 'LineWidth', 0.8);

        text(mean(x_limits), 0.5 * (z_edges(layer_index) + z_edges(layer_index + 1)), ...
            sprintf('%s  %.1f mm', geometry_runs(geometry_index).cfg_snapshot.layer_names{layer_index}, layer_thickness_mm(layer_index)), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.15 0.15 0.15]);
    end

    plot(0, 0, 'ko', 'MarkerFaceColor', [0.95 0.78 0.18], 'MarkerSize', 8);
    text(0, -1.4, 'Source', 'HorizontalAlignment', 'center', 'FontSize', 10);

    detector_positions = double(geometry_runs(geometry_index).cfg_snapshot.detector_rho_mm(:));
    for detector_index = 1:numel(detector_positions)
        plot(detector_positions(detector_index), 0, 'ks', 'MarkerFaceColor', [0.20 0.20 0.20], 'MarkerSize', 5);
    end

    plot([0, detector_positions(end)], [-0.4, -0.4], 'k-', 'LineWidth', 1.0);
    text(0.5 * detector_positions(end), -2.0, sprintf('Detectors 5 to %d mm', round(detector_positions(end))), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
    set(gca, 'YDir', 'reverse');
    xlim(x_limits);
    ylim([-3, z_edges(end)]);
    title(sprintf('Simplified geometry: %.0f mm extracerebral', double(geometry_runs(geometry_index).extracerebral_thickness_mm)));
end

saveas(f, fullfile(outdir, 'paper_fig1_simplified_geometry.png'));
close(f);

%% 3. Draw the DCS fit-fraction reference figure
% Purpose:
%   Show how DCS sensitivity, CoV, and CNR change as more of the correlation
%   decay curve is included in the BFI fit.
%
% Variables to notice:
%   D               - selected DCS block, usually the 850 nm result
%   fit_percentage  - x-axis showing how much of the decay is fitted

has_dcs_reference = isfield(reference_run.results, 'DCS_850') || isfield(reference_run.results, 'DCS');

if has_dcs_reference
    if isfield(reference_run.results, 'DCS_850')
        D = reference_run.results.DCS_850;
    else
        D = reference_run.results.DCS;
    end

    rho = reference_run.results.detector_rho_mm(:);
    fit_percentage = 100 .* D.fit_fractions(:);
    panel_labels = {'(a)', '(b)', '(c)'};
    metric_names = {'sensitivity', 'cov', 'cnr'};
    y_labels = {'Cerebral sensitivity (%)', 'CoV of BFI at 10 Hz', 'Contrast-to-noise ratio of BFI at 10 Hz'};

    f = figure('Color', 'w', 'Position', [100 100 1320 360]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for metric_index = 1:numel(metric_names)
        nexttile;
        hold on;
        data = D.(metric_names{metric_index});

        for detector_index = 1:numel(rho)
            if strcmp(metric_names{metric_index}, 'sensitivity')
                y_values = 100 .* data(detector_index, :);
            else
                y_values = data(detector_index, :);
            end

            plot(fit_percentage, y_values, '-', 'Color', sd_colors(detector_index, :), ...
                'LineWidth', 2.0, 'Marker', 'o', 'MarkerSize', 4.5, ...
                'MarkerFaceColor', sd_colors(detector_index, :), ...
                'DisplayName', sprintf('%d mm', round(rho(detector_index))));
        end

        set(gca, 'XLim', [min(fit_percentage), max(fit_percentage)], 'XTick', fit_percentage(:).');
        xlabel('Percentage of decay fit (%)');
        ylabel(y_labels{metric_index});
        set(gca, 'Box', 'on', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.5);
        text(gca, 0.03, 0.95, panel_labels{metric_index}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 11);

        if strcmp(metric_names{metric_index}, 'cov') || strcmp(metric_names{metric_index}, 'cnr')
            set(gca, 'YScale', 'log');
        end
    end

    lgd = legend('Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', numel(rho), 'Box', 'off');
    lgd.Layout.Tile = 'north';
    lgd.FontSize = 9;
    saveas(f, fullfile(outdir, 'paper_fig2_simplified_dcs_fit_fraction.png'));
    close(f);
end

%% 4. Draw the SCOS exposure reference figure
% Purpose:
%   Show how SCOS sensitivity, CoV, and CNR change with exposure time.
%
% Variables to notice:
%   S     - selected SCOS block, usually the 850 nm result
%   texp  - camera exposure times

has_scos_reference = isfield(reference_run.results, 'SCOS_850') || isfield(reference_run.results, 'SCOS');

if has_scos_reference
    if isfield(reference_run.results, 'SCOS_850')
        S = reference_run.results.SCOS_850;
    else
        S = reference_run.results.SCOS;
    end

    texp = S.exposure_s(:);
    rho = reference_run.results.detector_rho_mm(:);
    panel_labels = {'(a)', '(b)', '(c)'};
    metric_names = {'sensitivity', 'cov', 'cnr'};
    y_labels = {'Cerebral sensitivity (%)', 'CoV of BFI at 10 Hz', 'Contrast-to-noise ratio of BFI at 10 Hz'};

    f = figure('Color', 'w', 'Position', [100 100 1320 360]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for metric_index = 1:numel(metric_names)
        nexttile;
        hold on;
        data = S.(metric_names{metric_index}).';

        for detector_index = 1:numel(rho)
            if strcmp(metric_names{metric_index}, 'sensitivity')
                y_values = 100 .* data(:, detector_index);
            else
                y_values = data(:, detector_index);
            end

            plot(texp, y_values, '-', 'Color', sd_colors(detector_index, :), ...
                'LineWidth', 2.0, 'DisplayName', sprintf('%d mm', round(rho(detector_index))));
        end

        set(gca, 'XScale', 'log');
        if strcmp(metric_names{metric_index}, 'cov') || strcmp(metric_names{metric_index}, 'cnr')
            set(gca, 'YScale', 'log');
        end

        xlabel('Exposure time (s)');
        ylabel(y_labels{metric_index});
        set(gca, 'Box', 'on', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.5);
        text(gca, 0.03, 0.95, panel_labels{metric_index}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 11);
    end

    lgd = legend('Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', numel(rho), 'Box', 'off');
    lgd.Layout.Tile = 'north';
    lgd.FontSize = 9;
    saveas(f, fullfile(outdir, 'paper_fig3_simplified_scos_exposure.png'));
    close(f);
end

%% 5. Draw the multi-geometry comparison figure
% Purpose:
%   Compare 850 nm DCS, 1064 nm DCS, and 850 nm SCOS across simplified
%   extracerebral thicknesses.
%
% Variables to notice:
%   dcs_850_*   - full-fit DCS metrics at 850 nm
%   dcs_1064_*  - full-fit DCS metrics at 1064 nm
%   scos_850_*  - best-exposure SCOS metrics at 850 nm

has_geometry_comparison = false;
for geometry_index = 1:numel(geometry_runs)
    run_results = geometry_runs(geometry_index).results;
    if isfield(run_results, 'DCS_850') || isfield(run_results, 'DCS') || isfield(run_results, 'DCS_1064') || isfield(run_results, 'SCOS_850') || isfield(run_results, 'SCOS')
        has_geometry_comparison = true;
        break;
    end
end

if numel(geometry_runs) > 1 && has_geometry_comparison
    all_sensitivity_values = [];
    all_cov_values = [];
    all_cnr_values = [];

    for geometry_index = 1:numel(geometry_runs)
        run_results = geometry_runs(geometry_index).results;

        dcs_850_sensitivity = [];
        dcs_850_cov = [];
        dcs_850_cnr = [];
        if isfield(run_results, 'DCS_850')
            D = run_results.DCS_850;
        elseif isfield(run_results, 'DCS')
            D = run_results.DCS;
        else
            D = [];
        end
        if ~isempty(D)
            full_fraction_index = find(abs(D.fit_fractions - 1) < 1e-9, 1);
            if isempty(full_fraction_index)
                full_fraction_index = numel(D.fit_fractions);
            end
            dcs_850_sensitivity = D.sensitivity(:, full_fraction_index);
            dcs_850_cov = D.cov(:, full_fraction_index);
            dcs_850_cnr = D.cnr(:, full_fraction_index);
        end

        dcs_1064_sensitivity = [];
        dcs_1064_cov = [];
        dcs_1064_cnr = [];
        if isfield(run_results, 'DCS_1064')
            D = run_results.DCS_1064;
            full_fraction_index = find(abs(D.fit_fractions - 1) < 1e-9, 1);
            if isempty(full_fraction_index)
                full_fraction_index = numel(D.fit_fractions);
            end
            dcs_1064_sensitivity = D.sensitivity(:, full_fraction_index);
            dcs_1064_cov = D.cov(:, full_fraction_index);
            dcs_1064_cnr = D.cnr(:, full_fraction_index);
        end

        scos_850_sensitivity = [];
        scos_850_cov = [];
        scos_850_cnr = [];
        if isfield(run_results, 'SCOS_850')
            S = run_results.SCOS_850;
        elseif isfield(run_results, 'SCOS')
            S = run_results.SCOS;
        else
            S = [];
        end
        if ~isempty(S)
            scos_850_sensitivity = S.best_sensitivity(:);
            scos_850_cov = S.best_cov(:);
            scos_850_cnr = S.best_cnr(:);
        end

        all_sensitivity_values = [all_sensitivity_values; dcs_850_sensitivity(:); dcs_1064_sensitivity(:); scos_850_sensitivity(:)]; %#ok<AGROW>
        all_cov_values = [all_cov_values; dcs_850_cov(:); dcs_1064_cov(:); scos_850_cov(:)]; %#ok<AGROW>
        all_cnr_values = [all_cnr_values; dcs_850_cnr(:); dcs_1064_cnr(:); scos_850_cnr(:)]; %#ok<AGROW>
    end

    all_sensitivity_values = all_sensitivity_values(isfinite(all_sensitivity_values));
    all_cov_values = all_cov_values(isfinite(all_cov_values) & all_cov_values > 0);
    all_cnr_values = all_cnr_values(isfinite(all_cnr_values) & all_cnr_values > 0);

    sensitivity_limits = [0, 1.05 * max([100 .* all_sensitivity_values; 1])];
    if isempty(all_cov_values)
        cov_limits = [1e-6, 1];
    else
        cov_limits = [10 ^ floor(log10(max(min(all_cov_values), 1e-6))), 10 ^ ceil(log10(max(all_cov_values)))];
    end
    if isempty(all_cnr_values)
        cnr_limits = [1e-6, 1];
    else
        cnr_limits = [10 ^ floor(log10(max(min(all_cnr_values), 1e-6))), 10 ^ ceil(log10(max(all_cnr_values)))];
    end

    panel_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)'};
    f = figure('Color', 'w', 'Position', [100 100 980 950]);
    tiledlayout(numel(geometry_runs), 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    for geometry_index = 1:numel(geometry_runs)
        run_entry = geometry_runs(geometry_index);
        run_results = run_entry.results;
        rho = run_results.detector_rho_mm(:);
        geometry_label = sprintf('Geometry: %.0f mm extracerebral', double(run_entry.extracerebral_thickness_mm));

        dcs_850_sensitivity = [];
        dcs_850_cov = [];
        dcs_850_cnr = [];
        if isfield(run_results, 'DCS_850')
            D = run_results.DCS_850;
        elseif isfield(run_results, 'DCS')
            D = run_results.DCS;
        else
            D = [];
        end
        if ~isempty(D)
            full_fraction_index = find(abs(D.fit_fractions - 1) < 1e-9, 1);
            if isempty(full_fraction_index)
                full_fraction_index = numel(D.fit_fractions);
            end
            dcs_850_sensitivity = D.sensitivity(:, full_fraction_index);
            dcs_850_cov = D.cov(:, full_fraction_index);
            dcs_850_cnr = D.cnr(:, full_fraction_index);
        end

        dcs_1064_sensitivity = [];
        dcs_1064_cov = [];
        dcs_1064_cnr = [];
        if isfield(run_results, 'DCS_1064')
            D = run_results.DCS_1064;
            full_fraction_index = find(abs(D.fit_fractions - 1) < 1e-9, 1);
            if isempty(full_fraction_index)
                full_fraction_index = numel(D.fit_fractions);
            end
            dcs_1064_sensitivity = D.sensitivity(:, full_fraction_index);
            dcs_1064_cov = D.cov(:, full_fraction_index);
            dcs_1064_cnr = D.cnr(:, full_fraction_index);
        end

        scos_850_sensitivity = [];
        scos_850_cov = [];
        scos_850_cnr = [];
        if isfield(run_results, 'SCOS_850')
            S = run_results.SCOS_850;
        elseif isfield(run_results, 'SCOS')
            S = run_results.SCOS;
        else
            S = [];
        end
        if ~isempty(S)
            scos_850_sensitivity = S.best_sensitivity(:);
            scos_850_cov = S.best_cov(:);
            scos_850_cnr = S.best_cnr(:);
        end

        nexttile;
        hold on;
        if ~isempty(dcs_850_sensitivity)
            plot(rho, 100 .* dcs_850_sensitivity, '-o', 'Color', dcs850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs850_color, 'DisplayName', '850 nm DCS');
        end
        if ~isempty(dcs_1064_sensitivity)
            plot(rho, 100 .* dcs_1064_sensitivity, '-o', 'Color', dcs1064_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs1064_color, 'DisplayName', '1064 nm DCS');
        end
        if ~isempty(scos_850_sensitivity)
            plot(rho, 100 .* scos_850_sensitivity, '-o', 'Color', scos850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', scos850_color, 'DisplayName', '850 nm SCOS');
        end
        ylabel('Cerebral sensitivity (%)');
        xlabel('SD separation (mm)');
        ylim(sensitivity_limits);
        set(gca, 'Box', 'on', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.5);
        text(gca, 0.03, 0.95, panel_labels{3 * (geometry_index - 1) + 1}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 11);
        title(geometry_label, 'FontWeight', 'bold', 'FontSize', 10);

        nexttile;
        hold on;
        if ~isempty(dcs_850_cov)
            plot(rho, dcs_850_cov, '-o', 'Color', dcs850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs850_color, 'DisplayName', '850 nm DCS');
        end
        if ~isempty(dcs_1064_cov)
            plot(rho, dcs_1064_cov, '-o', 'Color', dcs1064_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs1064_color, 'DisplayName', '1064 nm DCS');
        end
        if ~isempty(scos_850_cov)
            plot(rho, scos_850_cov, '-o', 'Color', scos850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', scos850_color, 'DisplayName', '850 nm SCOS');
        end
        ylabel('CoV of BFI at 10 Hz');
        xlabel('SD separation (mm)');
        set(gca, 'YScale', 'log');
        ylim(cov_limits);
        set(gca, 'Box', 'on', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.5);
        text(gca, 0.03, 0.95, panel_labels{3 * (geometry_index - 1) + 2}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 11);

        nexttile;
        hold on;
        if ~isempty(dcs_850_cnr)
            plot(rho, dcs_850_cnr, '-o', 'Color', dcs850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs850_color, 'DisplayName', '850 nm DCS');
        end
        if ~isempty(dcs_1064_cnr)
            plot(rho, dcs_1064_cnr, '-o', 'Color', dcs1064_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', dcs1064_color, 'DisplayName', '1064 nm DCS');
        end
        if ~isempty(scos_850_cnr)
            plot(rho, scos_850_cnr, '-o', 'Color', scos850_color, 'LineWidth', 2.0, 'MarkerSize', 5, 'MarkerFaceColor', scos850_color, 'DisplayName', '850 nm SCOS');
        end
        ylabel('Contrast-to-noise ratio of BFI at 10 Hz');
        xlabel('SD separation (mm)');
        set(gca, 'YScale', 'log');
        ylim(cnr_limits);
        set(gca, 'Box', 'on', 'LineWidth', 0.8, 'FontName', 'Arial', 'FontSize', 10, ...
            'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.88 0.88 0.88], 'GridAlpha', 0.5);
        text(gca, 0.03, 0.95, panel_labels{3 * (geometry_index - 1) + 3}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontWeight', 'bold', 'FontSize', 11);
    end

    lgd = legend('Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4, 'Box', 'off');
    lgd.Layout.Tile = 'north';
    lgd.FontSize = 9;
    lgd.Title.String = 'Curves: modality/wavelength. Rows: geometry.';
    lgd.Title.FontWeight = 'normal';
    saveas(f, fullfile(outdir, 'paper_fig8_simplified_geometry_comparison.png'));
    close(f);
end
end
