function results = MC_Neurophotonic_Assignment_Load_And_Analyze(user_cfg)
%% GPU PMC Neurophotonic Assignment - Recompute Analysis and Make Figures
% This beginner-facing file does not launch photons.
%
% It asks students to select their completed analysis file:
%   MC_Neurophotonic_Assignment_Analysis.m
%
% Then it loads stored photon-transport summaries, recomputes the DCS/SCOS
% metrics with the selected student analysis code, and saves PNG figures in:
%   results_png

%% 1. Locate the assignment folder
% The saved data file and plotting code live beside this file.
if nargin < 1
    user_cfg = struct();
end

clc;
close all;

assignment_folder = fileparts(mfilename('fullpath'));
if isfield(user_cfg, 'data_folder')
    data_folder = user_cfg.data_folder;
else
    data_folder = fullfile(assignment_folder, 'Data_Files');
end
addpath(assignment_folder);

%% 2. Choose the student analysis file
% In normal use this opens a file picker. Automated tests can pass
% user_cfg.analysis_file to skip the picker.
if isfield(user_cfg, 'analysis_file')
    analysis_file = user_cfg.analysis_file;
else
    [analysis_name, analysis_folder] = uigetfile( ...
        {'*.m', 'MATLAB analysis files (*.m)'}, ...
        'Select your completed MC_Neurophotonic_Assignment_Analysis.m', ...
        fullfile(assignment_folder, 'MC_Neurophotonic_Assignment_Analysis.m'));

    if isequal(analysis_name, 0)
        error('No analysis file was selected.');
    end

    analysis_file = fullfile(analysis_folder, analysis_name);
end

if ~exist(analysis_file, 'file')
    error('Analysis file not found: %s', analysis_file);
end

[analysis_folder, analysis_function_name] = fileparts(analysis_file);
addpath(analysis_folder);
clear(analysis_function_name);
analysis_function = str2func(analysis_function_name);

%% 3. Choose the saved data source
% Prefer assignment-local raw transport files in Data_Files. The old solved
% results file is used only as a fallback manifest when raw files are absent.
results_file = '';

if isfield(user_cfg, 'data_file')
    results_file = user_cfg.data_file;
elseif isfield(user_cfg, 'results_file')
    results_file = user_cfg.results_file;
else
    discovered_transport_files = discover_data_folder_transport_files(data_folder);
    if ~isempty(discovered_transport_files)
        user_cfg.transport_files = discovered_transport_files;
    elseif exist(fullfile(data_folder, 'dcs_scos_study_results.mat'), 'file')
        results_file = fullfile(data_folder, 'dcs_scos_study_results.mat');
    elseif exist(fullfile(assignment_folder, 'dcs_scos_study_results.mat'), 'file')
        results_file = fullfile(assignment_folder, 'dcs_scos_study_results.mat');
    end
end

if ~isempty(results_file) && ~exist(results_file, 'file')
    error('Saved data file not found: %s', results_file);
end

%% 4. Choose the PNG output folder
% All figures are written to results_png so students know where to look.
if isfield(user_cfg, 'output_dir')
    results_png_folder = user_cfg.output_dir;
else
    results_png_folder = fullfile(assignment_folder, 'results_png');
end

if ~exist(results_png_folder, 'dir')
    mkdir(results_png_folder);
end

%% 5. Load raw transport inputs and recompute analysis
% Finished result fields in MAT files are not trusted as answers here. When
% present, they are used only to find stored raw transport files.
if isempty(results_file)
    loaded_data = struct();
else
    loaded_data = load(results_file);
end

analysis_inputs = collect_analysis_inputs(loaded_data, results_file, assignment_folder, data_folder, user_cfg);
cfg = get_plot_cfg_from_sources(loaded_data, analysis_inputs);
cfg.output_dir = results_png_folder;

fprintf('\n============================================================\n');
fprintf('GPU PMC Neurophotonic Assignment - Recompute Student Analysis\n');
fprintf('============================================================\n');
fprintf('Selected analysis file:\n%s\n', analysis_file);
if isempty(results_file)
    fprintf('Loaded raw transport files from:\n%s\n', data_folder);
else
    fprintf('Loaded data file:\n%s\n', results_file);
end
fprintf('Preferred assignment data folder:\n%s\n', data_folder);
fprintf('Saving PNG figures to:\n%s\n', results_png_folder);

results = recompute_student_results(cfg, analysis_inputs, analysis_function, analysis_file);

%% 6. Make the student PNG figures
% The plotting file creates all paper-style figures used by the assignment.
MC_Neurophotonic_Assignment_Plots(cfg, results);

%% 7. Confirm the expected figures exist
% The expected list is built from the recomputed result content.
expected_png_files = expected_pngs_from_results(results);

missing_png_files = {};
for file_index = 1:numel(expected_png_files)
    png_path = fullfile(results_png_folder, expected_png_files{file_index});
    if ~exist(png_path, 'file')
        missing_png_files{end + 1} = expected_png_files{file_index}; %#ok<AGROW>
    end
end

if ~isempty(missing_png_files)
    error('The following expected PNG files were not created: %s', strjoin(missing_png_files, ', '));
end

fprintf('\nCreated student PNG figures:\n');
for file_index = 1:numel(expected_png_files)
    fprintf('  %s\n', fullfile(results_png_folder, expected_png_files{file_index}));
end
end

function analysis_inputs = collect_analysis_inputs(loaded_data, results_file, assignment_folder, data_folder, user_cfg)
if isfield(user_cfg, 'transport_files')
    analysis_inputs = collect_inputs_from_transport_files(user_cfg.transport_files, assignment_folder, data_folder);
elseif isfield(loaded_data, 'raw_runs')
    analysis_inputs = collect_inputs_from_raw_runs(loaded_data.raw_runs);
elseif isfield(loaded_data, 'results') && isfield(loaded_data.results, 'geometry_runs')
    analysis_inputs = collect_inputs_from_result_manifest(loaded_data.results, results_file, assignment_folder, data_folder);
else
    analysis_inputs = collect_inputs_from_loaded_data(loaded_data);
end

if isempty(analysis_inputs)
    error(['No raw transport inputs were found. The analysis cannot be recomputed from finished results alone. ', ...
        'Use MAT files that contain cfg_850/transport_850, cfg_1064/transport_1064, or a raw_runs struct.']);
end
end

function transport_files = discover_data_folder_transport_files(data_folder)
transport_files = {};

if ~exist(data_folder, 'dir')
    return;
end

file_info = dir(fullfile(data_folder, '**', 'paper_gpu_transport_outputs_*.mat'));
if isempty(file_info)
    file_info = dir(fullfile(data_folder, '**', '*transport*.mat'));
end

if isempty(file_info)
    return;
end

full_paths = arrayfun(@(entry) fullfile(entry.folder, entry.name), file_info, 'UniformOutput', false);
transport_files = sort(full_paths(:)).';
end

function cfg = get_plot_cfg_from_sources(loaded_data, analysis_inputs)
if ~isempty(fieldnames(loaded_data))
    cfg = get_plot_cfg_from_loaded_data(loaded_data);
    return;
end

if isempty(analysis_inputs)
    error('No data source was found for plotting configuration.');
end

cfg = get_run_cfg_snapshot(analysis_inputs(1));
if isempty(fieldnames(cfg))
    error('Raw transport inputs did not provide a plotting configuration.');
end
end

function cfg = get_plot_cfg_from_loaded_data(loaded_data)
if isfield(loaded_data, 'cfg')
    cfg = loaded_data.cfg;
elseif isfield(loaded_data, 'cfg_850')
    cfg = loaded_data.cfg_850;
elseif isfield(loaded_data, 'cfg_1064')
    cfg = loaded_data.cfg_1064;
else
    error('The saved data file must contain cfg, cfg_850, or cfg_1064.');
end
end

function analysis_inputs = collect_inputs_from_result_manifest(saved_results, results_file, assignment_folder, data_folder)
geometry_runs = saved_results.geometry_runs;
analysis_inputs = initialize_analysis_inputs(numel(geometry_runs));

for geometry_index = 1:numel(geometry_runs)
    saved_run = geometry_runs(geometry_index);
    analysis_inputs(geometry_index).extracerebral_thickness_mm = saved_run.extracerebral_thickness_mm;
    analysis_inputs(geometry_index).cfg_snapshot = saved_run.cfg_snapshot;

    if isfield(saved_run, 'transport_file_850') && ~isempty(saved_run.transport_file_850)
        transport_file = resolve_existing_file(saved_run.transport_file_850, results_file, assignment_folder, data_folder);
        if ~isempty(transport_file)
            analysis_inputs(geometry_index).transport_file_850 = transport_file;
        end
    end

    if isfield(saved_run, 'transport_file_1064') && ~isempty(saved_run.transport_file_1064)
        transport_file = resolve_existing_file(saved_run.transport_file_1064, results_file, assignment_folder, data_folder);
        if ~isempty(transport_file)
            analysis_inputs(geometry_index).transport_file_1064 = transport_file;
        end
    end
end

analysis_inputs = remove_empty_analysis_inputs(analysis_inputs);
end

function analysis_inputs = collect_inputs_from_raw_runs(raw_runs)
analysis_inputs = initialize_analysis_inputs(numel(raw_runs));

for geometry_index = 1:numel(raw_runs)
    raw_run = raw_runs(geometry_index);
    analysis_inputs(geometry_index).extracerebral_thickness_mm = get_optional_field(raw_run, 'extracerebral_thickness_mm', []);
    analysis_inputs(geometry_index).cfg_snapshot = get_optional_field(raw_run, 'cfg_snapshot', struct());
    analysis_inputs(geometry_index) = copy_optional_field(raw_run, analysis_inputs(geometry_index), 'cfg_850');
    analysis_inputs(geometry_index) = copy_optional_field(raw_run, analysis_inputs(geometry_index), 'transport_850');
    analysis_inputs(geometry_index) = copy_optional_field(raw_run, analysis_inputs(geometry_index), 'cfg_1064');
    analysis_inputs(geometry_index) = copy_optional_field(raw_run, analysis_inputs(geometry_index), 'transport_1064');
end

analysis_inputs = remove_empty_analysis_inputs(analysis_inputs);
end

function analysis_inputs = collect_inputs_from_loaded_data(loaded_data)
analysis_inputs = initialize_analysis_inputs(1);
if isfield(loaded_data, 'cfg')
    analysis_inputs(1).cfg_snapshot = loaded_data.cfg;
end

if isfield(loaded_data, 'transport') && isfield(loaded_data, 'cfg')
    analysis_inputs(1).cfg_850 = loaded_data.cfg;
    analysis_inputs(1).transport_850 = loaded_data.transport;
end

if isfield(loaded_data, 'transport_850') && isfield(loaded_data, 'cfg_850')
    analysis_inputs(1).cfg_snapshot = loaded_data.cfg_850;
    analysis_inputs(1).cfg_850 = loaded_data.cfg_850;
    analysis_inputs(1).transport_850 = loaded_data.transport_850;
end

if isfield(loaded_data, 'transport_1064') && isfield(loaded_data, 'cfg_1064')
    analysis_inputs(1).cfg_snapshot = loaded_data.cfg_1064;
    analysis_inputs(1).cfg_1064 = loaded_data.cfg_1064;
    analysis_inputs(1).transport_1064 = loaded_data.transport_1064;
end

analysis_inputs = remove_empty_analysis_inputs(analysis_inputs);
end

function analysis_inputs = collect_inputs_from_transport_files(transport_files, assignment_folder, data_folder)
if ischar(transport_files) || isstring(transport_files)
    transport_files = cellstr(transport_files);
end

analysis_inputs = initialize_analysis_inputs(0);

for file_index = 1:numel(transport_files)
    transport_file = char(transport_files{file_index});
    if ~isfile(transport_file)
        candidate_path = fullfile(data_folder, transport_file);
        if isfile(candidate_path)
            transport_file = candidate_path;
        else
            transport_file = fullfile(assignment_folder, transport_file);
        end
    end

    if ~isfile(transport_file)
        error('Transport file not found: %s', char(transport_files{file_index}));
    end

    run_input = collect_input_from_transport_file(transport_file);
    analysis_inputs = merge_transport_file_input(analysis_inputs, run_input);
end

analysis_inputs = remove_empty_analysis_inputs(analysis_inputs);
end

function run_input = collect_input_from_transport_file(transport_file)
run_input = initialize_analysis_inputs(1);
file_variables = whos('-file', transport_file);
variable_names = {file_variables.name};

if any(strcmp(variable_names, 'cfg_850')) && any(strcmp(variable_names, 'transport_850'))
    raw_data = load(transport_file, 'cfg_850');
    run_input.cfg_snapshot = raw_data.cfg_850;
    run_input.cfg_850 = raw_data.cfg_850;
    run_input.transport_file_850 = transport_file;
end

if any(strcmp(variable_names, 'cfg_1064')) && any(strcmp(variable_names, 'transport_1064'))
    raw_data = load(transport_file, 'cfg_1064');
    run_input.cfg_snapshot = raw_data.cfg_1064;
    run_input.cfg_1064 = raw_data.cfg_1064;
    run_input.transport_file_1064 = transport_file;
end

if any(strcmp(variable_names, 'cfg')) && any(strcmp(variable_names, 'transport'))
    raw_data = load(transport_file, 'cfg');
    run_input.cfg_snapshot = raw_data.cfg;
    run_input.cfg_850 = raw_data.cfg;
    run_input.transport_file_850 = transport_file;
end
end

function analysis_inputs = initialize_analysis_inputs(number_of_runs)
empty_input = struct( ...
    'extracerebral_thickness_mm', [], ...
    'cfg_snapshot', struct(), ...
    'transport_file_850', '', ...
    'cfg_850', [], ...
    'transport_850', [], ...
    'transport_file_1064', '', ...
    'cfg_1064', [], ...
    'transport_1064', []);
analysis_inputs = repmat(empty_input, number_of_runs, 1);
end

function analysis_inputs = merge_transport_file_input(analysis_inputs, run_input)
thickness_mm = get_run_thickness(run_input);
match_index = [];

for geometry_index = 1:numel(analysis_inputs)
    if abs(double(get_run_thickness(analysis_inputs(geometry_index))) - double(thickness_mm)) < 1e-6
        match_index = geometry_index;
        break;
    end
end

if isempty(match_index)
    analysis_inputs(end + 1, 1) = run_input;
else
    analysis_inputs(match_index) = merge_single_run_input(analysis_inputs(match_index), run_input);
end
end

function merged_input = merge_single_run_input(merged_input, new_input)
field_names = fieldnames(new_input);

for field_index = 1:numel(field_names)
    field_name = field_names{field_index};
    field_value = new_input.(field_name);
    if is_nonempty_value(field_value)
        merged_input.(field_name) = field_value;
    end
end
end

function path_out = resolve_existing_file(path_in, results_file, assignment_folder, data_folder)
path_in = char(path_in);

if isfile(path_in)
    path_out = path_in;
    return;
end

[~, name, ext] = fileparts(path_in);
relative_tail = path_tail_after_folder(path_in, 'paper_gpu_outputs');
candidate_paths = { ...
    fullfile(data_folder, relative_tail), ...
    fullfile(data_folder, 'paper_gpu_outputs', relative_tail), ...
    fullfile(fileparts(results_file), [name, ext]), ...
    fullfile(data_folder, [name, ext]), ...
    fullfile(assignment_folder, [name, ext])};

path_out = '';
for candidate_index = 1:numel(candidate_paths)
    if isfile(candidate_paths{candidate_index})
        path_out = candidate_paths{candidate_index};
        return;
    end
end

if exist(data_folder, 'dir')
    recursive_matches = dir(fullfile(data_folder, '**', [name, ext]));
    if isscalar(recursive_matches)
        path_out = fullfile(recursive_matches(1).folder, recursive_matches(1).name);
    end
end
end

function relative_tail = path_tail_after_folder(path_in, folder_name)
relative_tail = '';
normalized_path = strrep(char(path_in), '/', filesep);
needle = [folder_name, filesep];
match_index = strfind(normalized_path, needle);

if ~isempty(match_index)
    relative_tail = normalized_path(match_index(end) + numel(needle):end);
end
end

function analysis_inputs = remove_empty_analysis_inputs(analysis_inputs)
keep_input = false(size(analysis_inputs));

for geometry_index = 1:numel(analysis_inputs)
    keep_input(geometry_index) = is_nonempty_value(analysis_inputs(geometry_index).transport_850) || ...
        is_nonempty_value(analysis_inputs(geometry_index).transport_1064) || ...
        is_nonempty_value(analysis_inputs(geometry_index).transport_file_850) || ...
        is_nonempty_value(analysis_inputs(geometry_index).transport_file_1064);
end

analysis_inputs = analysis_inputs(keep_input);
end

function results = recompute_student_results(cfg, analysis_inputs, analysis_function, analysis_file)
results = struct();
results.study_type = 'student_recomputed_from_raw_transport';
results.geometry_extracerebral_mm_list = zeros(numel(analysis_inputs), 1, 'single');
results.geometry_runs = repmat(struct( ...
    'extracerebral_thickness_mm', single(0), ...
    'cfg_snapshot', struct(), ...
    'transport_file_850', '', ...
    'analysis_file_850', '', ...
    'transport_file_1064', '', ...
    'analysis_file_1064', '', ...
    'results', struct()), numel(analysis_inputs), 1);

for geometry_index = 1:numel(analysis_inputs)
    input_run = analysis_inputs(geometry_index);
    geometry_results = struct();
    geometry_results.extracerebral_thickness_mm = get_run_thickness(input_run);
    geometry_results.detector_rho_mm = get_run_detector_rho(input_run);
    geometry_results.flow_change_fraction = get_run_flow_change(input_run);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Recomputing geometry %.0f mm extracerebral\n', double(geometry_results.extracerebral_thickness_mm));
    fprintf('------------------------------------------------------------\n');

    if has_case_input(input_run, '850')
        fprintf('\nStudent analysis on 850 nm transport\n');
        [case_cfg, case_transport] = load_case_input(input_run, '850');
        analysis_850 = analysis_function(case_cfg, case_transport);
        if isfield(analysis_850, 'DCS')
            geometry_results.DCS_850 = analysis_850.DCS;
        end
        if isfield(analysis_850, 'SCOS')
            geometry_results.SCOS_850 = analysis_850.SCOS;
        end
    end

    if has_case_input(input_run, '1064')
        fprintf('\nStudent analysis on 1064 nm transport\n');
        [case_cfg, case_transport] = load_case_input(input_run, '1064');
        analysis_1064 = analysis_function(case_cfg, case_transport);
        if isfield(analysis_1064, 'DCS')
            geometry_results.DCS_1064 = analysis_1064.DCS;
        end
        if isfield(analysis_1064, 'SCOS')
            geometry_results.SCOS_1064 = analysis_1064.SCOS;
        end
    end

    results.geometry_extracerebral_mm_list(geometry_index) = single(geometry_results.extracerebral_thickness_mm);
    results.geometry_runs(geometry_index).extracerebral_thickness_mm = single(geometry_results.extracerebral_thickness_mm);
    results.geometry_runs(geometry_index).cfg_snapshot = get_run_cfg_snapshot(input_run);
    results.geometry_runs(geometry_index).transport_file_850 = input_run.transport_file_850;
    results.geometry_runs(geometry_index).analysis_file_850 = analysis_file;
    results.geometry_runs(geometry_index).transport_file_1064 = input_run.transport_file_1064;
    results.geometry_runs(geometry_index).analysis_file_1064 = analysis_file;
    results.geometry_runs(geometry_index).results = geometry_results;
end

results.reference_geometry_index = get_reference_geometry_index(cfg, results.geometry_runs);
reference_run = results.geometry_runs(results.reference_geometry_index);
results.reference_extracerebral_thickness_mm = reference_run.extracerebral_thickness_mm;
results.reference_cfg_snapshot = reference_run.cfg_snapshot;
results.detector_rho_mm = reference_run.results.detector_rho_mm;
results.flow_change_fraction = reference_run.results.flow_change_fraction;

if isfield(reference_run.results, 'DCS_850')
    results.DCS = reference_run.results.DCS_850;
end

if isfield(reference_run.results, 'SCOS_850')
    results.SCOS = reference_run.results.SCOS_850;
end

if isfield(reference_run.results, 'DCS_1064')
    results.DCS_1064 = reference_run.results.DCS_1064;
end

if isfield(reference_run.results, 'SCOS_1064')
    results.SCOS_1064 = reference_run.results.SCOS_1064;
end

results.transport_file_850 = reference_run.transport_file_850;
results.analysis_file_850 = reference_run.analysis_file_850;
results.transport_file_1064 = reference_run.transport_file_1064;
results.analysis_file_1064 = reference_run.analysis_file_1064;
end

function tf = has_case_input(input_run, wavelength_label)
switch wavelength_label
    case '850'
        tf = is_nonempty_value(input_run.transport_850) || is_nonempty_value(input_run.transport_file_850);
    case '1064'
        tf = is_nonempty_value(input_run.transport_1064) || is_nonempty_value(input_run.transport_file_1064);
    otherwise
        tf = false;
end
end

function [case_cfg, case_transport] = load_case_input(input_run, wavelength_label)
switch wavelength_label
    case '850'
        if is_nonempty_value(input_run.transport_850)
            case_cfg = input_run.cfg_850;
            case_transport = input_run.transport_850;
        else
            raw_data = load(input_run.transport_file_850, 'cfg_850', 'transport_850');
            case_cfg = raw_data.cfg_850;
            case_transport = raw_data.transport_850;
        end
    case '1064'
        if is_nonempty_value(input_run.transport_1064)
            case_cfg = input_run.cfg_1064;
            case_transport = input_run.transport_1064;
        else
            raw_data = load(input_run.transport_file_1064, 'cfg_1064', 'transport_1064');
            case_cfg = raw_data.cfg_1064;
            case_transport = raw_data.transport_1064;
        end
    otherwise
        error('Unsupported wavelength label: %s', wavelength_label);
end
end

function reference_index = get_reference_geometry_index(cfg, geometry_runs)
if isfield(cfg, 'reference_extracerebral_thickness_mm')
    reference_thickness = double(cfg.reference_extracerebral_thickness_mm);
else
    reference_thickness = 15;
end

thickness_values = arrayfun(@(entry) double(entry.extracerebral_thickness_mm), geometry_runs);
[~, reference_index] = min(abs(thickness_values - reference_thickness));
end

function expected_png_files = expected_pngs_from_results(results)
expected_png_files = {'paper_fig1_simplified_geometry.png'};
reference_run = results.geometry_runs(results.reference_geometry_index);

if isfield(reference_run.results, 'DCS_850') || isfield(reference_run.results, 'DCS')
    expected_png_files{end + 1} = 'paper_fig2_simplified_dcs_fit_fraction.png';
end

if isfield(reference_run.results, 'SCOS_850') || isfield(reference_run.results, 'SCOS')
    expected_png_files{end + 1} = 'paper_fig3_simplified_scos_exposure.png';
end

if numel(results.geometry_runs) > 1
    expected_png_files{end + 1} = 'paper_fig8_simplified_geometry_comparison.png';
end
end

function value = get_optional_field(input_struct, field_name, default_value)
if isfield(input_struct, field_name)
    value = input_struct.(field_name);
else
    value = default_value;
end
end

function output_struct = copy_optional_field(input_struct, output_struct, field_name)
if isfield(input_struct, field_name)
    output_struct.(field_name) = input_struct.(field_name);
end
end

function cfg_snapshot = get_run_cfg_snapshot(input_run)
if is_nonempty_value(input_run.cfg_snapshot)
    cfg_snapshot = input_run.cfg_snapshot;
elseif is_nonempty_value(input_run.cfg_850)
    cfg_snapshot = input_run.cfg_850;
elseif is_nonempty_value(input_run.cfg_1064)
    cfg_snapshot = input_run.cfg_1064;
else
    cfg_snapshot = struct();
end
end

function thickness_mm = get_run_thickness(input_run)
if is_nonempty_value(input_run.extracerebral_thickness_mm)
    thickness_mm = input_run.extracerebral_thickness_mm;
else
    cfg_snapshot = get_run_cfg_snapshot(input_run);
    if isfield(cfg_snapshot, 'current_extracerebral_thickness_mm')
        thickness_mm = cfg_snapshot.current_extracerebral_thickness_mm;
    elseif isfield(cfg_snapshot, 'layer_thickness_mm')
        thickness_mm = single(sum(cfg_snapshot.layer_thickness_mm(1:3)));
    else
        thickness_mm = single(0);
    end
end
end

function detector_rho_mm = get_run_detector_rho(input_run)
cfg_snapshot = get_run_cfg_snapshot(input_run);

if isfield(cfg_snapshot, 'detector_rho_mm')
    detector_rho_mm = cfg_snapshot.detector_rho_mm(:);
else
    detector_rho_mm = [];
end
end

function flow_change_fraction = get_run_flow_change(input_run)
cfg_snapshot = get_run_cfg_snapshot(input_run);

if isfield(cfg_snapshot, 'BFi_baseline_cm2_s') && isfield(cfg_snapshot, 'BFi_perturbed_cm2_s') && isfield(cfg_snapshot, 'brain_layer_index')
    baseline_brain_bfi = double(cfg_snapshot.BFi_baseline_cm2_s(cfg_snapshot.brain_layer_index));
    perturbed_brain_bfi = double(cfg_snapshot.BFi_perturbed_cm2_s(cfg_snapshot.brain_layer_index));
    flow_change_fraction = (perturbed_brain_bfi - baseline_brain_bfi) / max(baseline_brain_bfi, eps);
else
    flow_change_fraction = NaN;
end
end

function tf = is_nonempty_value(value)
tf = ~(isempty(value) || (isstruct(value) && isempty(fieldnames(value))));
end
