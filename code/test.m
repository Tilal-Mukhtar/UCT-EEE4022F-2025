% Test Script
% Author: Tilal Zaheer Mukhtar
% Date: May 2025

% Reset workspace
close all;
clear all;

% Datasets
datasets = ["DAP_2010-10-06_18-00-05_002_Umoya_P872_55212_outbound" 
            "DAP_2010-10-06_18-01-54_002_Umoya_P873_50446_outbound"
            "DAP_2010-10-06_18-14-21_002_Umoya_P872_55336_inbound"
            "DAP_2010-10-09_06-50-40_008_Umoya_P873_02570_outbound"
            "DAP_2010-10-09_06-55-26_010_Umoya_P874_03468_outbound"
            "DAP_2010-10-09_06-58-05_012_Umoya_P874_03085_inbound"];

% Selected Datasets
dataset_numbers = [1, 2, 3, 4, 5, 6];

% Test parameters
initial_cptwl = 64;
initial_cptwl_step = 6;
overlap_factor = 0.8;
max_images = 30;

cptwl_candidates = [32:4:128];
ref_cptwl = 64;
ds_tolerance_factor = 0.1;

for dataset_number = dataset_numbers
    radar_data = load('..\datasets\' + datasets(dataset_number));

    hrr_profiles = radar_data.sb_HRR.G1.HRR_NoMC_calib.';
    effective_prf =  1/radar_data.sb_HRR.G1.Pattern_time;
    
    range_bin_count = size(hrr_profiles,2);
    hrr_profile_count = size(hrr_profiles,1);
    range_axis = radar_data.sb_HRR.G1.xaxis_downrange_m;
    
    % Find optimal middle profiles
    [optimal_middle_profiles, ic_array] = ic_based_atws.find_optimal_middle_profiles(hrr_profiles, initial_cptwl, overlap_factor, max_images);
    
    % Save the image contrast array
    ic_array_indices = find(ic_array > 0);

    writematrix(transpose(["Middle Profiles", ic_array_indices; ...
        "Image Contrast", ic_array(ic_array_indices); ...
        "Selected", ismember(ic_array_indices, optimal_middle_profiles)]), ...
        "..\results\ic_based_atws\dataset" + dataset_number + "\middle_profile_selection_image_contrast.csv");

    ic_array_plot = figure;
    plot(ic_array_indices, ic_array(ic_array_indices));
    hold on;
    plot(optimal_middle_profiles, ic_array(optimal_middle_profiles), LineStyle="none", Marker=".", MarkerSize=10);
    hold off;
    xlabel('HRR Profile Number');
    ylabel('Image Contrast');
    legend('','Selected Middle Profiles');
    saveas(ic_array_plot, "..\results\ic_based_atws\dataset" + dataset_number + ...
            "\middle_profile_selection_image_contrast.fig");

    % Find optimal cptwls
    [optimal_cptwls, ic_matrix] = ic_based_atws.find_optimal_cptwls(hrr_profiles, optimal_middle_profiles, initial_cptwl, initial_cptwl_step);
    
    % Save the image contrast matrix
    writematrix(ic_matrix, "..\results\ic_based_atws\dataset" + dataset_number + ...
        "\cptwl_selection_image_contrast.csv");
    
    % Save the results
    image_contrast_array = zeros(size(optimal_middle_profiles));
    for i = 1:size(optimal_middle_profiles, 2)
        % Form the ISAR image
        isar_image = tmc.form_isar_image(hrr_profiles, optimal_middle_profiles(i), optimal_cptwls(i));
        
        % Calculate the image contrast
        image_contrast = ic_based_atws.get_image_contrast(isar_image);
        image_contrast_array(i) = image_contrast;
    
        % Create Plot
        doppler_axis = (-floor(optimal_cptwls(i)/2):ceil(optimal_cptwls(i)/2)-1) * effective_prf/optimal_cptwls(i);
        isar_plot = create_isar_plot(isar_image, range_axis, doppler_axis, dataset_number, optimal_middle_profiles(i), optimal_cptwls(i), image_contrast);

        % Save the figure
        saveas(isar_plot, "..\results\ic_based_atws\dataset" + dataset_number + ...
            "\mp" + optimal_middle_profiles(i) + "cptwl" + optimal_cptwls(i) + ".fig");
    end
    
    writematrix(transpose(["Middle Profiles", optimal_middle_profiles; ...
        "CPTWLs", optimal_cptwls; ...
        "Image Contrast", image_contrast_array]), ...
        "..\results\ic_based_atws\dataset" + dataset_number + "\result.csv");
    
    ic_matrix_indices = find(ic_matrix(1, :));
    ic_array_plot = figure;
    plot(ic_matrix_indices, ic_matrix(1, ic_matrix_indices), LineStyle="none",Marker="square", MarkerSize=10);
    hold on;
    plot(optimal_cptwls(1), image_contrast_array(1), LineStyle="none", Marker=".", MarkerSize=10);
    hold off;
    xlabel('CPTWL');
    ylabel('Image Contrast');
    subtitle("Middle Profile = " + optimal_middle_profiles(1));
    legend('','Selected CPTWL');
    saveas(ic_array_plot, "..\results\ic_based_atws\dataset" + dataset_number + ...
            "\cptwl_selection_image_contrast.fig");

    % Find optimal middle profiles
    [optimal_middle_profiles, ds_array, ds_array_indices] = ds_and_ic_based_atws.find_optimal_middle_profiles(hrr_profiles, initial_cptwl, overlap_factor, effective_prf, max_images);

    writematrix(transpose(["Middle Profiles", ds_array_indices; ...
        "Doppler Spread (Hz)", ds_array; ...
        "Selected", ismember(ds_array_indices, optimal_middle_profiles)]), ...
        "..\results\ds_and_ic_based_atws\dataset" + dataset_number + "\middle_profile_selection_doppler_spread.csv");

    ds_array_plot = figure;
    plot(ds_array_indices, ds_array);
    hold on;
    plot(optimal_middle_profiles, ds_array(ismember(ds_array_indices, optimal_middle_profiles)), LineStyle="none", Marker=".", MarkerSize=10);
    hold off;
    xlabel('HRR Profile Number');
    ylabel('Doppler Spread (Hz)');
    legend('','Selected Middle Profiles');
    saveas(ds_array_plot, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
            "\middle_profile_selection_doppler_spread.fig");
    

    % Find optimal cptwls
    [optimal_cptwls, ds_matrix, ic_matrix, ds_tolerance_logical_matrix] = ds_and_ic_based_atws.find_optimal_cptwls(hrr_profiles, optimal_middle_profiles, cptwl_candidates, ref_cptwl, effective_prf, ds_tolerance_factor);

    % Save doppler spread and image contrast matrices
    writematrix(ds_matrix, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
        "\cptwl_selection_doppler_spread.csv");
    writematrix(ic_matrix, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
        "\cptwl_selection_image_contrast.csv");

    % Save the results
    image_contrast_array = zeros(size(optimal_middle_profiles));
    for i = 1:size(optimal_middle_profiles, 2)
        % Form the ISAR image
        isar_image = tmc.form_isar_image(hrr_profiles, optimal_middle_profiles(i), optimal_cptwls(i));
        
        % Calculate the image contrast
        image_contrast = ds_and_ic_based_atws.get_image_contrast(isar_image);
        image_contrast_array(i) = image_contrast;

        % Create Plot
        doppler_axis = (-floor(optimal_cptwls(i)/2):ceil(optimal_cptwls(i)/2)-1) * effective_prf/optimal_cptwls(i);
        isar_plot = create_isar_plot(isar_image, range_axis, doppler_axis, dataset_number, optimal_middle_profiles(i), optimal_cptwls(i), image_contrast);

        % Save the figure
        saveas(isar_plot, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
            "\mp" + optimal_middle_profiles(i) + "cptwl" + optimal_cptwls(i) + ".fig");
    end

    writematrix(transpose(["Middle Profiles", optimal_middle_profiles; ...
        "CPTWLs", optimal_cptwls; ...
        "Image Contrast", image_contrast_array; ...
        "Doppler Spread (Hz)", ds_array(ismember(ds_array_indices, optimal_middle_profiles))]), ...
        "..\results\ds_and_ic_based_atws\dataset" + dataset_number + "\result.csv");

    ds_array_plot = figure;
    plot(cptwl_candidates, ds_matrix(1, :));
    hold on;
    plot(cptwl_candidates(ds_tolerance_logical_matrix(1, :)), ds_matrix(1, ds_tolerance_logical_matrix(1, :)), LineStyle="none", Marker=".", MarkerSize=10);
    hold off;
    xlabel('CPTWL');
    ylabel('Doppler Spread (Hz)');
    subtitle("Middle Profile = " + optimal_middle_profiles(1));
    legend('','Within Doppler Spread Tolerance Bound');
    saveas(ds_array_plot, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
            "\cptwl_selection_doppler_spread.fig");

    ic_array_plot = figure;
    plot(cptwl_candidates, ic_matrix(1, :));
    hold on;
    plot(optimal_cptwls(1), image_contrast_array(1), LineStyle="none", Marker=".", MarkerSize=10);
    hold off;
    xlabel('CPTWL');
    ylabel('Image Contrast');
    subtitle("Middle Profile = " + optimal_middle_profiles(1));
    legend('','Selected CPTWL');
    saveas(ic_array_plot, "..\results\ds_and_ic_based_atws\dataset" + dataset_number + ...
            "\cptwl_selection_image_contrast.fig");
end

function isar_plot = create_isar_plot(isar_image, range_axis, doppler_axis, dataset, middle_profile, cptwl, image_contrast)
    % Normalise the ISAR image
    isar_image = abs(isar_image)./max(max(abs(isar_image)));

    % Convert the ISAR image to a dB scale
    isar_image = 20*log10(abs(isar_image));
    
    % Create plot
    isar_plot = imagesc(range_axis, doppler_axis, isar_image);
    xlabel('Range (m)');
    ylabel('Doppler (Hz)');
    subtitle("Dataset = " + dataset + " Middle Profile = " + middle_profile + ...
        " CPTWL = " + cptwl + " IC = " + image_contrast)
    colormap('jet');
    colorbar;
    axis xy;
    clim([-50 0]);
end
