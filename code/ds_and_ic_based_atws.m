% Doppler Spread and Image Contrast based Automatic Time Window Selection Algorithm
% Author: Tilal Zaheer Mukhtar
% Date: May 2025

classdef ds_and_ic_based_atws
    % Class containing static methods implementing the Doppler Spread and Image Contrast 
    % based Automatic Time Window Selection Algorithm
     
    methods (Static)
        function [optimal_middle_profiles, ds_array, spectrogram_middle_profiles] = find_optimal_middle_profiles(hrr_profiles, cptwl, overlap_factor, freq_resolution, max_profiles)       
            % Form a Spectrogram image from the range profiles
            [spectrogram, spectrogram_middle_profiles] = ds_and_ic_based_atws.form_spectrogram(hrr_profiles, cptwl, overlap_factor);
 
            % Calculate the Doppler Spread for each frame of the Spectrogram
            ds_array = ds_and_ic_based_atws.get_spectrogram_doppler_spread(spectrogram, freq_resolution/cptwl);
            
            % Find the local maximums for the Doppler Spread
            [~, ds_peak_indices] = findpeaks(ds_array, SortStr="descend", MinPeakDistance = floor(0.5/(1-overlap_factor)), NPeaks=max_profiles);
        
            % Convert the peak indices to middle profiles
            optimal_middle_profiles = sort(spectrogram_middle_profiles(ds_peak_indices));
        end
        
        function [optimal_cptwls, ds_matrix, ic_matrix, ds_tolerance_logical_matrix] = find_optimal_cptwls(hrr_profiles, middle_profiles, cptwl_candidates, ref_cptwl, freq_resolution, ds_tolerance_factor)
            % Number of middle profiles to find the optimal CPTWL for
            middle_profile_count = size(middle_profiles, 2);
        
            % Number of CPTWL candidates to test for each middle profile
            cptwl_candidate_count = size(cptwl_candidates, 2);
            
            % Loop through each CPTWL candidate for each middle Profile
            ds_matrix = zeros(middle_profile_count, cptwl_candidate_count);
            ic_matrix = zeros(middle_profile_count, cptwl_candidate_count);
            for i = 1:middle_profile_count
                for j = 1:cptwl_candidate_count
                    % Form the ISAR Image
                    isar_image = tmc.form_isar_image(hrr_profiles, middle_profiles(i), cptwl_candidates(j));
        
                    % Calculate the Doppler Spread
                    ds_matrix(i, j) = ds_and_ic_based_atws.get_isar_doppler_spread(isar_image, freq_resolution/cptwl_candidates(j));
        
                    % Calculate the Image Contrast
                    ic_matrix(i, j) = ds_and_ic_based_atws.get_image_contrast(isar_image);
                end
            end
            
            % Define the reference Doppler Spread matrix
            ref_cptwl_ds_matrix = repmat(ds_matrix(:, cptwl_candidates == ref_cptwl), 1, size(cptwl_candidates, 2));
        
            % Define the Doppler Spread tolerance matrix
            ds_tolerance_logical_matrix = (ds_matrix > (1-ds_tolerance_factor)*ref_cptwl_ds_matrix) ...
                & (ds_matrix < (1+ds_tolerance_factor)*ref_cptwl_ds_matrix);
            
            % Find the indices with the maximum Image Contrasts that are within the 
            % Doppler Spread tolerance
            [~, ic_max_indices] = max(ic_matrix .* ds_tolerance_logical_matrix, [], 2);
        
            % Convert the maximum indices to CPTWL values
            optimal_cptwls = cptwl_candidates(ic_max_indices);
        end

        function image_contrast = get_image_contrast(isar_image)
            % Calculate ISAR intensity matrix
            isar_intensity = abs(isar_image).^2;
            
            % Calculate ISAR image contrast
            image_contrast = sqrt(mean((isar_intensity - mean(isar_intensity, "all")).^2, "all"))/mean(isar_intensity, "all");
        end
        
        function [spectrogram, middle_profiles] = form_spectrogram(hrr_profiles, cptwl, overlap_factor)
            % Number of range profiles
            hrr_profile_count = size(hrr_profiles, 1);
            
            % Determine step size and total number of steps
            step_size = floor(cptwl*(1-overlap_factor));
            step_count = floor(hrr_profile_count/step_size);
        
            % Align all range profiles
            [range_aligned_hrr_profiles, ~, ~] = tmc.haywood_range_align(hrr_profiles, 1);
            
            % Loop through each spectrogram frame
            spectrogram = zeros(cptwl, step_count);
            middle_profiles = zeros(1, step_count);
            for i = 0:step_count-1
                middle_profiles(i+1) = i*step_size + floor(cptwl/2) + 1;
                
                % Extract subset of range profiles
                start = max(middle_profiles(i+1) - floor(cptwl/2), 1);                         
                stop = min(middle_profiles(i+1) + ceil(cptwl/2) - 1, hrr_profile_count);
                hrr_profiles_subset = range_aligned_hrr_profiles(start:stop, :);
        
                % Create window vector
                hamming_window = hamming(stop - start + 1);
        
                % Fill the spectrogram frame
                spectrogram(:, i+1) = fftshift(fft(sum(hrr_profiles_subset, 2).*hamming_window, cptwl));
            end
        end

        function doppler_spread = get_spectrogram_doppler_spread(spectrogram, freq_resolution)
            % Define function that calculates the Doppler Spread of a Spectrogram image
            k_spec = 1;
            lambda = repmat(mean(spectrogram, 1) + k_spec*std(spectrogram, 1), size(spectrogram, 1), 1);
            doppler_spread = sum((abs(spectrogram)>lambda), 1) * freq_resolution;
        end
        
        function doppler_spread = get_isar_doppler_spread(isar_image, freq_resolution)
            % Number of range profiles
            hrr_profile_count = size(isar_image, 1);
            
            % Number of Doppler Spread threshold candidates to test
            ds_threshold_candidate_count = 100;
            
            % Convert ISAR image to Doppler Profile
            isar_intensity = abs(isar_image).^2;
            isar_doppler_profile = sum(20*log10(isar_intensity), 2);
            
            % Create vector of Doppler Spread threshold candidates
            ds_threshold_candidates = linspace(min(isar_doppler_profile), max(isar_doppler_profile), ds_threshold_candidate_count);
            
            % Calculate the Doppler Spread for each threshold candidate
            doppler_profile_matrix = repmat(isar_doppler_profile, 1, ds_threshold_candidate_count);
            ds_threshold_candidate_matrix = repmat(ds_threshold_candidates, hrr_profile_count, 1);
            threshold_candidate_ds = sum(doppler_profile_matrix > ds_threshold_candidate_matrix, 1);
            
            % Fit the data points to a hyperbolic curve
            hyperbola = @(x, xData) 1 ./(x(1)*xData + x(2)) + x(3);
            x0 = [1, 1, 1];
            xData = (ds_threshold_candidates - min(isar_doppler_profile)) / (max(isar_doppler_profile) - min(isar_doppler_profile));
            yData = threshold_candidate_ds / threshold_candidate_ds(1);
            options=optimset('Display', 'off');
            x = lsqcurvefit(hyperbola, x0, xData, yData, [], [], options);
            a = x(1); b = x(2); c = x(3);

            % Determine the optimal Doppler Spread threshold
            optimal_threshold = min(isar_doppler_profile) + ((sqrt(a) - b)/a) * (max(isar_doppler_profile) - min(isar_doppler_profile));
            
            % Determine the Doppler Spread
            doppler_spread = sum(isar_doppler_profile > optimal_threshold, 1) * freq_resolution;
        end
    end
end