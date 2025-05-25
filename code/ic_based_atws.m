% Image Contrast based Automatic Time Window Selection Algorithm
% Author: Tilal Zaheer Mukhtar
% Date: May 2025

classdef ic_based_atws
    % Class containing static methods implementing the Image Contrast based Automatic Time
    % Window Selection Algorithm
    
    methods (Static)
        function [optimal_middle_profiles, ic_array] = find_optimal_middle_profiles(hrr_profiles, cptwl, overlap_factor, max_profiles)
            % Extract the dimensions of the HRR profile matrix
            hrr_profile_count = size(hrr_profiles, 1);
            
            % Calculate the loop parameters according to the specified CPTWL and overlap
            % factor
            start = floor(cptwl/2) + 1;                           
            stop = (hrr_profile_count - 1 - floor(cptwl/2)) + mod(cptwl, 2);
            step = floor(cptwl*(1-overlap_factor));
            
            % Calculate image contrast for all middle profile values tested
            ic_array = zeros(1, size(hrr_profiles, 1));
            for middle_profile = start:step:stop
                % Form the ISAR image
                isar_image = tmc.form_isar_image(hrr_profiles, middle_profile, cptwl);
        
                % Calculate the image contrast
                ic_array(middle_profile) = ic_based_atws.get_image_contrast(isar_image);
            end
        
            % Find the local maximums for image contrast
            middle_profiles = find(ic_array);
            [~, optimal_middle_profile_indices] = findpeaks(ic_array(middle_profiles), SortStr="descend", MinPeakDistance = floor(0.5*cptwl/step), NPeaks=max_profiles);
            optimal_middle_profiles = sort(middle_profiles(optimal_middle_profile_indices));
        end
        
        function [optimal_cptwls, ic_matrix] = find_optimal_cptwls(hrr_profiles, middle_profiles, initial_cptwl, initial_step_size)
            % Extract the dimensions of the HRR profile matrix 
            middle_profile_count = size(middle_profiles, 2);

            % Loop through each middle profile
            optimal_cptwls = zeros(1, middle_profile_count);
            ic_matrix = zeros(1, middle_profile_count);
            for i = 1:middle_profile_count
                % Calculate the ISAR image contrast at the initial CPTWL and for one step
                % in the positive and negative directions
                negative_step_ic = ic_based_atws.get_image_contrast( ...
                    tmc.form_isar_image(hrr_profiles, middle_profiles(i), 2*(initial_cptwl-initial_step_size)));
                positive_step_ic = ic_based_atws.get_image_contrast( ...
                    tmc.form_isar_image(hrr_profiles, middle_profiles(i), 2*(initial_cptwl+initial_step_size)));
            
                % Determine the direction of increasing ISAR image contrast
                if positive_step_ic >= negative_step_ic
                    step_direction = 1;
                else
                    step_direction = -1;
                end
                
                % Set initial values
                intial_image_contrast = ic_based_atws.get_image_contrast( ...
                    tmc.form_isar_image(hrr_profiles, middle_profiles(i), initial_cptwl));
                cptwl = initial_cptwl;
                previous_image_contrast = intial_image_contrast;
                ic_matrix(i, cptwl) = intial_image_contrast;
                
                % Find the optimal CPTWL
                for step_size = 0:initial_step_size-1
                    while ic_matrix(i, cptwl) >= previous_image_contrast
                        previous_image_contrast = ic_matrix(i, cptwl);
                        % Update the stored image contrast and optimal CPTWL values
                        cptwl = cptwl + 2*step_direction*(initial_step_size-step_size);
                        
                        % Calculate the image contrast
                        ic_matrix(i, cptwl) = ic_based_atws.get_image_contrast( ...
                            tmc.form_isar_image(hrr_profiles, middle_profiles(i), cptwl));
                    end
                    % Revert the stored image contrast and optimal CPTWL values to
                    % previous iteration
                    cptwl = cptwl - 2*step_direction*(initial_step_size-step_size);
                end
                optimal_cptwls(i) = cptwl;
            end
        end

        function image_contrast = get_image_contrast(isar_image)
            % Calculate ISAR intensity matrix
            isar_intensity = abs(isar_image).^2;
            
            % Calculate ISAR image contrast
            image_contrast = sqrt(mean((isar_intensity - mean(isar_intensity, "all")).^2, "all"))/mean(isar_intensity, "all");
        end
    end
end
