% Translational Motion Compensation Algorithms
% Author: Tilal Zaheer Mukhtar
% Date: May 2025

classdef tmc
    % Class containing static methods implementing Haywood's Range Alignment and Yuan's 
    % Autofocus algorithms
    
    methods (Static)
        function [aligned_hrr_profiles, integer_shifts, smooth_shifts] = haywood_range_align(hrr_profiles, ref_hrr_profile_num)
            % Extract the dimensions of the HRR profile matrix
            hrr_profile_count = size(hrr_profiles, 1);
            hrr_profile_axis = 0:hrr_profile_count-1;
            range_bin_count = size(hrr_profiles, 2);
            range_bin_axis = 0:range_bin_count-1;

            % Extract the reference HRR profile
            ref_hrr_profile = hrr_profiles(ref_hrr_profile_num, :);

            % Calculate the cross correlation of each HRR profile with the reference 
            % profile
            cross_correlations = abs(ifft(fft(abs(ref_hrr_profile), [], 2) .* ...
                conj(fft(abs(hrr_profiles), [], 2)), [], 2));
        
            % Calculate the integer shifts that must be applied range alignment by finding 
            % the peak of each cross correlation
            [~, integer_shifts] = max(cross_correlations, [], 2);

            % Centre the peak of the reference HRR profile on the
            [~, ref_range_bin] = max(ref_hrr_profile);
            integer_shifts = (integer_shifts - integer_shifts(ref_hrr_profile_num)) + ...
                (floor((range_bin_count-1)/2) - ref_range_bin);

            % Unwrap the integer shifts
            integer_shifts = unwrap(integer_shifts, range_bin_count-1);
                 
            % Replace values more than 3 scaled median absolute deviations from the median
            % with a linear interpolation of neighbouring nonoutlier values
            smooth_shifts = filloutliers(integer_shifts, "linear", "median");

            % Fit a 2nd order polynomial to the shift values
            polynomial_coefficients = polyfit(hrr_profile_axis, smooth_shifts, 2);
            smooth_shifts = transpose(polyval(polynomial_coefficients, hrr_profile_axis));
                      
            % Convert the time shifts for each HRR profile to a vector of phase shifts
            phase_shifts = exp(-j*2*pi * smooth_shifts .* range_bin_axis/range_bin_count);
        
            % Align the range profiles by applying the phase shift vector to the HRR 
            % profiles
            aligned_hrr_profiles = ifft(phase_shifts .* fft(hrr_profiles, [], 2), [], 2);
        end
        
        function [phase_adjusted_hrr_profiles, phase_differences] = haywood_autofocus(hrr_profiles, ref_hrr_profile_num)           
            % Determine the mean square value of each range bin over all range profiles
            range_bin_mean_square = mean(abs(hrr_profiles).^2, 1);
        
            % Determine the variance of each range bin over all range profiles
            range_bin_variance = var(abs(hrr_profiles), 0, 1);
        
            % Select the ranges bins with mean square values greater than the average mean
            % square value as candidates for the dominant scatterer
            dominant_scatterers_candidates = range_bin_mean_square > mean(range_bin_mean_square);

            % Select the scatterer with the minimum variance over all HRR profiles as the
            % dominant scatter
            dominant_scatterer = find(range_bin_variance == min(range_bin_variance(dominant_scatterers_candidates)));
            
            % Calculate the vector of phase differences of all range profiles wit the 
            % reference phase profile at the chosen range bin 
            phase_differences = angle(hrr_profiles(:, dominant_scatterer)) - angle(hrr_profiles(ref_hrr_profile_num, dominant_scatterer));

            % Calculate the phase adjusted range profiles
            phase_adjusted_hrr_profiles = exp(-j*phase_differences) .* hrr_profiles;
        end

        function isar_image = form_isar_image(hrr_profiles, middle_profile, cptwl)
            % Extract the dimensions of the HRR profile matrix
            hrr_profile_count = size(hrr_profiles, 1); 
            
            % Extract the specified subset of HRR profiles
            start = max(middle_profile - floor(cptwl/2), 1);                         
            stop = min(middle_profile + ceil(cptwl/2) - 1, hrr_profile_count);
            hrr_profiles_subset = hrr_profiles(start:stop, :);
        
            % Create a Hamming window to apply to the HRR profiles
            hamming_window = hamming(stop - start + 1);
        
            % Create the ISAR image
            [range_aligned_profiles, ~, ~] = tmc.haywood_range_align(hrr_profiles_subset, 1);
            [phase_adjusted_profiles, ~] = tmc.haywood_autofocus(range_aligned_profiles, 1);
            isar_image = fftshift(fft((phase_adjusted_profiles.*hamming_window), [], 1), 1);
        end
    end
end
