% Reads the output files from the tuning cases and extracts the one with the smallest posterior rmse
% Creates an output file with the optimal inflation and localization

% A list of available correction methods; this was originally also found in the subroutine programs
DO_NO_CORRECTION = 11;
DO_POSSIBLE_CORRECTION = 12;
DO_VAR_ONLY_CORRECTION = 13;
DO_IMPOSSIBLE_CORRECTION = 14;
DO_NONLIN_CORRECTION = 19;
DO_NORMAL_CORRECTION = 27;
DO_NONLINEAR_INTERPOLATION_CORRECTION = 29;

% Specify a particular algorithm to explore because this controls some parameters
correction_method = DO_POSSIBLE_CORRECTION;

% Open a file for this correction method
fid = fopen(strcat('tuning_out_', num2str(correction_method)), 'r');
%fid = fopen(strcat('tuning_out_revised', num2str(correction_method)), 'r');

fid_out = fopen(strcat('best_tuning_out_', num2str(correction_method)), 'a+');

% First line has correction method and obs_discard_width
str1 = fgetl(fid);
c_method = str2num(str1(18:22));
obs_discard_width = str2num(str1(41:end));

% Second and thire line have default info
str2 = fgetl(fid);
str3 = fgetl(fid);


obs_intvl_vals = [5 10 15 30 60];
time_err_sd_vals = [0.00001, 0.0125, 0.025, 0.05, 0.1, 0.2];

% Loop through the cases; defined by obs_intvl and time_err_sd
for obs_intvl_indx = 1:5
   obs_intvl = obs_intvl_vals(obs_intvl_indx);

   for time_err_sd_indx = 1:6
      time_err_sd = time_err_sd_vals(time_err_sd_indx);
      % Don't do values that are too big
      if(time_err_sd <= obs_intvl / 100)
         str4 = fgetl(fid);

         %  Make sure that vales are the same
         r_obs_intvl = str2num(str4(10:13));
         r_time_err_sd = str2num(str4(27:end));
         if(r_obs_intvl ~= obs_intvl) stop; end
         if(r_time_err_sd ~= time_err_sd) stop; end
 
         % Read the 49 experiment lines
         for i = 1:49
            str_obs = fgetl(fid);
            obs(i, 1:5) = str2num(str_obs);
            str_res = fgetl(fid);
            if(correction_method < 20)
               res(i, 1:8) = str2num(str_res);
            else
               res(i, 1:7) = str2num(str_res);
               res(i, 8) = -99;
            end

         end
         % Find the row with the minimum values
         [minval, min_indx] = min(res(:, 5));

         % Write out the minimum value row with the appropriate header
         fprintf(fid_out, '%f %f %f %f %f %f %f %f %d %d\n', obs_intvl, time_err_sd, res(min_indx, :));

      end
   end
end



fclose(fid);


