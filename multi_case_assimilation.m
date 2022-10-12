% Runs multiple cases of the linear assimilation algorithms


% A list of available correction methods; this was originally also found in the subroutine programs
DO_NO_CORRECTION = 11;
DO_POSSIBLE_CORRECTION = 12;
DO_VAR_ONLY_CORRECTION = 13;
DO_IMPOSSIBLE_CORRECTION = 14;
DO_NONLIN_CORRECTION = 19;
DO_NORMAL_CORRECTION = 27;
DO_INTERP_CORRECTION = 29;




% Specify a particular algorithm to explore because this controls some parameters
correction_method = DO_INTERP_CORRECTION;

% Open a file for this correction method
fid = fopen(strcat('multi_case_out_', num2str(correction_method)), 'a+');

% Standard will be to do 10 trials
num_trials = 10;

if(correction_method == DO_POSSIBLE_CORRECTION)
   obs_discard_width = 10;
else
   obs_discard_width = -1;
end

% Open the input file with optimal cases for posterior error
in_file_str = ['best_tuning_out_', num2str(correction_method)];
r = load(in_file_str);

fprintf(fid, 'correction_method %i   obs_discard_width %i \n', correction_method, obs_discard_width);
fprintf(fid, 'CAREFUL: DELTA_T 0.01  MODEL_SIZE 40  FORCING 8  ENSEMBLE SIZE 80 \n');
fprintf(fid, 'CAREFUL: num_times 1100  skipping 100  obs_err_var 1.0  \n');

% Loop through to do all the cases for this correction method
for i = 1:size(r, 1)
      obs_intvl = r(i, 1);
      time_err_sd = r(i, 2);
      inflate = r(i, 3);
      half_width = r(i, 4); 

      fprintf(fid, 'obs_intvl %i   time_err_sd %f \n', obs_intvl, time_err_sd);

      % Seed of 1 was used for the tuning runs, need num_trials independent seeds
      for seed = 2:num_trials + 1
         %[inflate, half_width, seed, correction_method, obs_intvl, time_err_sd, obs_discard_width]

         if(correction_method <= 14)
            sub_assimilate(fid, inflate, half_width, seed, correction_method, obs_intvl, time_err_sd, obs_discard_width);
         elseif(correction_method == DO_NONLIN_CORRECTION)
            sub_non_linear_assim(fid, inflate, half_width, seed, obs_intvl, time_err_sd, obs_discard_width);
         elseif(correction_method == DO_NORMAL_CORRECTION)
            num_interp = 1;
            sub_normal_time_err(fid, inflate, half_width, seed, obs_intvl, time_err_sd, ...
               obs_discard_width, num_interp);
         elseif(correction_method == DO_INTERP_CORRECTION)
            num_interp = 1;
            sub_time_err_sd_interp_assim(fid, inflate, half_width, seed, obs_intvl, time_err_sd, ...
               obs_discard_width, num_interp);
         end

      end

end

fclose(fid);

