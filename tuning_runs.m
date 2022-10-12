% Branched from mulit_case_assimilation.m 13 April 2020
% Runs the tuning cases for a given method
% Outputs a file with all the cases plus a file with the best case from prior and posterior RMSEs


% Runs cases to find best inflation and sd for each case for a given method


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
fid = fopen(strcat('tuning_out_', num2str(correction_method)), 'a+');

% Set up the optimal obs_discard_width for DO_POSSIBLE_CORRECTION
if(correction_method == DO_POSSIBLE_CORRECTION)
   obs_discard_width = 10;
else
   obs_discard_width = -1;
end

fprintf(fid, 'correction_method %i   obs_discard_width %i \n', correction_method, obs_discard_width);
fprintf(fid, 'CAREFUL: DELTA_T 0.01  MODEL_SIZE 40  FORCING 8  ENSEMBLE SIZE 80 \n');
fprintf(fid, 'CAREFUL: num_times 1100  skipping 100  obs_err_var 1.0  \n');


% Seed for the tuning case is 1 
% Need to use distinct seeds for the other 10 trials
seed = 1;

obs_intvl_vals = [5 10 15 30 60];
time_err_sd_vals = [0.00001, 0.0125, 0.025, 0.05, 0.1, 0.2];

for obs_intvl_indx = 1:5
   obs_intvl = obs_intvl_vals(obs_intvl_indx);

   for time_err_sd_indx = 1:6
      time_err_sd = time_err_sd_vals(time_err_sd_indx);
      % Don't do values that are too big
      if(time_err_sd <= obs_intvl / 100)
         fprintf(fid, 'obs_intvl %i   time_err_sd %f \n', obs_intvl, time_err_sd);

         % Explore over a range of inflations and localizations
         poss_inflation = [1 1.02 1.04 1.08 1.16 1.32 1.64];
         poss_loc = [0.125 0.15 0.175 0.2 0.25 0.4 1000000];

         for inf_indx = 1:7
            inflate = poss_inflation(inf_indx);
            for loc_indx = 1:7
               half_width = poss_loc(loc_indx);

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
   
      end 
   end
end

fclose(fid);

