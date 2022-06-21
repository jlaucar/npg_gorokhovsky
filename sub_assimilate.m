% Subroutinized version branched from assimilate.m on 25 March 2020
% Arguements will include the localization and inflation, a random number seed
% The Correction algorithm choice, obs_intvl, time_err_sd, and obs_discard_width

function sub_assimilate(fid, inflate, half_width, seed, correction_method, obs_intvl, time_err_sd, obs_discard_width)

% Moving on from linear model case for Elia's algorithm to dynamical estimate the time offset statistics
% given knowledge of the time error variance and assuming that the offset is the same for all observations
% at a given time but that it can vary with observation time.
% This version is designed to work with the Lorenz-96 (or 63) class of models.

% Set the random number generator to repeat with a given series
% seed = 2 was the base value for the exploratory runs
rng(seed);

% Correction method algorithm selection
DO_NO_CORRECTION = 11;
DO_POSSIBLE_CORRECTION = 12;
DO_VAR_ONLY_CORRECTION = 13;
DO_IMPOSSIBLE_CORRECTION = 14;

%fprintf(fid, 'correction_method  %i \n', correction_method);

%-----------------------------------------------------------
% Model details; The standard L96 Timestep is DELTA_T = 0.05; THESE ARE ALL FIXED FOR ALL EXPERIMENTS
global DELTA_T
DELTA_T = 0.01;
global MODEL_SIZE
MODEL_SIZE = 40;
FORCING = 8;

%fprintf(fid, 'DELTA_T %f  MODEL_SIZE %i   FORCING %f  \n', DELTA_T, MODEL_SIZE, FORCING);
%-----------------------------------------------------------

%-----------------------------------------------------------
% Time domain details; fixed for initial tuning cases
num_times = 1100;
% Skip some times at the beginning for spin-up
t_skip = 100;

%fprintf(fid, 'num_times %i  skipping %i \n', num_times, t_skip);
%-----------------------------------------------------------

%-----------------------------------------------------------
% Define the observing system characteristics
% Define the frequency of observations in terms of the delta_t
%obs_intvl = 5; Now specified as an argument

% Eventually will want correlated observation errors, but not quite yet; Chose not to explore this
obs_err_var = 1.0;
E = diag(ones(MODEL_SIZE, 1), 0) * obs_err_var;
Einv = inv(E);

% time_err_sd is now specified as an input argument
%time_err_sd = 0.0125;
time_err_var = time_err_sd^2;

%fprintf(fid, 'obs_intvl %i    obs_err_var %f    time_err_sd %f \n', obs_intvl, obs_err_var, time_err_sd);
%-----------------------------------------------------------

%-----------------------------------------------------------
% Ensemble solution details
ens_size = 80;
% SPECIFY THE WINDOW WITH for discarding nearby obs; with of 0 means just discard this variable
% A negative value means don't discard anything
%obs_discard_width = -1; This is now specified as an argument

%fprintf(fid, 'ens_size %i    obs_discard_width %i  \n', ens_size, obs_discard_width);

%-----------------------------------------------------------

%-----------------------------------------------------------
% Storage for the assimilation 
prior(1:MODEL_SIZE, 1:ens_size, 1:num_times + 1) = -99999;
post(1:MODEL_SIZE, 1:ens_size, 1:num_times + 1) = -99999;
obs_increments(1:ens_size) = -99999;
prior_mn(1:MODEL_SIZE, 1:num_times + 1) = -99999;
post_mn(1:MODEL_SIZE, 1:num_times + 1) = -99999;
prior_spread(1:MODEL_SIZE, 1:num_times+1) = -99999;
post_spread(1:MODEL_SIZE, 1:num_times+1) = -99999;
prior_err(1:MODEL_SIZE, 1:num_times+1) = -99999;
post_err(1:MODEL_SIZE, 1:num_times+1) = -99999;
mn_offset(1:num_times + 1) = -99999;
lim_mn_offset(1:num_times + 1) = -99999;
impossible_mn_offset(1:num_times + 1) = -99999;

%-----------------------------------------------------------

%-----------------------------------------------------------------------------------------

% Storage that is only changed as part of the truth run generation
fr_truth(1:MODEL_SIZE, 1:(num_times + 1) * obs_intvl + 1) = -99999;
truth(1:MODEL_SIZE, 1:num_times + 2) = -999999;
obs(1:MODEL_SIZE, 1:num_times + 1) = -99999;
true_obs(1:MODEL_SIZE, 1:num_times + 1) = -99999;
lin_true_obs(1:MODEL_SIZE, 1:num_times + 1) = -99999;
real_offset(1:num_times + 1) = -99999;
da_time(1:num_times + 1) = -99999;
actual_time(1:num_times + 1) = -99999;

% Set the initial condition to spin-up the truth initial condition
fr_truth(1:MODEL_SIZE, 1) = 0.0;
fr_truth(1, 1) = 1.0;

% Advance the fr_truth for 1000 basic time steps to let things spin up
% In subroutine version, want a different IC for each seed and don't want sequences to overlap
temp_time = 0.0;
for i = 1:seed * num_times * obs_intvl;
   [fr_truth(:, 1), temp_time] = lorenz_96_adv_1step(fr_truth(:, 1)', temp_time, FORCING);
end

% Generate the full resolution truth for all assimilation steps
for i = 1:(num_times + 1)*obs_intvl
      [fr_truth(:, i+1), temp_time] = lorenz_96_adv_1step(fr_truth(:, i)', temp_time, FORCING);
end

% Generate the actual subsampled truth at the observation times
truth = fr_truth(:, 1:obs_intvl:end);

%=====================================================================================================
% Generate a time series of observations assuming a different time error for each assimilation time
% Note that the first time is just assumed to have no time error and should be used for spin-up
real_offset(1) = 0;
violations = 0;

% Can't let things go beyond the surrouding truth times
max_offset = obs_intvl * DELTA_T;
% First observation is just 0 for all state variables
true_obs(:, 1) = 0; obs(:, 1) = 0; lin_true_obs(:, 1) = 0;
for t = 2:num_times + 1
   % Compute the time offset
   real_offset(t) = randn * time_err_sd;

   % If the offset is more than the available intvl space it is an error
   % Could change this to bound the errors by bounding the range
   if(abs(real_offset(t)) > max_offset) 
      %'Time offset bigger than max_offset'
      violations = violations + 1;
      %real_offset(t)
      real_offset(t) = sign(real_offset(t)) * max_offset / 2;
      %real_offset(t)
   end

   % Now compute the truth at this time by linear interpolation in the full resolution truth
   if(real_offset(t) < 0) 
      % What to do starting from the previous truth time
      neg_t_ratio = real_offset(t) / DELTA_T;
      % Get the ratio moving forward from the previous time which is obs_intvl times DELTA_T earlier
      t_ratio = obs_intvl + neg_t_ratio; 
      base_bot_indx = (t - 2) * obs_intvl + 1;
   else
      t_ratio = real_offset(t) / DELTA_T;
      base_bot_indx = (t - 1) * obs_intvl + 1;
   end

   full_steps = floor(t_ratio);
   part_steps = t_ratio - full_steps;
   bot_indx = base_bot_indx + full_steps;
   top_indx = bot_indx + 1;
   true_obs(:, t) = fr_truth(:, bot_indx) + part_steps * (fr_truth(:, top_indx) - fr_truth(:, bot_indx));

   % Compute a linearized approximation to the true_obs to keep track of how nonlinear things are
   [slope, stemp_time] = lorenz_96_adv_1step(truth(:, t)', -99, FORCING);
   lin_true_obs(:, t) = truth(:, t) + real_offset(t) * slope';

   % Compute the noisy observations, diagonal error variance matrix for now
   for j = 1:MODEL_SIZE
      obs(j, t) = true_obs(j, t) + randn * sqrt(obs_err_var);
   end

end

% Compute the actual time of the offset observations
for i = 1:num_times + 1
   da_time(i) = obs_intvl * DELTA_T * (i - 1) + DELTA_T;
   actual_time(i) = da_time(i) + real_offset(i);
end

% Standard deviation of difference between true_obs and truth gives a feeling for the magnitude of the errors
% Treat these matrices as vectors to get the full series
sd_temp = sqrt(mean(mean(((true_obs(1:MODEL_SIZE, t_skip:num_times) - truth(1:MODEL_SIZE, t_skip:num_times)) .^2))));

% Look at how linear things are around range of time perturbation
% Interested in the relative difference between the differences lin_diff / (true_obs - truth);
rel_lin_diff = abs(lin_true_obs(:, t_skip:num_times) - true_obs(:, t_skip:num_times)) ./ ...
   abs(truth(:, t_skip:num_times) - true_obs(:, t_skip:num_times));
% This is an unstable statistic because of cases where denominator gets small; look at median
% NOte that differences can come from the RK4 vs the raw time tendency here
sprintf('%i  %g %f %f %f \n', violations, sd_temp, median(median(rel_lin_diff)), mean(mean(rel_lin_diff)), max(max(rel_lin_diff)))
fprintf(fid, '%i  %g %f %f %f \n', violations, sd_temp, median(median(rel_lin_diff)), mean(mean(rel_lin_diff)), max(max(rel_lin_diff)));

% Plot to confirm that offset observations are working
% Turn off plotting for subroutine
if(1 == 2)
figure(1)
   plot(da_time(1:num_times), truth(37, 1:num_times), 'b*', 'markersize', 14, 'linewidth', 3);
   hold on
   plot(actual_time(1:num_times), true_obs(37, 1:num_times), 'r*', 'markersize', 14, 'linewidth', 3)
   plot(actual_time(1:num_times), lin_true_obs(37, 1:num_times), 'g*', 'markersize', 14, 'linewidth', 3)
   plot(actual_time(1:num_times), obs(37, 1:num_times), 'k*', 'markersize', 14, 'linewidth', 3)

   fr_time(1:size(fr_truth, 2)) = (1:size(fr_truth, 2)) * DELTA_T;
   plot(fr_time(1:size(fr_truth, 2)), fr_truth(37, 1:size(fr_truth, 2)), 'g-*');
end

%----------------------------------------------------------------

% Initial ensemble prior is at time 1; Shouldn't be changed between trials
for i = 1:MODEL_SIZE
   for j = 1:ens_size
      prior(i, j, 1) = truth(i, 1) + randn * 0.1;
   end
end


%=====================================================================================
% Assimilation times over parameters start here

% Explore over a range of inflations and localizations
%poss_inflation = [1 1.02 1.04 1.08 1.16 1.32 1.64];
%poss_loc = [0.125 0.15 0.175 0.2 0.25 0.4 1000000];

% inflate and half_width are now input parameters




% Do the assimilation model adavance loop
for t = 1:num_times

   % At this point can test Elia's algorithm; Start here with diagonal obs error covariance
   prior_mn(:, t) = mean(prior(:, :, t), 2);

   % Do the prior inflation here
   for j = 1:ens_size
      prior(:, j, t) = sqrt(inflate) * (prior(:, j, t) - prior_mn(:, t)) + prior_mn(:, t);
   end

   % This is an approximate prior spread; Need to clean this up
   for i = 1:MODEL_SIZE
      prior_spread(i, t) = std(prior(i, :, t));
   end

   % Generate a vector of the slope, just the time tendency hacked into the model advance; Use the ensemble mean
   [slope, stemp_time] = lorenz_96_adv_1step(prior_mn(:, t)', -99, FORCING);
   % v is Elia's notation for the slope
   v = slope';

   % Compute the estimate of the mean of the time offset distribution: Equation 16 in Elia's draft
   mn_offset_factor = v' * (Einv + Einv') / (  2 * v' * Einv *v + 2 / time_err_var);

   % SHOULD EXPLORE THE ADDITIONAL CORRECTION FACTOR THAT USES PRIOR COVARIANCE MATRIX (Elia's eqn 30)

   % Get the difference between the obs and prior ensemble mean
   d = obs(:, t) - prior_mn(:, t);   

   % Compute the estimate of the mean of the time offset distribution (final step from Elia's eqn 16)
   mn_offset(t) = mn_offset_factor * d;

   % OR REPLACE WITH IMPOSSIBLE ILLUSTRATION OF DIFFERENCE FROM TRUTH
   true_d = obs(:, t) - truth(:, t);   
   impossible_mn_offset(t) = mn_offset_factor * true_d;

   % Initial posterior guess for the assimilation is the prior
   post(:, :, t) = prior(:, :, t);

   % Compute the variance for the time offset, sigma_tau for the two correction cases
   % Impossible case has only this term, Elia draft eqn. 17
   if(correction_method == DO_IMPOSSIBLE_CORRECTION)
      var_offset = 1.0 / (v' * Einv * v + 1 / time_err_var);
   end
   if(correction_method == DO_POSSIBLE_CORRECTION)
      % Need to compute prior covariance for observations; for identity obs this is just covariance of prior ensemble
      prior_cov = cov(prior(:, :, t)');
      %var_offset = 1.0 / (v' * Einv * v + 1 / time_err_var);
      % Revised variant to replace Elia's derivation
      var_offset = 1.0 / (v' * inv(E + prior_cov) * v + 1 / time_err_var);
% WARNING: THIS IS INCONSISTENT WITH EQUATION 30 IN PAPER CHECK THAT EQN
      %%var_offset_corr = (v' * (Einv + Einv')  * prior_cov * (Einv + Einv') * v) / ...
         %%(2 * v' * Einv * v + 2 / time_err_var) ^ 2;
      %%var_offset = var_offset + var_offset_corr;
   end

   % Every state variable is observed, do them sequentially
   for i = 1:MODEL_SIZE

      % TRYING TO DEAL WITH BIAS ISSUE THAT LEADS TO DIVERGENCE
      % JUST USE THE DIFFERENCES FROM OTHER VARIABLES TO GET TIME CORRECTION FOR THIS ONE
      % Copy the vector of differences to a temporary 
      ddd = d;

      for k = -obs_discard_width:obs_discard_width
         this_indx = i + k;
         if(this_indx > MODEL_SIZE) this_indx = this_indx - MODEL_SIZE ; end
         if(this_indx < 1) this_indx = MODEL_SIZE + this_indx; end
         ddd(this_indx) = 0;
      end
      lim_mn_offset(t) = mn_offset_factor * ddd;
      
      % The observations are identity for now
      obs_prior = post(i, :, t);

      % Get an estimated time error correction forward operator test
      if(correction_method == DO_IMPOSSIBLE_CORRECTION)
         obs_prior = obs_prior + impossible_mn_offset(t) * v(i);
      elseif(correction_method == DO_POSSIBLE_CORRECTION)
         % May be a need to spinup in a less abrupt way
         %if(t < t_skip)
         obs_prior = obs_prior + lim_mn_offset(t) * v(i);
      end

      % Can include the info about the time offset error variance
      if(correction_method == DO_NO_CORRECTION) 
         corr_obs_err_var = obs_err_var;
      elseif(correction_method == DO_VAR_ONLY_CORRECTION)
         % This was an error in the original implementation, paper equations are inaccurate
         %corr_obs_err_var = obs_err_var + time_err_var * v(i);
         corr_obs_err_var = obs_err_var + time_err_var * v(i) * v(i);
      elseif(correction_method == DO_POSSIBLE_CORRECTION | correction_method == DO_IMPOSSIBLE_CORRECTION) 
         % This was an error in the original implementation, paper equations are inaccurate
         %%%corr_obs_err_var = obs_err_var + var_offset * v(i);
         corr_obs_err_var = obs_err_var + var_offset * v(i)*v(i);
      end

      [obs_increments, err] = obs_increment_eakf(obs_prior, obs(i, t), corr_obs_err_var);

      % Update the state ensembles with these increments
      for k = 1: MODEL_SIZE
         % Compute distance for localization
         loc_dist = abs(i - k) / MODEL_SIZE;
         if(loc_dist > 0.5) loc_dist = 1 - loc_dist; end
         loc = comp_cov_factor(loc_dist, half_width);

         % Only do increment computations if localization is not tiny
         if(loc > 0.0000001) 
            [state_incs, Rxy] = get_state_increments(post(k, :, t), obs_prior, obs_increments);
            post(k, :, t) = post(k, :, t) + loc * state_incs;
         end
      end
   end

   %  Compute the posetrior diagnotics
   post_mn(:, t) = mean(post(:, :, t), 2);
   for i = 1:MODEL_SIZE
      post_spread(i, t) = std(post(i, :, t));
   end

   % Advance the posterior to get the priors at the next time
   for j = 1:ens_size
      prior(:, j, t + 1) = post(:, j, t);
      for k = 1:obs_intvl
         [prior(:, j, t + 1), temp_time] = lorenz_96_adv_1step(prior(:, j, t + 1)', temp_time, FORCING);
      end
   end
end

%==================================================================================

% Compute the prior and posterior rmse over all state variables
for i = 1:MODEL_SIZE
   prior_err(i, :) = prior_mn(i, :) - truth(i, 1:num_times + 1);
   post_err(i, :) = post_mn(i, :) - truth(i, 1:num_times + 1);
end

% Compute the prior RMSE in a way consistent with the DART matlab routines
temp_rmse_sec = sqrt(sum(prior_err.^2, 1));
prior_rmse = sum(temp_rmse_sec(t_skip:num_times)) / (num_times - t_skip + 1) / sqrt(MODEL_SIZE);
temp_spread_sec = sqrt(sum(prior_spread.^2, 1));
prior_mn_spread = sum(temp_spread_sec(t_skip:num_times)) / (num_times - t_skip + 1) / sqrt(MODEL_SIZE);

% TEST OF DART consistent computation of spread

temp_rmse_sec = sqrt(sum(post_err.^2, 1));
post_rmse = sum(temp_rmse_sec(t_skip:num_times)) / (num_times - t_skip + 1) / sqrt(MODEL_SIZE);
temp_spread_sec = sqrt(sum(post_spread.^2, 1));
post_mn_spread = sum(temp_spread_sec(t_skip:num_times)) / (num_times - t_skip + 1) / sqrt(MODEL_SIZE);


% Turn off graphics for long runs
if(1 == 2)
   % Do the plots of errors for a few ensemble members
   for i = 1:min(4, MODEL_SIZE)
      figure(2)
      subplot(2, 2, i)
      % Plot the prior ensemble members
      for j = 1:ens_size
         plot(squeeze(prior(i, j, 1:num_times)) - truth(i, 1:num_times)', 'g')
         hold on
      end
   
      % Plot the differences of the prior and posterior from the truth for subset of variables
      prior_mn_err = prior_mn(i, 1:num_times) - truth(i, 1:num_times);
      plot(prior_mn_err(1:num_times), 'k', 'linewidth', 3);
      post_mn_err = post_mn(i, 1:num_times) - truth(i, 1:num_times);
      plot(post_mn_err(1:num_times), 'b', 'linewidth', 3);
   
      % Look at a series of overlaid mean errors
      figure(4)
      plot(prior_mn_err(1:num_times), 'k')
      hold on
   
      figure(5)
      subplot(2, 2, i)
      for j = 1:ens_size
         plot(squeeze(prior(i, j, 1:num_times)), 'g')
         hold on
      end
      plot(prior_mn(i, 1:num_times), 'b', 'linewidth', 2);
      plot(truth(i, 1:num_times), 'k', 'linewidth', 2);
   end
end

% Compute the rmse of the possible and impossible estimates of the time offset
rmse_offset = sqrt(mean((mn_offset(t_skip:num_times) - real_offset(t_skip:num_times)) .^ 2));
rmse_impossible_offset = sqrt(mean((impossible_mn_offset(t_skip:num_times) - real_offset(t_skip:num_times)) .^ 2));

%ccc = corrcoef(mn_offset(t_skip:num_times), real_offset(t_skip:num_times));
%correlation = ccc(1, 2)

%ccc = corrcoef(impossible_mn_offset(t_skip:num_times), real_offset(t_skip:num_times));
%impossible_corr = ccc(1, 2)

% Print as a single string
sprintf('%f %f %f %f %f %f %f %f \n', inflate, half_width, prior_rmse, prior_mn_spread, post_rmse, post_mn_spread, rmse_offset, rmse_impossible_offset)
fprintf(fid, '%f %f %f %f %f %f %f %f \n', inflate, half_width, prior_rmse, prior_mn_spread, post_rmse, post_mn_spread, rmse_offset, rmse_impossible_offset);


