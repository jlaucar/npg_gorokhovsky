% Read in the trials files
% Partition data into nice arrays
% Plots of RMSE and RMSE of estimates of the time offset

num_trials = 10;
lines_per_case = 2*num_trials + 1;

% Create an array with the following dimensions ddd(num_cases, 12)
% THe 12 values are (obs_intvl, time_err_sd, method, trial, infl, loc, prior_rmse, ...
%   prior_spread, post_rmse, post_spread, offset_rmse, perfect_offset_rmes)

% Loop through the 5 different methods
rmse_line = 0;
mean_line = 0;
for method_indx = 1:5
   method_vals = [11 12 13 14 29];
   method = method_vals(method_indx);

   % Open the file of the multicase output for this method
   fid = fopen(['multi_case_out_', num2str(method)], 'r');

   % First line has correction method and obs_discard_width
   str1 = fgetl(fid);
   c_method = str2num(str1(18:22));
   obs_discard_width = str2num(str1(41:end));

   % Second and thire line have default info
   str2 = fgetl(fid);
   str3 = fgetl(fid);

   % Loop over the number of cases for this method type
   obs_intvl_vals = [5 10 15 30 60];
   time_err_sd_vals = [0.00001, 0.0125, 0.025, 0.05, 0.1, 0.2];
   
   % Loop through the cases; defined by obs_intvl and time_err_sd
   for obs_intvl_indx = 1:5
      obs_intvl = obs_intvl_vals(obs_intvl_indx);
   
      for time_err_sd_indx = 1:6
         time_err_sd = time_err_sd_vals(time_err_sd_indx);
         % Don't do values that are too big
         if(time_err_sd <= obs_intvl / 100)

            % Read in the case header line
            str4 = fgetl(fid);

            %  Make sure that vales are the same
            r_obs_intvl = str2num(str4(10:13));
            r_time_err_sd = str2num(str4(27:end));
            if(r_obs_intvl ~= obs_intvl) stop; end
            if(r_time_err_sd ~= time_err_sd) stop; end

            % Read the 10 Trials
            for i = 1:10
               % First line has info about the truth trajctory and observations
               str_obs = fgetl(fid);
               obs(i, 1:5) = str2num(str_obs);
               % Second line has the results
               str_res = fgetl(fid);
               if(method == 29)
                  res(i, 1:7) = str2num(str_res);
               else
                  res(i, 1:8) = str2num(str_res);
               end

               % Stuff the values into the array for this entry
               rmse_line = rmse_line + 1;
               rmse_data(rmse_line, 1) = obs_intvl;
               rmse_data(rmse_line, 2) = time_err_sd;
               rmse_data(rmse_line, 3) = method;
               rmse_data(rmse_line, 4) = i;
               if(method == 14)
                  rmse_data(rmse_line, 5:12) = res(i, 1:8);
               else
                  rmse_data(rmse_line, 5:11) = res(i, 1:7);
                  rmse_data(rmse_line, 12) = -99;
               end

            end    % Of loop over 10 trials

            % Create the trial mean array, too
            mean_line = mean_line + 1;
            rmse_mean(mean_line, 1) = obs_intvl;
            rmse_mean(mean_line, 2) = time_err_sd;
            rmse_mean(mean_line, 3) = method;
            rmse_mean(mean_line, 4) = 0;
            for i = 5:12
               rmse_mean(mean_line, i) = mean(rmse_data(rmse_line-num_trials + 1 : rmse_line, i));
            end

         end   % Of if on legal case combination
      end    % Of time error sd loop
   end   % Of obs_interval loop
end   % Of method loop


for oi_index = 1:5
   oi_vals = [5 10 15 30 60];
   oi = oi_vals(oi_index);
   figure(oi_index);

   % Customize the axes
   if(oi_index == 1)
      ax = [0.005, 0.06, 0.18, 0.30];
      axt = [0.005, 0.06, 0.0, 0.025];
      off_val = 0.03;
      xt = [0.00625, 0.0125, 0.025, 0.05];
   elseif(oi_index ==2)
      ax = [0.005, 0.12, 0.2, 0.8];
      axt = [0.005, 0.12, 0.0, 0.03];
      off_val = 0.04;
      xt = [0.00625, 0.0125, 0.025, 0.05 0.1];
   elseif(oi_index ==3)
      ax = [0.005, 0.12, 0.3, 1.2];
      axt = [0.005, 0.12, 0.0, 0.04];
      off_val = 0.04;
      xt = [0.00625, 0.0125, 0.025, 0.05 0.1];
   elseif(oi_index ==4)
      ax = [0.005, 0.26, 0.7, 2.5];
      axt = [0.005, 0.26, 0.0, 0.15];
      off_val = 0.06;
      xt = [0.00625 0.0125, 0.025, 0.05 0.1 0.2];
   else
      ax = [0.005, 0.25, 1.8, 3.3];
      axt = [0.005, 0.25, 0, 0.2];
      off_val = 0.05;
      xt = [0.00625, 0.0125, 0.025, 0.05 0.1 0.2];
   end

   for method_indx = 1:5
      method_vals = [11 12 13 14 29];
      method = method_vals(method_indx);
      method_colors = ['g', 'b', 'r', 'k', 'm'];

      %--------------- Plot of average of trials -------
      indx = find(rmse_mean(:, 1) == oi & rmse_mean(:, 3) == method);
      y = rmse_mean(indx, :);
      %h = plot(y(:, 2), y(:, 7), '*');
      %set(h, 'color', method_colors(method_indx));
      %set(h, 'linewidth', 3);
      %hold on


      %--------------- Plot of all trials ------------
      
      for trial = 1:num_trials
         figure(oi_index);
         indx = find(rmse_data(:, 1) == oi & rmse_data(:, 3) == method & rmse_data(:, 4) == trial);
         x = rmse_data(indx, :);

         % Avoid that nasty log 0 problem by moving the 0 results have a value of 0.0625
         % Will need to note this in caption and relabel the tick mark
         x(:, 2) = max(x(:, 2), 0.00625);

         oi_offset = (method_indx - 3) * off_val;
         x_offset = x(:, 2) * exp(oi_offset);

         h = semilogx(x_offset, x(:, 7), '*');
         axis(ax);
         set(h, 'color', method_colors(method_indx));
         set(h, 'linewidth', 2);

         if(trial == 1) lh(oi_index, method_indx) = h; end

         set(gca, 'XTICK', xt);
         set(gca, 'linewidth', 2);
         set(gca, 'fontsize', 16);
         ylabel 'RMSE of Ensemble Mean';
         xlabel 'Time Error Standard Deviation';
         title(['Analysis Period ', num2str(oi * 0.01)]);
         hold on

         % Can also do separate plots of the errors in the estimate of the time offset
         figure(5 + oi_index);
         if(method == 14) 
            % For the impossible method, plot the impossible rmse of time offset
            hh = semilogx(x_offset, x(:, 12), '*'); 
         else
            hh = semilogx(x_offset, x(:, 11), '*'); 
         end
         axis(axt);
         set(hh, 'color', method_colors(method_indx));
         set(hh, 'linewidth', 2);

         if(trial == 1) lhh(oi_index, method_indx) = hh; end

         set(gca, 'XTICK', xt);
         set(gca, 'linewidth', 2);
         set(gca, 'fontsize', 16);
         xlabel 'Time Error Standard Deviation';
         ylabel 'Ensemble Mean Time Offset RMSE';
         title(['Analysis Period ', num2str(oi * 0.01)]);
         hold on

      end
      %--------------- Plot of all trials ------------

   end
   
end

% Slap on some legends; RMSE has 5 methods
for oi_index = 1:5
   legend(lh(oi_index, :), 'NOCORRECTION', 'LINEAR', 'VARONLY', 'IMPOSSIBLE', 'NONLINEAR', 'location', 'nw');
   % Slap on some legends; RMSE of time offset has 4 methods
   legend(lhh(oi_index, :), 'NOCORRECTION', 'LINEAR', 'VARONLY', 'IMPOSSIBLE', 'NONLINEAR', 'location', 'nw');
end

% Fix the 0 label
for i = 1:10 
   figure(i)
   xtl = get(gca, 'XTICKLABEL');
   xtl(1) = {'0.0'};
   set(gca, 'XTICKLABEL', xtl);
end

for i = 1:5
   figure(i);
   print(strcat('RMSE_Period_', num2str(obs_intvl_vals(i))), '-dpng');

   figure(i + 5);
   print(strcat('RMSE_Tau_', num2str(obs_intvl_vals(i))), '-dpng');
end




