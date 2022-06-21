function [state_incs, Rxy] = get_state_increments(state_ens, obs_ens, obs_incs)
%% get_state_increments Computes state increments given observation increments and
% the state and obs prior ensembles

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id: get_state_increments.m 11694 2017-06-02 21:39:25Z nancy@ucar.edu $

% Compute state variance and covariance
covar = cov(state_ens, obs_ens);

Rxy = covar(1, 2);

state_incs = obs_incs * covar(1, 2) / covar(2, 2);

% <next few lines under version control, do not edit>
% $URL: https://svn-dares-dart.cgd.ucar.edu/DART/releases/Manhattan/documentation/DART_LAB/matlab/private/get_state_increments.m $
% $Revision: 11694 $
% $Date: 2017-06-02 15:39:25 -0600 (Fri, 02 Jun 2017) $
