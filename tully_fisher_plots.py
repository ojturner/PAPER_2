import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
from astropy.io import fits
from scipy.stats import binned_statistic
from matplotlib import rc
import matplotlib.ticker as ticker
from lmfit import Model, Parameters
from matplotlib.colors import LogNorm

global_beta = 2.5

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', weight='bold')
rc('text', usetex=True)
rc('axes', linewidth=2)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# read in the tables
table_cresci_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/cresci_2009.cat')
table_dynamo_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/dynamo_kinematics.txt')
table_massiv_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/massiv_kinematics_all.txt')
table_massiv_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/massiv_kinematics_rot.txt')
table_massiv_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/massiv_kinematics_disp.txt')
table_swinbank_low_redshift_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_low_redshift_kinematics_all.txt')
table_swinbank_low_redshift_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_low_redshift_kinematics_rot.txt')
table_swinbank_low_redshift_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_low_redshift_kinematics_disp.txt')
table_swinbank_high_redshift_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_high_redshift_kinematics_all.txt')
table_swinbank_high_redshift_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_high_redshift_kinematics_rot.txt')
table_swinbank_high_redshift_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/swinbank_high_redshift_kinematics_disp.txt')
table_sigma_low_redshift_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_low_redshift_kinematics_all.txt')
table_sigma_low_redshift_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_low_redshift_kinematics_rot.txt')
table_sigma_low_redshift_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_low_redshift_kinematics_disp.txt')
table_sigma_high_redshift_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_high_redshift_kinematics_all.txt')
table_sigma_high_redshift_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_high_redshift_kinematics_rot.txt')
table_sigma_high_redshift_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/sigma_high_redshift_kinematics_disp.txt')
table_amaze_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/gnerucci_clean.txt')
table_kross_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/kross_for_paper_2_all.cat')
table_kross_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/kross_for_paper_2_rot.cat')
table_kross_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/kross_for_paper_2_disp.cat')
table_kds_all = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/v_over_sigma_evolution_plots/TABLES/joint_isolated_rotator_properties_new_mass_PRE_REFEREE.cat')
table_kds_rot = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/v_over_sigma_evolution_plots/TABLES/tf_joint_vc_gt_1_new_mass.cat')
table_kds_disp = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/v_over_sigma_evolution_plots/TABLES/tf_joint_vc_lt_1_new_mass.cat')


# define the convenience functions used throughout this script

def line(x, slope, intercept):
    "line"
    return slope * x + intercept

def tf_fit(mass_10,
           log_v,
           log_v_error,
           slope,
           slope_error):

    """
    Perturb by errors each time - fit a fixed slope relation to the points
    in the TF plane to look at the change in normalisation of the relation

    mass_10 - the log mass - 10
    log_v - the log velocity in kms-1
    log_v_error - average error on the velocity to perturb by
    """
        
    tf_x = []
    for entry in mass_10:
        tf_x.append(np.random.normal(entry,0.2))
    tf_x = np.array(tf_x)

    tf_y = []
    for val,err in zip(log_v,log_v_error):
        tf_y.append(np.random.normal(val,err))
    tf_y = np.array(tf_y)

    fit_slope = np.random.normal(slope,slope_error)

    # define the fractional error as a better weighting scheme
    velocity = np.power(10, log_v)
    velocity_error = (velocity * log_v_error) / 0.434
    fractional_error = velocity_error / velocity

    mod = Model(line)
    pars = Parameters()
    pars.add('slope', value=fit_slope, vary=False)
    pars.add('intercept', value=2.14, vary=True)
    result = mod.fit(tf_y, pars, weights=(1./fractional_error), x=tf_x)
    return result.best_values, result.redchi

def tf_fit_1000(mass_10,
                log_v,
                log_v_error,
                slope,
                slope_error):

    """
    Run the tf_fit 1000 times and record the parameter distribution
    return the median, 16th and 84th percentiles.

    mass_10 - the log mass - 10
    log_v - the log velocity in kms-1
    log_v_error - average error on the velocity to perturb by
    """

    param = []
    redchi = []
    for i in range(100):
        parameters, chi = tf_fit(mass_10,
                                 log_v,
                                 log_v_error,
                                 slope,
                                 slope_error)
        param.append(parameters['intercept'])
        redchi.append(chi)

    number = len(log_v)

    return [np.median(param),
            np.percentile(param,16),
            np.percentile(param,84),
            np.median(redchi),
            number]

def ang_fit(mass_10,
            log_ang,
            log_ang_error,
            slope,
            slope_error):

    """
    Perturb by errors each time - fit a fixed slope relation to the points
    in the TF plane to look at the change in normalisation of the relation

    mass_10 - the log mass - 10
    log_v - the log velocity in kms-1
    log_v_error - average error on the velocity to perturb by
    """
        
    tf_x = []
    for entry in mass_10:
        tf_x.append(np.random.normal(entry,0.2))
    tf_x = np.array(tf_x)

    tf_y = []
    for val,err in zip(log_ang,log_ang_error):
        tf_y.append(np.random.normal(val,err))
    tf_y = np.array(tf_y)

    fit_slope = np.random.normal(slope,slope_error)

    # define the fractional error weights
    ang_mom = np.power(10, log_ang)
    ang_error = (ang_mom * log_ang_error) / 0.434
    fractional_error = ang_error / ang_mom

    mod = Model(line)
    pars = Parameters()
    pars.add('slope', value=fit_slope, vary=False)
    pars.add('intercept', value=2.80, vary=True)
    result = mod.fit(tf_y, pars, weights=(1./fractional_error), x=tf_x)
    return result.best_values, result.redchi

def ang_fit_1000(mass_10,
                 log_ang,
                 log_ang_error,
                 slope,
                 slope_error):

    """
    Run the tf_fit 1000 times and record the parameter distribution
    return the median, 16th and 84th percentiles.

    mass_10 - the log mass - 10
    log_v - the log velocity in kms-1
    log_v_error - average error on the velocity to perturb by
    """

    param = []
    redchi = []
    for i in range(100):
        parameters, chi = ang_fit(mass_10,
                                  log_ang,
                                  log_ang_error,
                                  slope,
                                  slope_error)
        param.append(parameters['intercept'])
        redchi.append(chi)

    number = len(log_ang)

    return [np.median(param),
            np.percentile(param,16),
            np.percentile(param,84),
            np.median(redchi),
            number]

def return_betas(mass,
                 velocity,
                 velocity_upper_error,
                 velocity_lower_error,
                 sigma,
                 sigma_upper_error,
                 sigma_lower_error,
                 radius,
                 beta=2.5,
                 alpha_velocity=0.28,
                 alpha_velocity_error=0.1,
                 alpha_ang_mom=0.55,
                 alpha_ang_mom_error=0.1):
    """
    Def: Spit out all the parameters required so that we don't
    need to do all the calculations many times

    Outputs: fit_mass
             log_velocity
             log_velocity_upper_error
             log_velocity_lower_error
             v_tot
             v_tot_lower_error
             v_tot_upper_error
             log_v_tot
             log_v_tot_lower_error
             log_v_tot_upper_error
             log_v_tot_av_error
             ang_mom
             ang_mom_lower_error
             ang_mom_upper_error
             log_ang_mom
             log_ang_mom_lower_error
             log_ang_mom_upper_error
             log_ang_mom_av_error
             ang_mom_tot
             ang_mom_tot_lower_error
             ang_mom_tot_upper_error
             log_ang_mom_tot
             log_ang_mom_tot_lower_error
             log_ang_mom_tot_upper_error
             log_ang_mom_tot_av_error
             log_velocity_mass_fit_params
             log_velocity_mass_fit_line
             log_velocity_mass_fit_line_lower
             log_velocity_mass_fit_line_upper
             log_v_tot_mass_fit_params
             log_v_tot_mass_fit_line
             log_v_tot_mass_fit_line_lower
             log_v_tot_mass_fit_line_upper
             log_ang_mom_mass_fit_params
             log_ang_mom_mass_fit_line
             log_ang_mom_mass_fit_line_lower
             log_ang_mom_mass_fit_line_upper
             log_ang_mom_tot_mass_fit_params
             log_ang_mom_tot_mass_fit_line
             log_ang_mom_tot_mass_fit_line_lower
             log_ang_mom_tot_mass_fit_line_lower


    """

    # define the mass baseline for reconstructing the fit lines

    fit_mass = np.arange(8.0,11.4,0.2)

    # define the mass for fitting

    mass_10 = mass - 10.10

    log_velocity = np.log10(velocity)
    log_velocity_upper_error = 0.434 * (1.0*velocity_upper_error/velocity)
    log_velocity_lower_error = 0.434 * (1.0*velocity_lower_error/velocity)
    log_velocity_av_error = 0.5 * (log_velocity_upper_error + log_velocity_lower_error)

    v_tot = np.sqrt((velocity)**2 + (beta/radius)*(sigma)**2)
    v_tot_upper_error = np.sqrt((1.0*velocity*velocity_upper_error/v_tot)**2 + ((beta/radius) * sigma * sigma_upper_error / v_tot)**2)
    v_tot_lower_error = np.sqrt((1.0*velocity*velocity_lower_error/v_tot)**2 + ((beta/radius) * sigma * sigma_lower_error / v_tot)**2)
    log_v_tot = np.log10(v_tot)
    log_v_tot_upper_error = 0.434 * (1.0*v_tot_upper_error/v_tot)
    log_v_tot_lower_error = 0.434 * (1.0*v_tot_lower_error/v_tot)
    log_v_tot_av_error = 0.5 * (log_v_tot_upper_error + log_v_tot_lower_error)

    ang_mom = 1.19 * radius * velocity
    ang_mom_upper_error = 1.19 * radius * velocity_upper_error
    ang_mom_lower_error = 1.19 * radius * velocity_lower_error
    log_ang_mom = np.log10(ang_mom)
    log_ang_mom_upper_error = 0.434 * (1.0*ang_mom_upper_error/ang_mom)
    log_ang_mom_lower_error = 0.434 * (1.0*ang_mom_lower_error/ang_mom)
    log_ang_mom_av_error = 0.5 * (log_ang_mom_upper_error + log_ang_mom_lower_error)

    ang_mom_tot = 1.19 * radius * v_tot
    ang_mom_tot_upper_error = 1.19 * radius * v_tot_upper_error
    ang_mom_tot_lower_error = 1.19 * radius * v_tot_lower_error
    log_ang_mom_tot = np.log10(ang_mom_tot)
    log_ang_mom_tot_upper_error = 0.434 * (1.0*ang_mom_tot_upper_error/ang_mom_tot)
    log_ang_mom_tot_lower_error = 0.434 * (1.0*ang_mom_tot_lower_error/ang_mom_tot)
    log_ang_mom_tot_av_error = 0.5 * (log_ang_mom_tot_upper_error + log_ang_mom_tot_lower_error)

    # dynamical mass and ratio 
    dynamical_mass = (2 * 3.086E19 * radius * (velocity * 1000)**2) / 1.326E20
    dynamical_mass_upper_error = np.sqrt((4 * radius * 4.6267E35 * velocity * velocity_upper_error)**2)/1.9989E30
    dynamical_mass_lower_error = np.sqrt((4 * radius * 4.6267E35 * velocity * velocity_lower_error)**2)/1.9989E30
    log_dynamical_mass = np.log10(dynamical_mass)
    log_dynamical_mass_upper_error = 0.434 * (dynamical_mass_upper_error/dynamical_mass)
    log_dynamical_mass_lower_error = 0.434 * (dynamical_mass_lower_error/dynamical_mass)
    dynamical_mass_ratio = dynamical_mass / np.power(10, mass)
    median_dynamical_mass_ratio = np.median(dynamical_mass_ratio)
    median_dynamical_mass_ratio_lower_error = np.percentile(dynamical_mass_ratio, 16)
    median_dynamical_mass_ratio_upper_error = np.percentile(dynamical_mass_ratio, 84)
    mean_dynamical_mass_ratio = np.mean(dynamical_mass_ratio)
    mean_dynamical_mass_ratio_error = np.std(dynamical_mass_ratio)/np.sqrt(len(dynamical_mass_ratio))

    # dynamical mass with sigma and ratio 
    dynamical_mass_with_sigma = (2 * radius * 3.086E19 * ((velocity * 1000)**2 + (beta/radius)*(sigma*1000)**2)) / 1.3267E20
    dynamical_mass_with_sigma_upper_error = np.sqrt((4 * radius * 4.6267E35 * velocity * velocity_upper_error)**2 + (4 * (beta/radius) * radius * 4.6267E35 * sigma * sigma_upper_error)**2)/1.9989E30
    dynamical_mass_with_sigma_lower_error = np.sqrt((4 * radius * 4.6267E35 * velocity * velocity_lower_error)**2 + (4 * (beta/radius) * radius * 4.6267E35 * sigma * sigma_lower_error)**2)/1.9989E30
    log_dynamical_mass_with_sigma = np.log10(dynamical_mass_with_sigma)
    log_dynamical_mass_with_sigma_upper_error = 0.434 * (dynamical_mass_with_sigma_upper_error/dynamical_mass_with_sigma)
    log_dynamical_mass_with_sigma_lower_error = 0.434 * (dynamical_mass_with_sigma_lower_error/dynamical_mass_with_sigma)
    dynamical_mass_with_sigma_ratio = dynamical_mass_with_sigma / np.power(10, mass)
    median_dynamical_mass_with_sigma_ratio = np.median(dynamical_mass_with_sigma_ratio)
    median_dynamical_mass_with_sigma_ratio_lower_error = np.percentile(dynamical_mass_with_sigma_ratio, 16)
    median_dynamical_mass_with_sigma_ratio_upper_error = np.percentile(dynamical_mass_with_sigma_ratio, 84)
    mean_dynamical_mass_with_sigma_ratio = np.mean(dynamical_mass_with_sigma_ratio)
    mean_dynamical_mass_with_sigma_ratio_error = np.std(dynamical_mass_with_sigma_ratio)/np.sqrt(len(dynamical_mass_with_sigma_ratio))

    # now use the defined methods to do the fitting to each of these relations

    log_velocity_mass_fit_params = tf_fit_1000(mass_10,
                                               log_velocity,
                                               log_velocity_av_error,
                                               alpha_velocity,
                                               alpha_velocity_error)


    log_velocity_mass_fit_line = log_velocity_mass_fit_params[0] + (alpha_velocity*(fit_mass - 10.10))
    log_velocity_mass_fit_line_lower = log_velocity_mass_fit_params[1] + (alpha_velocity*(fit_mass - 10.10))
    log_velocity_mass_fit_line_upper = log_velocity_mass_fit_params[2] + (alpha_velocity*(fit_mass - 10.10))

    log_v_tot_mass_fit_params = tf_fit_1000(mass_10,
                                            log_v_tot,
                                            log_v_tot_av_error,
                                            alpha_velocity,
                                            alpha_velocity_error)

    log_v_tot_mass_fit_line = log_v_tot_mass_fit_params[0] + (alpha_velocity*(fit_mass - 10.10))
    log_v_tot_mass_fit_line_lower = log_v_tot_mass_fit_params[1] + (alpha_velocity*(fit_mass - 10.10))
    log_v_tot_mass_fit_line_upper = log_v_tot_mass_fit_params[2] + (alpha_velocity*(fit_mass - 10.10))

    log_ang_mom_mass_fit_params = ang_fit_1000(mass_10,
                                               log_ang_mom,
                                               log_ang_mom_av_error,
                                               alpha_ang_mom,
                                               alpha_ang_mom_error)

    log_ang_mom_mass_fit_line = log_ang_mom_mass_fit_params[0] + (alpha_ang_mom*(fit_mass - 10.10))
    log_ang_mom_mass_fit_line_lower = log_ang_mom_mass_fit_params[1] + (alpha_ang_mom*(fit_mass - 10.10))
    log_ang_mom_mass_fit_line_upper = log_ang_mom_mass_fit_params[2] + (alpha_ang_mom*(fit_mass - 10.10))

    log_ang_mom_tot_mass_fit_params = ang_fit_1000(mass_10,
                                                   log_ang_mom_tot,
                                                   log_ang_mom_tot_av_error,
                                                   alpha_ang_mom,
                                                   alpha_ang_mom_error)

    log_ang_mom_tot_mass_fit_line = log_ang_mom_tot_mass_fit_params[0] + (alpha_ang_mom*(fit_mass - 10.10))
    log_ang_mom_tot_mass_fit_line_lower = log_ang_mom_tot_mass_fit_params[1] + (alpha_ang_mom*(fit_mass - 10.10))
    log_ang_mom_tot_mass_fit_line_upper = log_ang_mom_tot_mass_fit_params[2] + (alpha_ang_mom*(fit_mass - 10.10))

    return { 'fit_mass': fit_mass,
             'log_velocity': log_velocity,
             'log_velocity_upper_error': log_velocity_upper_error,
             'log_velocity_lower_error': log_velocity_lower_error,
             'v_tot': v_tot,
             'v_tot_lower_error': v_tot_lower_error,
             'v_tot_upper_error': v_tot_upper_error,
             'log_v_tot': log_v_tot,
             'log_v_tot_lower_error': log_v_tot_lower_error,
             'log_v_tot_upper_error': log_v_tot_upper_error,
             'log_v_tot_av_error': log_v_tot_av_error,
             'ang_mom': ang_mom,
             'ang_mom_lower_error': ang_mom_lower_error,
             'ang_mom_upper_error': ang_mom_upper_error,
             'log_ang_mom': log_ang_mom,
             'log_ang_mom_lower_error': log_ang_mom_lower_error,
             'log_ang_mom_upper_error': log_ang_mom_upper_error,
             'log_ang_mom_av_error': log_ang_mom_av_error,
             'ang_mom_tot': ang_mom_tot,
             'ang_mom_tot_lower_error': ang_mom_tot_lower_error,
             'ang_mom_tot_upper_error': ang_mom_tot_upper_error,
             'log_ang_mom_tot': log_ang_mom_tot,
             'log_ang_mom_tot_lower_error': log_ang_mom_tot_lower_error,
             'log_ang_mom_tot_upper_error': log_ang_mom_tot_upper_error,
             'log_ang_mom_tot_av_error': log_ang_mom_tot_av_error,
             'log_velocity_mass_fit_params': log_velocity_mass_fit_params,
             'log_velocity_mass_fit_line': log_velocity_mass_fit_line,
             'log_velocity_mass_fit_line_lower': log_velocity_mass_fit_line_lower,
             'log_velocity_mass_fit_line_upper': log_velocity_mass_fit_line_upper,
             'log_v_tot_mass_fit_params': log_v_tot_mass_fit_params,
             'log_v_tot_mass_fit_line': log_v_tot_mass_fit_line,
             'log_v_tot_mass_fit_line_lower': log_v_tot_mass_fit_line_lower,
             'log_v_tot_mass_fit_line_upper': log_v_tot_mass_fit_line_upper,
             'log_ang_mom_mass_fit_params': log_ang_mom_mass_fit_params,
             'log_ang_mom_mass_fit_line': log_ang_mom_mass_fit_line,
             'log_ang_mom_mass_fit_line_lower': log_ang_mom_mass_fit_line_lower,
             'log_ang_mom_mass_fit_line_upper': log_ang_mom_mass_fit_line_upper,
             'log_ang_mom_tot_mass_fit_params': log_ang_mom_tot_mass_fit_params,
             'log_ang_mom_tot_mass_fit_line': log_ang_mom_tot_mass_fit_line,
             'log_ang_mom_tot_mass_fit_line_lower': log_ang_mom_tot_mass_fit_line_lower,
             'log_ang_mom_tot_mass_fit_line_upper': log_ang_mom_tot_mass_fit_line_upper,
             'median_dynamical_mass_ratio': median_dynamical_mass_ratio,
             'median_dynamical_mass_ratio_lower_error': median_dynamical_mass_ratio_lower_error,
             'median_dynamical_mass_ratio_upper_error': median_dynamical_mass_ratio_upper_error,
             'median_dynamical_mass_with_sigma_ratio': median_dynamical_mass_with_sigma_ratio,
             'median_dynamical_mass_with_sigma_ratio_lower_error': median_dynamical_mass_with_sigma_ratio_lower_error,
             'median_dynamical_mass_with_sigma_ratio_upper_error': median_dynamical_mass_with_sigma_ratio_upper_error}

# some lines from the literature
kassin_mass = np.arange(8.0, 11.6, 0.1)
kassin_s0 = (0.29*kassin_mass)
ten_index = np.argmin(abs(10 - kassin_mass))
kassin_s0 = (kassin_s0 - kassin_s0[ten_index]) + 1.89

tiley_v = np.arange(1.0, 3.0, 0.1)
tiley_mass = (4.7 * (tiley_v - 2.25)) + 10.0

tiley_z0_v = np.arange(1.0, 3.0, 0.1)
tiley_z0_mass = (3.68 * (tiley_v - 2.09)) + 9.83

reyes_tf_mass = np.arange(8.0,11.4,0.2)
log_reyes_tf_velocity = 2.142 + (0.278*(reyes_tf_mass - 10.1))
log_fall_ang_mom = 2.89 + (0.51*(reyes_tf_mass - 10.1))

log_fall_ang_mom_ev = 2.83 + (0.6*(reyes_tf_mass - 10.1))
ang_mom_ev = np.log10(np.power(10,log_fall_ang_mom_ev) * (1.9)**(-0.67))

# going to automatically write out to each of the evolution.txt files
# will allow for quick exploration of different values of beta and 
# different choices for the fitting, whilst still evaluating all
# initialise here the lists which will be used to create these files

velocity_evolution_list = [['# Survey',
                           'redshift',
                           'N_all',
                           'chi_all',
                           'beta_all',
                           'delta_beta_all',
                           'delta_beta_all_lower_error',
                           'delta_beta_all_upper_error',
                           'N_rot',
                           'chi_rot',
                           'beta_rot',
                           'delta_beta_rot',
                           'delta_beta_rot_lower_error',
                           'delta_beta_rot_upper_error',
                           'N_disp',
                           'chi_disp',
                           'beta_disp',
                           'delta_beta_disp',
                           'delta_beta_disp_lower_error',
                           'delta_beta_disp_upper_error']]

v_tot_evolution_list = [['# Survey',
                       'redshift',
                       'N_all',
                       'chi_all',
                       'beta_all',
                       'delta_beta_all',
                       'delta_beta_all_lower_error',
                       'delta_beta_all_upper_error',
                       'N_rot',
                       'chi_rot',
                       'beta_rot',
                       'delta_beta_rot',
                       'delta_beta_rot_lower_error',
                       'delta_beta_rot_upper_error',
                       'N_disp',
                       'chi_disp',
                       'beta_disp',
                       'delta_beta_disp',
                       'delta_beta_disp_lower_error',
                       'delta_beta_disp_upper_error']]

ang_mom_evolution_list = [['# Survey',
                       'redshift',
                       'N_all',
                       'chi_all',
                       'beta_all',
                       'delta_beta_all',
                       'delta_beta_all_lower_error',
                       'delta_beta_all_upper_error',
                       'N_rot',
                       'chi_rot',
                       'beta_rot',
                       'delta_beta_rot',
                       'delta_beta_rot_lower_error',
                       'delta_beta_rot_upper_error',
                       'N_disp',
                       'chi_disp',
                       'beta_disp',
                       'delta_beta_disp',
                       'delta_beta_disp_lower_error',
                       'delta_beta_disp_upper_error']]

ang_mom_tot_evolution_list = [['# Survey',
                           'redshift',
                           'N_all',
                           'chi_all',
                           'beta_all',
                           'delta_beta_all',
                           'delta_beta_all_lower_error',
                           'delta_beta_all_upper_error',
                           'N_rot',
                           'chi_rot',
                           'beta_rot',
                           'delta_beta_rot',
                           'delta_beta_rot_lower_error',
                           'delta_beta_rot_upper_error',
                           'N_disp',
                           'chi_disp',
                           'beta_disp',
                           'delta_beta_disp',
                           'delta_beta_disp_lower_error',
                           'delta_beta_disp_upper_error']]

dyn_mass_evolution_list = [['# Survey',
                            'redshift',
                            'N_all',
                            'dyn_ratio_all',
                            'dyn_ratio_all_lower_error',
                            'dyn_ratio_all_upper_error',
                            'dyn_ratio_with_sigma_all',
                            'dyn_ratio_with_sigma_all_lower_error',
                            'dyn_ratio_with_sigma_all_upper_error',
                            'N_rot',
                            'dyn_ratio_rot',
                            'dyn_ratio_rot_lower_error',
                            'dyn_ratio_rot_upper_error',
                            'dyn_ratio_with_sigma_rot',
                            'dyn_ratio_with_sigma_rot_lower_error',
                            'dyn_ratio_with_sigma_rot_upper_error',
                            'N_disp',
                            'dyn_ratio_disp',
                            'dyn_ratio_disp_lower_error',
                            'dyn_ratio_disp_upper_error',
                            'dyn_ratio_with_sigma_disp',
                            'dyn_ratio_with_sigma_disp_lower_error',
                            'dyn_ratio_with_sigma_disp_upper_error']]

# after each of the survey evaluations, population a separate list
# with all of the entries required. Also define here the redshift zero
# betas that are being compared against

velocity_beta_zero = 2.02
velocity_beta_zero_error = 0.02
v_tot_beta_zero = 2.02
v_tot_beta_zero_error = 0.02
ang_mom_beta_zero = 2.70
ang_mom_beta_zero_error = 0.02
ang_mom_tot_beta_zero = 2.70
ang_mom_tot_beta_zero_error = 0.02

# DYNAMICAL PROPERTIES FOR CRESCI

cresci_all_velocity = table_cresci_all['Vmax']
cresci_all_velocity_upper_error = table_cresci_all['Vmax_err_lower']
cresci_all_velocity_lower_error = table_cresci_all['Vmax_err_upper']
cresci_all_mass = table_cresci_all['Mstar']
cresci_all_radius = table_cresci_all['Rd']
cresci_all_sigma = table_cresci_all['sigma']
cresci_all_sigma_upper_error = table_cresci_all['sigma_upper_error']
cresci_all_sigma_lower_error = table_cresci_all['sigma_lower_error']

cresci_all_dictionary = return_betas(cresci_all_mass,
                                    cresci_all_velocity,
                                    cresci_all_velocity_upper_error,
                                    cresci_all_velocity_lower_error,
                                    cresci_all_sigma,
                                    cresci_all_sigma_upper_error,
                                    cresci_all_sigma_lower_error,
                                    cresci_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = cresci_all_dictionary['fit_mass']
cresci_all_log_velocity = cresci_all_dictionary['log_velocity']
cresci_all_log_velocity_upper_error = cresci_all_dictionary['log_velocity_upper_error']
cresci_all_log_velocity_lower_error = cresci_all_dictionary['log_velocity_lower_error']
cresci_all_v_tot = cresci_all_dictionary['v_tot']
cresci_all_v_tot_lower_error = cresci_all_dictionary['v_tot_lower_error']
cresci_all_v_tot_upper_error = cresci_all_dictionary['v_tot_upper_error']
cresci_all_log_v_tot = cresci_all_dictionary['log_v_tot']
cresci_all_log_v_tot_lower_error = cresci_all_dictionary['log_v_tot_lower_error']
cresci_all_log_v_tot_upper_error = cresci_all_dictionary['log_v_tot_upper_error']
cresci_all_log_v_tot_av_error = cresci_all_dictionary['log_v_tot_av_error']
cresci_all_ang_mom = cresci_all_dictionary['ang_mom']
cresci_all_ang_mom_lower_error = cresci_all_dictionary['ang_mom_lower_error']
cresci_all_ang_mom_upper_error = cresci_all_dictionary['ang_mom_upper_error']
cresci_all_log_ang_mom = cresci_all_dictionary['log_ang_mom']
cresci_all_log_ang_mom_lower_error = cresci_all_dictionary['log_ang_mom_lower_error']
cresci_all_log_ang_mom_upper_error = cresci_all_dictionary['log_ang_mom_upper_error']
cresci_all_log_ang_mom_av_error = cresci_all_dictionary['log_ang_mom_av_error']
cresci_all_ang_mom_tot = cresci_all_dictionary['ang_mom_tot']
cresci_all_ang_mom_tot_lower_error = cresci_all_dictionary['ang_mom_tot_lower_error']
cresci_all_ang_mom_tot_upper_error = cresci_all_dictionary['ang_mom_tot_upper_error']
cresci_all_log_ang_mom_tot = cresci_all_dictionary['log_ang_mom_tot']
cresci_all_log_ang_mom_tot_lower_error = cresci_all_dictionary['log_ang_mom_tot_lower_error']
cresci_all_log_ang_mom_tot_upper_error = cresci_all_dictionary['log_ang_mom_tot_upper_error']
cresci_all_log_ang_mom_tot_av_error = cresci_all_dictionary['log_ang_mom_tot_av_error']
cresci_all_log_velocity_mass_fit_params = cresci_all_dictionary['log_velocity_mass_fit_params']
print 'cresci_all velocity vs. mass fit params: %s' % cresci_all_log_velocity_mass_fit_params
cresci_all_log_velocity_mass_fit_line = cresci_all_dictionary['log_velocity_mass_fit_line']
cresci_all_log_velocity_mass_fit_line_lower = cresci_all_dictionary['log_velocity_mass_fit_line_lower']
cresci_all_log_velocity_mass_fit_line_upper = cresci_all_dictionary['log_velocity_mass_fit_line_upper']
cresci_all_log_v_tot_mass_fit_params = cresci_all_dictionary['log_v_tot_mass_fit_params']
print 'cresci_all v_tot vs. mass fit params: %s' % cresci_all_log_v_tot_mass_fit_params
cresci_all_log_v_tot_mass_fit_line = cresci_all_dictionary['log_v_tot_mass_fit_line']
cresci_all_log_v_tot_mass_fit_line_lower = cresci_all_dictionary['log_v_tot_mass_fit_line_lower']
cresci_all_log_v_tot_mass_fit_line_upper = cresci_all_dictionary['log_v_tot_mass_fit_line_upper']
cresci_all_log_ang_mom_mass_fit_params = cresci_all_dictionary['log_ang_mom_mass_fit_params']
print 'cresci_all ang_mom vs. mass fit params: %s' % cresci_all_log_ang_mom_mass_fit_params
cresci_all_log_ang_mom_mass_fit_line = cresci_all_dictionary['log_ang_mom_mass_fit_line']
cresci_all_log_ang_mom_mass_fit_line_lower = cresci_all_dictionary['log_ang_mom_mass_fit_line_lower']
cresci_all_log_ang_mom_mass_fit_line_upper = cresci_all_dictionary['log_ang_mom_mass_fit_line_upper']
cresci_all_log_ang_mom_tot_mass_fit_params = cresci_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'cresci_all ang_mom_tot vs. mass fit params: %s' % cresci_all_log_ang_mom_tot_mass_fit_params
cresci_all_log_ang_mom_tot_mass_fit_line = cresci_all_dictionary['log_ang_mom_tot_mass_fit_line']
cresci_all_log_ang_mom_tot_mass_fit_line_lower = cresci_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
cresci_all_log_ang_mom_tot_mass_fit_line_upper = cresci_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
cresci_all_median_dynamical_mass_ratio = cresci_all_dictionary['median_dynamical_mass_ratio']
cresci_all_median_dynamical_mass_ratio_lower_error = cresci_all_dictionary['median_dynamical_mass_ratio_lower_error']
cresci_all_median_dynamical_mass_ratio_upper_error = cresci_all_dictionary['median_dynamical_mass_ratio_upper_error']
cresci_all_median_dynamical_mass_with_sigma_ratio = cresci_all_dictionary['median_dynamical_mass_with_sigma_ratio']
cresci_all_median_dynamical_mass_with_sigma_ratio_lower_error = cresci_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
cresci_all_median_dynamical_mass_with_sigma_ratio_upper_error = cresci_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['cresci',
                                2.0,
                                cresci_all_log_velocity_mass_fit_params[4],
                                cresci_all_log_velocity_mass_fit_params[3],
                                cresci_all_log_velocity_mass_fit_params[0],
                                cresci_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((cresci_all_log_velocity_mass_fit_params[0] - cresci_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((cresci_all_log_velocity_mass_fit_params[2] - cresci_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
v_tot_evolution_list.append(['cresci',
                                2.0,
                                cresci_all_log_v_tot_mass_fit_params[4],
                                cresci_all_log_v_tot_mass_fit_params[3],
                                cresci_all_log_v_tot_mass_fit_params[0],
                                cresci_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((cresci_all_log_v_tot_mass_fit_params[0] - cresci_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((cresci_all_log_v_tot_mass_fit_params[2] - cresci_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_evolution_list.append(['cresci',
                                2.0,
                                cresci_all_log_ang_mom_mass_fit_params[4],
                                cresci_all_log_ang_mom_mass_fit_params[3],
                                cresci_all_log_ang_mom_mass_fit_params[0],
                                cresci_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((cresci_all_log_ang_mom_mass_fit_params[0] - cresci_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((cresci_all_log_ang_mom_mass_fit_params[2] - cresci_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_tot_evolution_list.append(['cresci',
                                2.0,
                                cresci_all_log_ang_mom_tot_mass_fit_params[4],
                                cresci_all_log_ang_mom_tot_mass_fit_params[3],
                                cresci_all_log_ang_mom_tot_mass_fit_params[0],
                                cresci_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((cresci_all_log_ang_mom_tot_mass_fit_params[0] - cresci_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((cresci_all_log_ang_mom_tot_mass_fit_params[2] - cresci_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
dyn_mass_evolution_list.append(['cresci',
                                2.0,
                                cresci_all_log_velocity_mass_fit_params[4],
                                cresci_all_median_dynamical_mass_ratio,
                                cresci_all_median_dynamical_mass_ratio_lower_error,
                                cresci_all_median_dynamical_mass_ratio_upper_error,
                                cresci_all_median_dynamical_mass_with_sigma_ratio,
                                cresci_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                cresci_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])



# DYNAMICAL PROPERTIES FOR SIGMA LOW REDSHIFT ALL

sigma_low_redshift_all_velocity = table_sigma_low_redshift_all['Vrot']
sigma_low_redshift_all_velocity_upper_error = table_sigma_low_redshift_all['Vrot_err']
sigma_low_redshift_all_velocity_lower_error = table_sigma_low_redshift_all['Vrot_err']
sigma_low_redshift_all_mass = table_sigma_low_redshift_all['Mstar']
sigma_low_redshift_all_radius = table_sigma_low_redshift_all['radius']
sigma_low_redshift_all_sigma = table_sigma_low_redshift_all['sigma']
sigma_low_redshift_all_sigma_upper_error = table_sigma_low_redshift_all['sigma_err']
sigma_low_redshift_all_sigma_lower_error = table_sigma_low_redshift_all['sigma_err']

sigma_low_redshift_all_dictionary = return_betas(sigma_low_redshift_all_mass,
                                    sigma_low_redshift_all_velocity,
                                    sigma_low_redshift_all_velocity_upper_error,
                                    sigma_low_redshift_all_velocity_lower_error,
                                    sigma_low_redshift_all_sigma,
                                    sigma_low_redshift_all_sigma_upper_error,
                                    sigma_low_redshift_all_sigma_lower_error,
                                    sigma_low_redshift_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_low_redshift_all_dictionary['fit_mass']
sigma_low_redshift_all_log_velocity = sigma_low_redshift_all_dictionary['log_velocity']
sigma_low_redshift_all_log_velocity_upper_error = sigma_low_redshift_all_dictionary['log_velocity_upper_error']
sigma_low_redshift_all_log_velocity_lower_error = sigma_low_redshift_all_dictionary['log_velocity_lower_error']
sigma_low_redshift_all_v_tot = sigma_low_redshift_all_dictionary['v_tot']
sigma_low_redshift_all_v_tot_lower_error = sigma_low_redshift_all_dictionary['v_tot_lower_error']
sigma_low_redshift_all_v_tot_upper_error = sigma_low_redshift_all_dictionary['v_tot_upper_error']
sigma_low_redshift_all_log_v_tot = sigma_low_redshift_all_dictionary['log_v_tot']
sigma_low_redshift_all_log_v_tot_lower_error = sigma_low_redshift_all_dictionary['log_v_tot_lower_error']
sigma_low_redshift_all_log_v_tot_upper_error = sigma_low_redshift_all_dictionary['log_v_tot_upper_error']
sigma_low_redshift_all_log_v_tot_av_error = sigma_low_redshift_all_dictionary['log_v_tot_av_error']
sigma_low_redshift_all_ang_mom = sigma_low_redshift_all_dictionary['ang_mom']
sigma_low_redshift_all_ang_mom_lower_error = sigma_low_redshift_all_dictionary['ang_mom_lower_error']
sigma_low_redshift_all_ang_mom_upper_error = sigma_low_redshift_all_dictionary['ang_mom_upper_error']
sigma_low_redshift_all_log_ang_mom = sigma_low_redshift_all_dictionary['log_ang_mom']
sigma_low_redshift_all_log_ang_mom_lower_error = sigma_low_redshift_all_dictionary['log_ang_mom_lower_error']
sigma_low_redshift_all_log_ang_mom_upper_error = sigma_low_redshift_all_dictionary['log_ang_mom_upper_error']
sigma_low_redshift_all_log_ang_mom_av_error = sigma_low_redshift_all_dictionary['log_ang_mom_av_error']
sigma_low_redshift_all_ang_mom_tot = sigma_low_redshift_all_dictionary['ang_mom_tot']
sigma_low_redshift_all_ang_mom_tot_lower_error = sigma_low_redshift_all_dictionary['ang_mom_tot_lower_error']
sigma_low_redshift_all_ang_mom_tot_upper_error = sigma_low_redshift_all_dictionary['ang_mom_tot_upper_error']
sigma_low_redshift_all_log_ang_mom_tot = sigma_low_redshift_all_dictionary['log_ang_mom_tot']
sigma_low_redshift_all_log_ang_mom_tot_lower_error = sigma_low_redshift_all_dictionary['log_ang_mom_tot_lower_error']
sigma_low_redshift_all_log_ang_mom_tot_upper_error = sigma_low_redshift_all_dictionary['log_ang_mom_tot_upper_error']
sigma_low_redshift_all_log_ang_mom_tot_av_error = sigma_low_redshift_all_dictionary['log_ang_mom_tot_av_error']
sigma_low_redshift_all_log_velocity_mass_fit_params = sigma_low_redshift_all_dictionary['log_velocity_mass_fit_params']
print 'sigma_low_redshift_all velocity vs. mass fit params: %s' % sigma_low_redshift_all_log_velocity_mass_fit_params
sigma_low_redshift_all_log_velocity_mass_fit_line = sigma_low_redshift_all_dictionary['log_velocity_mass_fit_line']
sigma_low_redshift_all_log_velocity_mass_fit_line_lower = sigma_low_redshift_all_dictionary['log_velocity_mass_fit_line_lower']
sigma_low_redshift_all_log_velocity_mass_fit_line_upper = sigma_low_redshift_all_dictionary['log_velocity_mass_fit_line_upper']
sigma_low_redshift_all_log_v_tot_mass_fit_params = sigma_low_redshift_all_dictionary['log_v_tot_mass_fit_params']
print 'sigma_low_redshift_all v_tot vs. mass fit params: %s' % sigma_low_redshift_all_log_v_tot_mass_fit_params
sigma_low_redshift_all_log_v_tot_mass_fit_line = sigma_low_redshift_all_dictionary['log_v_tot_mass_fit_line']
sigma_low_redshift_all_log_v_tot_mass_fit_line_lower = sigma_low_redshift_all_dictionary['log_v_tot_mass_fit_line_lower']
sigma_low_redshift_all_log_v_tot_mass_fit_line_upper = sigma_low_redshift_all_dictionary['log_v_tot_mass_fit_line_upper']
sigma_low_redshift_all_log_ang_mom_mass_fit_params = sigma_low_redshift_all_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_low_redshift_all ang_mom vs. mass fit params: %s' % sigma_low_redshift_all_log_ang_mom_mass_fit_params
sigma_low_redshift_all_log_ang_mom_mass_fit_line = sigma_low_redshift_all_dictionary['log_ang_mom_mass_fit_line']
sigma_low_redshift_all_log_ang_mom_mass_fit_line_lower = sigma_low_redshift_all_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_low_redshift_all_log_ang_mom_mass_fit_line_upper = sigma_low_redshift_all_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params = sigma_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_low_redshift_all ang_mom_tot vs. mass fit params: %s' % sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params
sigma_low_redshift_all_log_ang_mom_tot_mass_fit_line = sigma_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_low_redshift_all_log_ang_mom_tot_mass_fit_line_lower = sigma_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_low_redshift_all_log_ang_mom_tot_mass_fit_line_upper = sigma_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_low_redshift_all_median_dynamical_mass_ratio = sigma_low_redshift_all_dictionary['median_dynamical_mass_ratio']
sigma_low_redshift_all_median_dynamical_mass_ratio_lower_error = sigma_low_redshift_all_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_low_redshift_all_median_dynamical_mass_ratio_upper_error = sigma_low_redshift_all_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio = sigma_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']


# DYNAMICAL PROPERTIES FOR SIGMA LOW REDSHIFT ROT

sigma_low_redshift_rot_velocity = table_sigma_low_redshift_rot['Vrot']
sigma_low_redshift_rot_velocity_upper_error = table_sigma_low_redshift_rot['Vrot_err']
sigma_low_redshift_rot_velocity_lower_error = table_sigma_low_redshift_rot['Vrot_err']
sigma_low_redshift_rot_mass = table_sigma_low_redshift_rot['Mstar']
sigma_low_redshift_rot_radius = table_sigma_low_redshift_rot['radius']
sigma_low_redshift_rot_sigma = table_sigma_low_redshift_rot['sigma']
sigma_low_redshift_rot_sigma_upper_error = table_sigma_low_redshift_rot['sigma_err']
sigma_low_redshift_rot_sigma_lower_error = table_sigma_low_redshift_rot['sigma_err']

sigma_low_redshift_rot_dictionary = return_betas(sigma_low_redshift_rot_mass,
                                    sigma_low_redshift_rot_velocity,
                                    sigma_low_redshift_rot_velocity_upper_error,
                                    sigma_low_redshift_rot_velocity_lower_error,
                                    sigma_low_redshift_rot_sigma,
                                    sigma_low_redshift_rot_sigma_upper_error,
                                    sigma_low_redshift_rot_sigma_lower_error,
                                    sigma_low_redshift_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_low_redshift_rot_dictionary['fit_mass']
sigma_low_redshift_rot_log_velocity = sigma_low_redshift_rot_dictionary['log_velocity']
sigma_low_redshift_rot_log_velocity_upper_error = sigma_low_redshift_rot_dictionary['log_velocity_upper_error']
sigma_low_redshift_rot_log_velocity_lower_error = sigma_low_redshift_rot_dictionary['log_velocity_lower_error']
sigma_low_redshift_rot_v_tot = sigma_low_redshift_rot_dictionary['v_tot']
sigma_low_redshift_rot_v_tot_lower_error = sigma_low_redshift_rot_dictionary['v_tot_lower_error']
sigma_low_redshift_rot_v_tot_upper_error = sigma_low_redshift_rot_dictionary['v_tot_upper_error']
sigma_low_redshift_rot_log_v_tot = sigma_low_redshift_rot_dictionary['log_v_tot']
sigma_low_redshift_rot_log_v_tot_lower_error = sigma_low_redshift_rot_dictionary['log_v_tot_lower_error']
sigma_low_redshift_rot_log_v_tot_upper_error = sigma_low_redshift_rot_dictionary['log_v_tot_upper_error']
sigma_low_redshift_rot_log_v_tot_av_error = sigma_low_redshift_rot_dictionary['log_v_tot_av_error']
sigma_low_redshift_rot_ang_mom = sigma_low_redshift_rot_dictionary['ang_mom']
sigma_low_redshift_rot_ang_mom_lower_error = sigma_low_redshift_rot_dictionary['ang_mom_lower_error']
sigma_low_redshift_rot_ang_mom_upper_error = sigma_low_redshift_rot_dictionary['ang_mom_upper_error']
sigma_low_redshift_rot_log_ang_mom = sigma_low_redshift_rot_dictionary['log_ang_mom']
sigma_low_redshift_rot_log_ang_mom_lower_error = sigma_low_redshift_rot_dictionary['log_ang_mom_lower_error']
sigma_low_redshift_rot_log_ang_mom_upper_error = sigma_low_redshift_rot_dictionary['log_ang_mom_upper_error']
sigma_low_redshift_rot_log_ang_mom_av_error = sigma_low_redshift_rot_dictionary['log_ang_mom_av_error']
sigma_low_redshift_rot_ang_mom_tot = sigma_low_redshift_rot_dictionary['ang_mom_tot']
sigma_low_redshift_rot_ang_mom_tot_lower_error = sigma_low_redshift_rot_dictionary['ang_mom_tot_lower_error']
sigma_low_redshift_rot_ang_mom_tot_upper_error = sigma_low_redshift_rot_dictionary['ang_mom_tot_upper_error']
sigma_low_redshift_rot_log_ang_mom_tot = sigma_low_redshift_rot_dictionary['log_ang_mom_tot']
sigma_low_redshift_rot_log_ang_mom_tot_lower_error = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_lower_error']
sigma_low_redshift_rot_log_ang_mom_tot_upper_error = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_upper_error']
sigma_low_redshift_rot_log_ang_mom_tot_av_error = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_av_error']
sigma_low_redshift_rot_log_velocity_mass_fit_params = sigma_low_redshift_rot_dictionary['log_velocity_mass_fit_params']
print 'sigma_low_redshift_rot velocity vs. mass fit params: %s' % sigma_low_redshift_rot_log_velocity_mass_fit_params
sigma_low_redshift_rot_log_velocity_mass_fit_line = sigma_low_redshift_rot_dictionary['log_velocity_mass_fit_line']
sigma_low_redshift_rot_log_velocity_mass_fit_line_lower = sigma_low_redshift_rot_dictionary['log_velocity_mass_fit_line_lower']
sigma_low_redshift_rot_log_velocity_mass_fit_line_upper = sigma_low_redshift_rot_dictionary['log_velocity_mass_fit_line_upper']
sigma_low_redshift_rot_log_v_tot_mass_fit_params = sigma_low_redshift_rot_dictionary['log_v_tot_mass_fit_params']
print 'sigma_low_redshift_rot v_tot vs. mass fit params: %s' % sigma_low_redshift_rot_log_v_tot_mass_fit_params
sigma_low_redshift_rot_log_v_tot_mass_fit_line = sigma_low_redshift_rot_dictionary['log_v_tot_mass_fit_line']
sigma_low_redshift_rot_log_v_tot_mass_fit_line_lower = sigma_low_redshift_rot_dictionary['log_v_tot_mass_fit_line_lower']
sigma_low_redshift_rot_log_v_tot_mass_fit_line_upper = sigma_low_redshift_rot_dictionary['log_v_tot_mass_fit_line_upper']
sigma_low_redshift_rot_log_ang_mom_mass_fit_params = sigma_low_redshift_rot_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_low_redshift_rot ang_mom vs. mass fit params: %s' % sigma_low_redshift_rot_log_ang_mom_mass_fit_params
sigma_low_redshift_rot_log_ang_mom_mass_fit_line = sigma_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line']
sigma_low_redshift_rot_log_ang_mom_mass_fit_line_lower = sigma_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_low_redshift_rot_log_ang_mom_mass_fit_line_upper = sigma_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_low_redshift_rot ang_mom_tot vs. mass fit params: %s' % sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params
sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_line = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_line_lower = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_line_upper = sigma_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_low_redshift_rot_median_dynamical_mass_ratio = sigma_low_redshift_rot_dictionary['median_dynamical_mass_ratio']
sigma_low_redshift_rot_median_dynamical_mass_ratio_lower_error = sigma_low_redshift_rot_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_low_redshift_rot_median_dynamical_mass_ratio_upper_error = sigma_low_redshift_rot_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio = sigma_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']


# DYNAMICAL PROPERTIES FOR SIGMA LOW REDSHIFT DISP

sigma_low_redshift_disp_velocity = table_sigma_low_redshift_disp['Vrot']
sigma_low_redshift_disp_velocity_upper_error = table_sigma_low_redshift_disp['Vrot_err']
sigma_low_redshift_disp_velocity_lower_error = table_sigma_low_redshift_disp['Vrot_err']
sigma_low_redshift_disp_mass = table_sigma_low_redshift_disp['Mstar']
sigma_low_redshift_disp_radius = table_sigma_low_redshift_disp['radius']
sigma_low_redshift_disp_sigma = table_sigma_low_redshift_disp['sigma']
sigma_low_redshift_disp_sigma_upper_error = table_sigma_low_redshift_disp['sigma_err']
sigma_low_redshift_disp_sigma_lower_error = table_sigma_low_redshift_disp['sigma_err']

sigma_low_redshift_disp_dictionary = return_betas(sigma_low_redshift_disp_mass,
                                    sigma_low_redshift_disp_velocity,
                                    sigma_low_redshift_disp_velocity_upper_error,
                                    sigma_low_redshift_disp_velocity_lower_error,
                                    sigma_low_redshift_disp_sigma,
                                    sigma_low_redshift_disp_sigma_upper_error,
                                    sigma_low_redshift_disp_sigma_lower_error,
                                    sigma_low_redshift_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_low_redshift_disp_dictionary['fit_mass']
sigma_low_redshift_disp_log_velocity = sigma_low_redshift_disp_dictionary['log_velocity']
sigma_low_redshift_disp_log_velocity_upper_error = sigma_low_redshift_disp_dictionary['log_velocity_upper_error']
sigma_low_redshift_disp_log_velocity_lower_error = sigma_low_redshift_disp_dictionary['log_velocity_lower_error']
sigma_low_redshift_disp_v_tot = sigma_low_redshift_disp_dictionary['v_tot']
sigma_low_redshift_disp_v_tot_lower_error = sigma_low_redshift_disp_dictionary['v_tot_lower_error']
sigma_low_redshift_disp_v_tot_upper_error = sigma_low_redshift_disp_dictionary['v_tot_upper_error']
sigma_low_redshift_disp_log_v_tot = sigma_low_redshift_disp_dictionary['log_v_tot']
sigma_low_redshift_disp_log_v_tot_lower_error = sigma_low_redshift_disp_dictionary['log_v_tot_lower_error']
sigma_low_redshift_disp_log_v_tot_upper_error = sigma_low_redshift_disp_dictionary['log_v_tot_upper_error']
sigma_low_redshift_disp_log_v_tot_av_error = sigma_low_redshift_disp_dictionary['log_v_tot_av_error']
sigma_low_redshift_disp_ang_mom = sigma_low_redshift_disp_dictionary['ang_mom']
sigma_low_redshift_disp_ang_mom_lower_error = sigma_low_redshift_disp_dictionary['ang_mom_lower_error']
sigma_low_redshift_disp_ang_mom_upper_error = sigma_low_redshift_disp_dictionary['ang_mom_upper_error']
sigma_low_redshift_disp_log_ang_mom = sigma_low_redshift_disp_dictionary['log_ang_mom']
sigma_low_redshift_disp_log_ang_mom_lower_error = sigma_low_redshift_disp_dictionary['log_ang_mom_lower_error']
sigma_low_redshift_disp_log_ang_mom_upper_error = sigma_low_redshift_disp_dictionary['log_ang_mom_upper_error']
sigma_low_redshift_disp_log_ang_mom_av_error = sigma_low_redshift_disp_dictionary['log_ang_mom_av_error']
sigma_low_redshift_disp_ang_mom_tot = sigma_low_redshift_disp_dictionary['ang_mom_tot']
sigma_low_redshift_disp_ang_mom_tot_lower_error = sigma_low_redshift_disp_dictionary['ang_mom_tot_lower_error']
sigma_low_redshift_disp_ang_mom_tot_upper_error = sigma_low_redshift_disp_dictionary['ang_mom_tot_upper_error']
sigma_low_redshift_disp_log_ang_mom_tot = sigma_low_redshift_disp_dictionary['log_ang_mom_tot']
sigma_low_redshift_disp_log_ang_mom_tot_lower_error = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_lower_error']
sigma_low_redshift_disp_log_ang_mom_tot_upper_error = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_upper_error']
sigma_low_redshift_disp_log_ang_mom_tot_av_error = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_av_error']
sigma_low_redshift_disp_log_velocity_mass_fit_params = sigma_low_redshift_disp_dictionary['log_velocity_mass_fit_params']
print 'sigma_low_redshift_disp velocity vs. mass fit params: %s' % sigma_low_redshift_disp_log_velocity_mass_fit_params
sigma_low_redshift_disp_log_velocity_mass_fit_line = sigma_low_redshift_disp_dictionary['log_velocity_mass_fit_line']
sigma_low_redshift_disp_log_velocity_mass_fit_line_lower = sigma_low_redshift_disp_dictionary['log_velocity_mass_fit_line_lower']
sigma_low_redshift_disp_log_velocity_mass_fit_line_upper = sigma_low_redshift_disp_dictionary['log_velocity_mass_fit_line_upper']
sigma_low_redshift_disp_log_v_tot_mass_fit_params = sigma_low_redshift_disp_dictionary['log_v_tot_mass_fit_params']
print 'sigma_low_redshift_disp v_tot vs. mass fit params: %s' % sigma_low_redshift_disp_log_v_tot_mass_fit_params
sigma_low_redshift_disp_log_v_tot_mass_fit_line = sigma_low_redshift_disp_dictionary['log_v_tot_mass_fit_line']
sigma_low_redshift_disp_log_v_tot_mass_fit_line_lower = sigma_low_redshift_disp_dictionary['log_v_tot_mass_fit_line_lower']
sigma_low_redshift_disp_log_v_tot_mass_fit_line_upper = sigma_low_redshift_disp_dictionary['log_v_tot_mass_fit_line_upper']
sigma_low_redshift_disp_log_ang_mom_mass_fit_params = sigma_low_redshift_disp_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_low_redshift_disp ang_mom vs. mass fit params: %s' % sigma_low_redshift_disp_log_ang_mom_mass_fit_params
sigma_low_redshift_disp_log_ang_mom_mass_fit_line = sigma_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line']
sigma_low_redshift_disp_log_ang_mom_mass_fit_line_lower = sigma_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_low_redshift_disp_log_ang_mom_mass_fit_line_upper = sigma_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_low_redshift_disp ang_mom_tot vs. mass fit params: %s' % sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params
sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_line = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_line_lower = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_line_upper = sigma_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_low_redshift_disp_median_dynamical_mass_ratio = sigma_low_redshift_disp_dictionary['median_dynamical_mass_ratio']
sigma_low_redshift_disp_median_dynamical_mass_ratio_lower_error = sigma_low_redshift_disp_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_low_redshift_disp_median_dynamical_mass_ratio_upper_error = sigma_low_redshift_disp_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio = sigma_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['sigma_low_redshift',
                                1.5,
                                sigma_low_redshift_all_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_all_log_velocity_mass_fit_params[3],
                                sigma_low_redshift_all_log_velocity_mass_fit_params[0],
                                sigma_low_redshift_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_all_log_velocity_mass_fit_params[0] - sigma_low_redshift_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_all_log_velocity_mass_fit_params[2] - sigma_low_redshift_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                sigma_low_redshift_rot_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_rot_log_velocity_mass_fit_params[3],
                                sigma_low_redshift_rot_log_velocity_mass_fit_params[0],
                                sigma_low_redshift_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_rot_log_velocity_mass_fit_params[0] - sigma_low_redshift_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_rot_log_velocity_mass_fit_params[2] - sigma_low_redshift_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                sigma_low_redshift_disp_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_disp_log_velocity_mass_fit_params[3],
                                sigma_low_redshift_disp_log_velocity_mass_fit_params[0],
                                sigma_low_redshift_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_disp_log_velocity_mass_fit_params[0] - sigma_low_redshift_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_disp_log_velocity_mass_fit_params[2] - sigma_low_redshift_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['sigma_low_redshift',
                                1.5,
                                sigma_low_redshift_all_log_v_tot_mass_fit_params[4],
                                sigma_low_redshift_all_log_v_tot_mass_fit_params[3],
                                sigma_low_redshift_all_log_v_tot_mass_fit_params[0],
                                sigma_low_redshift_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_all_log_v_tot_mass_fit_params[0] - sigma_low_redshift_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_all_log_v_tot_mass_fit_params[2] - sigma_low_redshift_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                sigma_low_redshift_rot_log_v_tot_mass_fit_params[4],
                                sigma_low_redshift_rot_log_v_tot_mass_fit_params[3],
                                sigma_low_redshift_rot_log_v_tot_mass_fit_params[0],
                                sigma_low_redshift_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_rot_log_v_tot_mass_fit_params[0] - sigma_low_redshift_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_rot_log_v_tot_mass_fit_params[2] - sigma_low_redshift_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                sigma_low_redshift_disp_log_v_tot_mass_fit_params[4],
                                sigma_low_redshift_disp_log_v_tot_mass_fit_params[3],
                                sigma_low_redshift_disp_log_v_tot_mass_fit_params[0],
                                sigma_low_redshift_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_low_redshift_disp_log_v_tot_mass_fit_params[0] - sigma_low_redshift_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_disp_log_v_tot_mass_fit_params[2] - sigma_low_redshift_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['sigma_low_redshift',
                                1.5,
                                sigma_low_redshift_all_log_ang_mom_mass_fit_params[4],
                                sigma_low_redshift_all_log_ang_mom_mass_fit_params[3],
                                sigma_low_redshift_all_log_ang_mom_mass_fit_params[0],
                                sigma_low_redshift_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_all_log_ang_mom_mass_fit_params[0] - sigma_low_redshift_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_all_log_ang_mom_mass_fit_params[2] - sigma_low_redshift_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                sigma_low_redshift_rot_log_ang_mom_mass_fit_params[4],
                                sigma_low_redshift_rot_log_ang_mom_mass_fit_params[3],
                                sigma_low_redshift_rot_log_ang_mom_mass_fit_params[0],
                                sigma_low_redshift_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_rot_log_ang_mom_mass_fit_params[0] - sigma_low_redshift_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_rot_log_ang_mom_mass_fit_params[2] - sigma_low_redshift_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                sigma_low_redshift_disp_log_ang_mom_mass_fit_params[4],
                                sigma_low_redshift_disp_log_ang_mom_mass_fit_params[3],
                                sigma_low_redshift_disp_log_ang_mom_mass_fit_params[0],
                                sigma_low_redshift_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_disp_log_ang_mom_mass_fit_params[0] - sigma_low_redshift_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_disp_log_ang_mom_mass_fit_params[2] - sigma_low_redshift_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['sigma_low_redshift',
                                1.5,
                                sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[4],
                                sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[3],
                                sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[0],
                                sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[0] - sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[2] - sigma_low_redshift_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[4],
                                sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[3],
                                sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0],
                                sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[2] - sigma_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[4],
                                sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[3],
                                sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0],
                                sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[2] - sigma_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['sigma_low_redshift',
                                1.5,
                                sigma_low_redshift_all_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_all_median_dynamical_mass_ratio,
                                sigma_low_redshift_all_median_dynamical_mass_ratio_lower_error,
                                sigma_low_redshift_all_median_dynamical_mass_ratio_upper_error,
                                sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio,
                                sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_low_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                sigma_low_redshift_rot_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_rot_median_dynamical_mass_ratio,
                                sigma_low_redshift_rot_median_dynamical_mass_ratio_lower_error,
                                sigma_low_redshift_rot_median_dynamical_mass_ratio_upper_error,
                                sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio,
                                sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                sigma_low_redshift_disp_log_velocity_mass_fit_params[4],
                                sigma_low_redshift_disp_median_dynamical_mass_ratio,
                                sigma_low_redshift_disp_median_dynamical_mass_ratio_lower_error,
                                sigma_low_redshift_disp_median_dynamical_mass_ratio_upper_error,
                                sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio,
                                sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error])

# DYNAMICAL PROPERTIES FOR SIGMA HIGH REDSHIFT ALL

sigma_high_redshift_all_velocity = table_sigma_high_redshift_all['Vrot']
sigma_high_redshift_all_velocity_upper_error = table_sigma_high_redshift_all['Vrot_err']
sigma_high_redshift_all_velocity_lower_error = table_sigma_high_redshift_all['Vrot_err']
sigma_high_redshift_all_mass = table_sigma_high_redshift_all['Mstar']
sigma_high_redshift_all_radius = table_sigma_high_redshift_all['radius']
sigma_high_redshift_all_sigma = table_sigma_high_redshift_all['sigma']
sigma_high_redshift_all_sigma_upper_error = table_sigma_high_redshift_all['sigma_err']
sigma_high_redshift_all_sigma_lower_error = table_sigma_high_redshift_all['sigma_err']

sigma_high_redshift_all_dictionary = return_betas(sigma_high_redshift_all_mass,
                                    sigma_high_redshift_all_velocity,
                                    sigma_high_redshift_all_velocity_upper_error,
                                    sigma_high_redshift_all_velocity_lower_error,
                                    sigma_high_redshift_all_sigma,
                                    sigma_high_redshift_all_sigma_upper_error,
                                    sigma_high_redshift_all_sigma_lower_error,
                                    sigma_high_redshift_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_high_redshift_all_dictionary['fit_mass']
sigma_high_redshift_all_log_velocity = sigma_high_redshift_all_dictionary['log_velocity']
sigma_high_redshift_all_log_velocity_upper_error = sigma_high_redshift_all_dictionary['log_velocity_upper_error']
sigma_high_redshift_all_log_velocity_lower_error = sigma_high_redshift_all_dictionary['log_velocity_lower_error']
sigma_high_redshift_all_v_tot = sigma_high_redshift_all_dictionary['v_tot']
sigma_high_redshift_all_v_tot_lower_error = sigma_high_redshift_all_dictionary['v_tot_lower_error']
sigma_high_redshift_all_v_tot_upper_error = sigma_high_redshift_all_dictionary['v_tot_upper_error']
sigma_high_redshift_all_log_v_tot = sigma_high_redshift_all_dictionary['log_v_tot']
sigma_high_redshift_all_log_v_tot_lower_error = sigma_high_redshift_all_dictionary['log_v_tot_lower_error']
sigma_high_redshift_all_log_v_tot_upper_error = sigma_high_redshift_all_dictionary['log_v_tot_upper_error']
sigma_high_redshift_all_log_v_tot_av_error = sigma_high_redshift_all_dictionary['log_v_tot_av_error']
sigma_high_redshift_all_ang_mom = sigma_high_redshift_all_dictionary['ang_mom']
sigma_high_redshift_all_ang_mom_lower_error = sigma_high_redshift_all_dictionary['ang_mom_lower_error']
sigma_high_redshift_all_ang_mom_upper_error = sigma_high_redshift_all_dictionary['ang_mom_upper_error']
sigma_high_redshift_all_log_ang_mom = sigma_high_redshift_all_dictionary['log_ang_mom']
sigma_high_redshift_all_log_ang_mom_lower_error = sigma_high_redshift_all_dictionary['log_ang_mom_lower_error']
sigma_high_redshift_all_log_ang_mom_upper_error = sigma_high_redshift_all_dictionary['log_ang_mom_upper_error']
sigma_high_redshift_all_log_ang_mom_av_error = sigma_high_redshift_all_dictionary['log_ang_mom_av_error']
sigma_high_redshift_all_ang_mom_tot = sigma_high_redshift_all_dictionary['ang_mom_tot']
sigma_high_redshift_all_ang_mom_tot_lower_error = sigma_high_redshift_all_dictionary['ang_mom_tot_lower_error']
sigma_high_redshift_all_ang_mom_tot_upper_error = sigma_high_redshift_all_dictionary['ang_mom_tot_upper_error']
sigma_high_redshift_all_log_ang_mom_tot = sigma_high_redshift_all_dictionary['log_ang_mom_tot']
sigma_high_redshift_all_log_ang_mom_tot_lower_error = sigma_high_redshift_all_dictionary['log_ang_mom_tot_lower_error']
sigma_high_redshift_all_log_ang_mom_tot_upper_error = sigma_high_redshift_all_dictionary['log_ang_mom_tot_upper_error']
sigma_high_redshift_all_log_ang_mom_tot_av_error = sigma_high_redshift_all_dictionary['log_ang_mom_tot_av_error']
sigma_high_redshift_all_log_velocity_mass_fit_params = sigma_high_redshift_all_dictionary['log_velocity_mass_fit_params']
print 'sigma_high_redshift_all velocity vs. mass fit params: %s' % sigma_high_redshift_all_log_velocity_mass_fit_params
sigma_high_redshift_all_log_velocity_mass_fit_line = sigma_high_redshift_all_dictionary['log_velocity_mass_fit_line']
sigma_high_redshift_all_log_velocity_mass_fit_line_lower = sigma_high_redshift_all_dictionary['log_velocity_mass_fit_line_lower']
sigma_high_redshift_all_log_velocity_mass_fit_line_upper = sigma_high_redshift_all_dictionary['log_velocity_mass_fit_line_upper']
sigma_high_redshift_all_log_v_tot_mass_fit_params = sigma_high_redshift_all_dictionary['log_v_tot_mass_fit_params']
print 'sigma_high_redshift_all v_tot vs. mass fit params: %s' % sigma_high_redshift_all_log_v_tot_mass_fit_params
sigma_high_redshift_all_log_v_tot_mass_fit_line = sigma_high_redshift_all_dictionary['log_v_tot_mass_fit_line']
sigma_high_redshift_all_log_v_tot_mass_fit_line_lower = sigma_high_redshift_all_dictionary['log_v_tot_mass_fit_line_lower']
sigma_high_redshift_all_log_v_tot_mass_fit_line_upper = sigma_high_redshift_all_dictionary['log_v_tot_mass_fit_line_upper']
sigma_high_redshift_all_log_ang_mom_mass_fit_params = sigma_high_redshift_all_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_high_redshift_all ang_mom vs. mass fit params: %s' % sigma_high_redshift_all_log_ang_mom_mass_fit_params
sigma_high_redshift_all_log_ang_mom_mass_fit_line = sigma_high_redshift_all_dictionary['log_ang_mom_mass_fit_line']
sigma_high_redshift_all_log_ang_mom_mass_fit_line_lower = sigma_high_redshift_all_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_high_redshift_all_log_ang_mom_mass_fit_line_upper = sigma_high_redshift_all_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params = sigma_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_high_redshift_all ang_mom_tot vs. mass fit params: %s' % sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params
sigma_high_redshift_all_log_ang_mom_tot_mass_fit_line = sigma_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_high_redshift_all_log_ang_mom_tot_mass_fit_line_lower = sigma_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_high_redshift_all_log_ang_mom_tot_mass_fit_line_upper = sigma_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_high_redshift_all_median_dynamical_mass_ratio = sigma_high_redshift_all_dictionary['median_dynamical_mass_ratio']
sigma_high_redshift_all_median_dynamical_mass_ratio_lower_error = sigma_high_redshift_all_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_high_redshift_all_median_dynamical_mass_ratio_upper_error = sigma_high_redshift_all_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio = sigma_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SIGMA HIGH REDSHIFT ROT

sigma_high_redshift_rot_velocity = table_sigma_high_redshift_rot['Vrot']
sigma_high_redshift_rot_velocity_upper_error = table_sigma_high_redshift_rot['Vrot_err']
sigma_high_redshift_rot_velocity_lower_error = table_sigma_high_redshift_rot['Vrot_err']
sigma_high_redshift_rot_mass = table_sigma_high_redshift_rot['Mstar']
sigma_high_redshift_rot_radius = table_sigma_high_redshift_rot['radius']
sigma_high_redshift_rot_sigma = table_sigma_high_redshift_rot['sigma']
sigma_high_redshift_rot_sigma_upper_error = table_sigma_high_redshift_rot['sigma_err']
sigma_high_redshift_rot_sigma_lower_error = table_sigma_high_redshift_rot['sigma_err']

sigma_high_redshift_rot_dictionary = return_betas(sigma_high_redshift_rot_mass,
                                    sigma_high_redshift_rot_velocity,
                                    sigma_high_redshift_rot_velocity_upper_error,
                                    sigma_high_redshift_rot_velocity_lower_error,
                                    sigma_high_redshift_rot_sigma,
                                    sigma_high_redshift_rot_sigma_upper_error,
                                    sigma_high_redshift_rot_sigma_lower_error,
                                    sigma_high_redshift_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_high_redshift_rot_dictionary['fit_mass']
sigma_high_redshift_rot_log_velocity = sigma_high_redshift_rot_dictionary['log_velocity']
sigma_high_redshift_rot_log_velocity_upper_error = sigma_high_redshift_rot_dictionary['log_velocity_upper_error']
sigma_high_redshift_rot_log_velocity_lower_error = sigma_high_redshift_rot_dictionary['log_velocity_lower_error']
sigma_high_redshift_rot_v_tot = sigma_high_redshift_rot_dictionary['v_tot']
sigma_high_redshift_rot_v_tot_lower_error = sigma_high_redshift_rot_dictionary['v_tot_lower_error']
sigma_high_redshift_rot_v_tot_upper_error = sigma_high_redshift_rot_dictionary['v_tot_upper_error']
sigma_high_redshift_rot_log_v_tot = sigma_high_redshift_rot_dictionary['log_v_tot']
sigma_high_redshift_rot_log_v_tot_lower_error = sigma_high_redshift_rot_dictionary['log_v_tot_lower_error']
sigma_high_redshift_rot_log_v_tot_upper_error = sigma_high_redshift_rot_dictionary['log_v_tot_upper_error']
sigma_high_redshift_rot_log_v_tot_av_error = sigma_high_redshift_rot_dictionary['log_v_tot_av_error']
sigma_high_redshift_rot_ang_mom = sigma_high_redshift_rot_dictionary['ang_mom']
sigma_high_redshift_rot_ang_mom_lower_error = sigma_high_redshift_rot_dictionary['ang_mom_lower_error']
sigma_high_redshift_rot_ang_mom_upper_error = sigma_high_redshift_rot_dictionary['ang_mom_upper_error']
sigma_high_redshift_rot_log_ang_mom = sigma_high_redshift_rot_dictionary['log_ang_mom']
sigma_high_redshift_rot_log_ang_mom_lower_error = sigma_high_redshift_rot_dictionary['log_ang_mom_lower_error']
sigma_high_redshift_rot_log_ang_mom_upper_error = sigma_high_redshift_rot_dictionary['log_ang_mom_upper_error']
sigma_high_redshift_rot_log_ang_mom_av_error = sigma_high_redshift_rot_dictionary['log_ang_mom_av_error']
sigma_high_redshift_rot_ang_mom_tot = sigma_high_redshift_rot_dictionary['ang_mom_tot']
sigma_high_redshift_rot_ang_mom_tot_lower_error = sigma_high_redshift_rot_dictionary['ang_mom_tot_lower_error']
sigma_high_redshift_rot_ang_mom_tot_upper_error = sigma_high_redshift_rot_dictionary['ang_mom_tot_upper_error']
sigma_high_redshift_rot_log_ang_mom_tot = sigma_high_redshift_rot_dictionary['log_ang_mom_tot']
sigma_high_redshift_rot_log_ang_mom_tot_lower_error = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_lower_error']
sigma_high_redshift_rot_log_ang_mom_tot_upper_error = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_upper_error']
sigma_high_redshift_rot_log_ang_mom_tot_av_error = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_av_error']
sigma_high_redshift_rot_log_velocity_mass_fit_params = sigma_high_redshift_rot_dictionary['log_velocity_mass_fit_params']
print 'sigma_high_redshift_rot velocity vs. mass fit params: %s' % sigma_high_redshift_rot_log_velocity_mass_fit_params
sigma_high_redshift_rot_log_velocity_mass_fit_line = sigma_high_redshift_rot_dictionary['log_velocity_mass_fit_line']
sigma_high_redshift_rot_log_velocity_mass_fit_line_lower = sigma_high_redshift_rot_dictionary['log_velocity_mass_fit_line_lower']
sigma_high_redshift_rot_log_velocity_mass_fit_line_upper = sigma_high_redshift_rot_dictionary['log_velocity_mass_fit_line_upper']
sigma_high_redshift_rot_log_v_tot_mass_fit_params = sigma_high_redshift_rot_dictionary['log_v_tot_mass_fit_params']
print 'sigma_high_redshift_rot v_tot vs. mass fit params: %s' % sigma_high_redshift_rot_log_v_tot_mass_fit_params
sigma_high_redshift_rot_log_v_tot_mass_fit_line = sigma_high_redshift_rot_dictionary['log_v_tot_mass_fit_line']
sigma_high_redshift_rot_log_v_tot_mass_fit_line_lower = sigma_high_redshift_rot_dictionary['log_v_tot_mass_fit_line_lower']
sigma_high_redshift_rot_log_v_tot_mass_fit_line_upper = sigma_high_redshift_rot_dictionary['log_v_tot_mass_fit_line_upper']
sigma_high_redshift_rot_log_ang_mom_mass_fit_params = sigma_high_redshift_rot_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_high_redshift_rot ang_mom vs. mass fit params: %s' % sigma_high_redshift_rot_log_ang_mom_mass_fit_params
sigma_high_redshift_rot_log_ang_mom_mass_fit_line = sigma_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line']
sigma_high_redshift_rot_log_ang_mom_mass_fit_line_lower = sigma_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_high_redshift_rot_log_ang_mom_mass_fit_line_upper = sigma_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_high_redshift_rot ang_mom_tot vs. mass fit params: %s' % sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params
sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_line = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_line_lower = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_line_upper = sigma_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_high_redshift_rot_median_dynamical_mass_ratio = sigma_high_redshift_rot_dictionary['median_dynamical_mass_ratio']
sigma_high_redshift_rot_median_dynamical_mass_ratio_lower_error = sigma_high_redshift_rot_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_high_redshift_rot_median_dynamical_mass_ratio_upper_error = sigma_high_redshift_rot_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio = sigma_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SIGMA HIGH REDSHIFT DISP

sigma_high_redshift_disp_velocity = table_sigma_high_redshift_disp['Vrot']
sigma_high_redshift_disp_velocity_upper_error = table_sigma_high_redshift_disp['Vrot_err']
sigma_high_redshift_disp_velocity_lower_error = table_sigma_high_redshift_disp['Vrot_err']
sigma_high_redshift_disp_mass = table_sigma_high_redshift_disp['Mstar']
sigma_high_redshift_disp_radius = table_sigma_high_redshift_disp['radius']
sigma_high_redshift_disp_sigma = table_sigma_high_redshift_disp['sigma']
sigma_high_redshift_disp_sigma_upper_error = table_sigma_high_redshift_disp['sigma_err']
sigma_high_redshift_disp_sigma_lower_error = table_sigma_high_redshift_disp['sigma_err']

sigma_high_redshift_disp_dictionary = return_betas(sigma_high_redshift_disp_mass,
                                    sigma_high_redshift_disp_velocity,
                                    sigma_high_redshift_disp_velocity_upper_error,
                                    sigma_high_redshift_disp_velocity_lower_error,
                                    sigma_high_redshift_disp_sigma,
                                    sigma_high_redshift_disp_sigma_upper_error,
                                    sigma_high_redshift_disp_sigma_lower_error,
                                    sigma_high_redshift_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = sigma_high_redshift_disp_dictionary['fit_mass']
sigma_high_redshift_disp_log_velocity = sigma_high_redshift_disp_dictionary['log_velocity']
sigma_high_redshift_disp_log_velocity_upper_error = sigma_high_redshift_disp_dictionary['log_velocity_upper_error']
sigma_high_redshift_disp_log_velocity_lower_error = sigma_high_redshift_disp_dictionary['log_velocity_lower_error']
sigma_high_redshift_disp_v_tot = sigma_high_redshift_disp_dictionary['v_tot']
sigma_high_redshift_disp_v_tot_lower_error = sigma_high_redshift_disp_dictionary['v_tot_lower_error']
sigma_high_redshift_disp_v_tot_upper_error = sigma_high_redshift_disp_dictionary['v_tot_upper_error']
sigma_high_redshift_disp_log_v_tot = sigma_high_redshift_disp_dictionary['log_v_tot']
sigma_high_redshift_disp_log_v_tot_lower_error = sigma_high_redshift_disp_dictionary['log_v_tot_lower_error']
sigma_high_redshift_disp_log_v_tot_upper_error = sigma_high_redshift_disp_dictionary['log_v_tot_upper_error']
sigma_high_redshift_disp_log_v_tot_av_error = sigma_high_redshift_disp_dictionary['log_v_tot_av_error']
sigma_high_redshift_disp_ang_mom = sigma_high_redshift_disp_dictionary['ang_mom']
sigma_high_redshift_disp_ang_mom_lower_error = sigma_high_redshift_disp_dictionary['ang_mom_lower_error']
sigma_high_redshift_disp_ang_mom_upper_error = sigma_high_redshift_disp_dictionary['ang_mom_upper_error']
sigma_high_redshift_disp_log_ang_mom = sigma_high_redshift_disp_dictionary['log_ang_mom']
sigma_high_redshift_disp_log_ang_mom_lower_error = sigma_high_redshift_disp_dictionary['log_ang_mom_lower_error']
sigma_high_redshift_disp_log_ang_mom_upper_error = sigma_high_redshift_disp_dictionary['log_ang_mom_upper_error']
sigma_high_redshift_disp_log_ang_mom_av_error = sigma_high_redshift_disp_dictionary['log_ang_mom_av_error']
sigma_high_redshift_disp_ang_mom_tot = sigma_high_redshift_disp_dictionary['ang_mom_tot']
sigma_high_redshift_disp_ang_mom_tot_lower_error = sigma_high_redshift_disp_dictionary['ang_mom_tot_lower_error']
sigma_high_redshift_disp_ang_mom_tot_upper_error = sigma_high_redshift_disp_dictionary['ang_mom_tot_upper_error']
sigma_high_redshift_disp_log_ang_mom_tot = sigma_high_redshift_disp_dictionary['log_ang_mom_tot']
sigma_high_redshift_disp_log_ang_mom_tot_lower_error = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_lower_error']
sigma_high_redshift_disp_log_ang_mom_tot_upper_error = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_upper_error']
sigma_high_redshift_disp_log_ang_mom_tot_av_error = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_av_error']
sigma_high_redshift_disp_log_velocity_mass_fit_params = sigma_high_redshift_disp_dictionary['log_velocity_mass_fit_params']
print 'sigma_high_redshift_disp velocity vs. mass fit params: %s' % sigma_high_redshift_disp_log_velocity_mass_fit_params
sigma_high_redshift_disp_log_velocity_mass_fit_line = sigma_high_redshift_disp_dictionary['log_velocity_mass_fit_line']
sigma_high_redshift_disp_log_velocity_mass_fit_line_lower = sigma_high_redshift_disp_dictionary['log_velocity_mass_fit_line_lower']
sigma_high_redshift_disp_log_velocity_mass_fit_line_upper = sigma_high_redshift_disp_dictionary['log_velocity_mass_fit_line_upper']
sigma_high_redshift_disp_log_v_tot_mass_fit_params = sigma_high_redshift_disp_dictionary['log_v_tot_mass_fit_params']
print 'sigma_high_redshift_disp v_tot vs. mass fit params: %s' % sigma_high_redshift_disp_log_v_tot_mass_fit_params
sigma_high_redshift_disp_log_v_tot_mass_fit_line = sigma_high_redshift_disp_dictionary['log_v_tot_mass_fit_line']
sigma_high_redshift_disp_log_v_tot_mass_fit_line_lower = sigma_high_redshift_disp_dictionary['log_v_tot_mass_fit_line_lower']
sigma_high_redshift_disp_log_v_tot_mass_fit_line_upper = sigma_high_redshift_disp_dictionary['log_v_tot_mass_fit_line_upper']
sigma_high_redshift_disp_log_ang_mom_mass_fit_params = sigma_high_redshift_disp_dictionary['log_ang_mom_mass_fit_params']
print 'sigma_high_redshift_disp ang_mom vs. mass fit params: %s' % sigma_high_redshift_disp_log_ang_mom_mass_fit_params
sigma_high_redshift_disp_log_ang_mom_mass_fit_line = sigma_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line']
sigma_high_redshift_disp_log_ang_mom_mass_fit_line_lower = sigma_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line_lower']
sigma_high_redshift_disp_log_ang_mom_mass_fit_line_upper = sigma_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line_upper']
sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'sigma_high_redshift_disp ang_mom_tot vs. mass fit params: %s' % sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params
sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_line = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line']
sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_line_lower = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_line_upper = sigma_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
sigma_high_redshift_disp_median_dynamical_mass_ratio = sigma_high_redshift_disp_dictionary['median_dynamical_mass_ratio']
sigma_high_redshift_disp_median_dynamical_mass_ratio_lower_error = sigma_high_redshift_disp_dictionary['median_dynamical_mass_ratio_lower_error']
sigma_high_redshift_disp_median_dynamical_mass_ratio_upper_error = sigma_high_redshift_disp_dictionary['median_dynamical_mass_ratio_upper_error']
sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio = sigma_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error = sigma_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error = sigma_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['sigma_high_redshift',
                                2.25,
                                sigma_high_redshift_all_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_all_log_velocity_mass_fit_params[3],
                                sigma_high_redshift_all_log_velocity_mass_fit_params[0],
                                sigma_high_redshift_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_all_log_velocity_mass_fit_params[0] - sigma_high_redshift_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_all_log_velocity_mass_fit_params[2] - sigma_high_redshift_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                sigma_high_redshift_rot_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_rot_log_velocity_mass_fit_params[3],
                                sigma_high_redshift_rot_log_velocity_mass_fit_params[0],
                                sigma_high_redshift_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_rot_log_velocity_mass_fit_params[0] - sigma_high_redshift_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_rot_log_velocity_mass_fit_params[2] - sigma_high_redshift_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                sigma_high_redshift_disp_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_disp_log_velocity_mass_fit_params[3],
                                sigma_high_redshift_disp_log_velocity_mass_fit_params[0],
                                sigma_high_redshift_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_disp_log_velocity_mass_fit_params[0] - sigma_high_redshift_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_disp_log_velocity_mass_fit_params[2] - sigma_high_redshift_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2)])
v_tot_evolution_list.append(['sigma_high_redshift',
                                2.25,
                                sigma_high_redshift_all_log_v_tot_mass_fit_params[4],
                                sigma_high_redshift_all_log_v_tot_mass_fit_params[3],
                                sigma_high_redshift_all_log_v_tot_mass_fit_params[0],
                                sigma_high_redshift_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_all_log_v_tot_mass_fit_params[0] - sigma_high_redshift_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_all_log_v_tot_mass_fit_params[2] - sigma_high_redshift_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                sigma_high_redshift_rot_log_v_tot_mass_fit_params[4],
                                sigma_high_redshift_rot_log_v_tot_mass_fit_params[3],
                                sigma_high_redshift_rot_log_v_tot_mass_fit_params[0],
                                sigma_high_redshift_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_rot_log_v_tot_mass_fit_params[0] - sigma_high_redshift_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_rot_log_v_tot_mass_fit_params[2] - sigma_high_redshift_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                sigma_high_redshift_disp_log_v_tot_mass_fit_params[4],
                                sigma_high_redshift_disp_log_v_tot_mass_fit_params[3],
                                sigma_high_redshift_disp_log_v_tot_mass_fit_params[0],
                                sigma_high_redshift_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((sigma_high_redshift_disp_log_v_tot_mass_fit_params[0] - sigma_high_redshift_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_disp_log_v_tot_mass_fit_params[2] - sigma_high_redshift_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['sigma_high_redshift',
                                2.25,
                                sigma_high_redshift_all_log_ang_mom_mass_fit_params[4],
                                sigma_high_redshift_all_log_ang_mom_mass_fit_params[3],
                                sigma_high_redshift_all_log_ang_mom_mass_fit_params[0],
                                sigma_high_redshift_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_all_log_ang_mom_mass_fit_params[0] - sigma_high_redshift_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_all_log_ang_mom_mass_fit_params[2] - sigma_high_redshift_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                sigma_high_redshift_rot_log_ang_mom_mass_fit_params[4],
                                sigma_high_redshift_rot_log_ang_mom_mass_fit_params[3],
                                sigma_high_redshift_rot_log_ang_mom_mass_fit_params[0],
                                sigma_high_redshift_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_rot_log_ang_mom_mass_fit_params[0] - sigma_high_redshift_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_rot_log_ang_mom_mass_fit_params[2] - sigma_high_redshift_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                sigma_high_redshift_disp_log_ang_mom_mass_fit_params[4],
                                sigma_high_redshift_disp_log_ang_mom_mass_fit_params[3],
                                sigma_high_redshift_disp_log_ang_mom_mass_fit_params[0],
                                sigma_high_redshift_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_disp_log_ang_mom_mass_fit_params[0] - sigma_high_redshift_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_disp_log_ang_mom_mass_fit_params[2] - sigma_high_redshift_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['sigma_high_redshift',
                                2.25,
                                sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[4],
                                sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[3],
                                sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[0],
                                sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[0] - sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[2] - sigma_high_redshift_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[4],
                                sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[3],
                                sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0],
                                sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[2] - sigma_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[4],
                                sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[3],
                                sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0],
                                sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[2] - sigma_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['sigma_high_redshift',
                                2.25,
                                sigma_high_redshift_all_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_all_median_dynamical_mass_ratio,
                                sigma_high_redshift_all_median_dynamical_mass_ratio_lower_error,
                                sigma_high_redshift_all_median_dynamical_mass_ratio_upper_error,
                                sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio,
                                sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_high_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                sigma_high_redshift_rot_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_rot_median_dynamical_mass_ratio,
                                sigma_high_redshift_rot_median_dynamical_mass_ratio_lower_error,
                                sigma_high_redshift_rot_median_dynamical_mass_ratio_upper_error,
                                sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio,
                                sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                sigma_high_redshift_disp_log_velocity_mass_fit_params[4],
                                sigma_high_redshift_disp_median_dynamical_mass_ratio,
                                sigma_high_redshift_disp_median_dynamical_mass_ratio_lower_error,
                                sigma_high_redshift_disp_median_dynamical_mass_ratio_upper_error,
                                sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio,
                                sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                sigma_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error])


# DYNAMICAL PROPERTIES FOR SWINBANK LOW REDSHIFT ALL

swinbank_low_redshift_all_velocity = table_swinbank_low_redshift_all['V']
swinbank_low_redshift_all_velocity_upper_error = table_swinbank_low_redshift_all['V_err']
swinbank_low_redshift_all_velocity_lower_error = table_swinbank_low_redshift_all['V_err']
swinbank_low_redshift_all_mass = table_swinbank_low_redshift_all['log(M*)']
swinbank_low_redshift_all_radius = table_swinbank_low_redshift_all['radius']
swinbank_low_redshift_all_sigma = table_swinbank_low_redshift_all['sig_in']
swinbank_low_redshift_all_sigma_upper_error = table_swinbank_low_redshift_all['sig_in_err']
swinbank_low_redshift_all_sigma_lower_error = table_swinbank_low_redshift_all['sig_in_err']

swinbank_low_redshift_all_dictionary = return_betas(swinbank_low_redshift_all_mass,
                                    swinbank_low_redshift_all_velocity,
                                    swinbank_low_redshift_all_velocity_upper_error,
                                    swinbank_low_redshift_all_velocity_lower_error,
                                    swinbank_low_redshift_all_sigma,
                                    swinbank_low_redshift_all_sigma_upper_error,
                                    swinbank_low_redshift_all_sigma_lower_error,
                                    swinbank_low_redshift_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_low_redshift_all_dictionary['fit_mass']
swinbank_low_redshift_all_log_velocity = swinbank_low_redshift_all_dictionary['log_velocity']
swinbank_low_redshift_all_log_velocity_upper_error = swinbank_low_redshift_all_dictionary['log_velocity_upper_error']
swinbank_low_redshift_all_log_velocity_lower_error = swinbank_low_redshift_all_dictionary['log_velocity_lower_error']
swinbank_low_redshift_all_v_tot = swinbank_low_redshift_all_dictionary['v_tot']
swinbank_low_redshift_all_v_tot_lower_error = swinbank_low_redshift_all_dictionary['v_tot_lower_error']
swinbank_low_redshift_all_v_tot_upper_error = swinbank_low_redshift_all_dictionary['v_tot_upper_error']
swinbank_low_redshift_all_log_v_tot = swinbank_low_redshift_all_dictionary['log_v_tot']
swinbank_low_redshift_all_log_v_tot_lower_error = swinbank_low_redshift_all_dictionary['log_v_tot_lower_error']
swinbank_low_redshift_all_log_v_tot_upper_error = swinbank_low_redshift_all_dictionary['log_v_tot_upper_error']
swinbank_low_redshift_all_log_v_tot_av_error = swinbank_low_redshift_all_dictionary['log_v_tot_av_error']
swinbank_low_redshift_all_ang_mom = swinbank_low_redshift_all_dictionary['ang_mom']
swinbank_low_redshift_all_ang_mom_lower_error = swinbank_low_redshift_all_dictionary['ang_mom_lower_error']
swinbank_low_redshift_all_ang_mom_upper_error = swinbank_low_redshift_all_dictionary['ang_mom_upper_error']
swinbank_low_redshift_all_log_ang_mom = swinbank_low_redshift_all_dictionary['log_ang_mom']
swinbank_low_redshift_all_log_ang_mom_lower_error = swinbank_low_redshift_all_dictionary['log_ang_mom_lower_error']
swinbank_low_redshift_all_log_ang_mom_upper_error = swinbank_low_redshift_all_dictionary['log_ang_mom_upper_error']
swinbank_low_redshift_all_log_ang_mom_av_error = swinbank_low_redshift_all_dictionary['log_ang_mom_av_error']
swinbank_low_redshift_all_ang_mom_tot = swinbank_low_redshift_all_dictionary['ang_mom_tot']
swinbank_low_redshift_all_ang_mom_tot_lower_error = swinbank_low_redshift_all_dictionary['ang_mom_tot_lower_error']
swinbank_low_redshift_all_ang_mom_tot_upper_error = swinbank_low_redshift_all_dictionary['ang_mom_tot_upper_error']
swinbank_low_redshift_all_log_ang_mom_tot = swinbank_low_redshift_all_dictionary['log_ang_mom_tot']
swinbank_low_redshift_all_log_ang_mom_tot_lower_error = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_lower_error']
swinbank_low_redshift_all_log_ang_mom_tot_upper_error = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_upper_error']
swinbank_low_redshift_all_log_ang_mom_tot_av_error = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_av_error']
swinbank_low_redshift_all_log_velocity_mass_fit_params = swinbank_low_redshift_all_dictionary['log_velocity_mass_fit_params']
print 'swinbank_low_redshift_all velocity vs. mass fit params: %s' % swinbank_low_redshift_all_log_velocity_mass_fit_params
swinbank_low_redshift_all_log_velocity_mass_fit_line = swinbank_low_redshift_all_dictionary['log_velocity_mass_fit_line']
swinbank_low_redshift_all_log_velocity_mass_fit_line_lower = swinbank_low_redshift_all_dictionary['log_velocity_mass_fit_line_lower']
swinbank_low_redshift_all_log_velocity_mass_fit_line_upper = swinbank_low_redshift_all_dictionary['log_velocity_mass_fit_line_upper']
swinbank_low_redshift_all_log_v_tot_mass_fit_params = swinbank_low_redshift_all_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_low_redshift_all v_tot vs. mass fit params: %s' % swinbank_low_redshift_all_log_v_tot_mass_fit_params
swinbank_low_redshift_all_log_v_tot_mass_fit_line = swinbank_low_redshift_all_dictionary['log_v_tot_mass_fit_line']
swinbank_low_redshift_all_log_v_tot_mass_fit_line_lower = swinbank_low_redshift_all_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_low_redshift_all_log_v_tot_mass_fit_line_upper = swinbank_low_redshift_all_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_low_redshift_all_log_ang_mom_mass_fit_params = swinbank_low_redshift_all_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_low_redshift_all ang_mom vs. mass fit params: %s' % swinbank_low_redshift_all_log_ang_mom_mass_fit_params
swinbank_low_redshift_all_log_ang_mom_mass_fit_line = swinbank_low_redshift_all_dictionary['log_ang_mom_mass_fit_line']
swinbank_low_redshift_all_log_ang_mom_mass_fit_line_lower = swinbank_low_redshift_all_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_low_redshift_all_log_ang_mom_mass_fit_line_upper = swinbank_low_redshift_all_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_low_redshift_all ang_mom_tot vs. mass fit params: %s' % swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params
swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_line = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_line_lower = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_line_upper = swinbank_low_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_low_redshift_all_median_dynamical_mass_ratio = swinbank_low_redshift_all_dictionary['median_dynamical_mass_ratio']
swinbank_low_redshift_all_median_dynamical_mass_ratio_lower_error = swinbank_low_redshift_all_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_low_redshift_all_median_dynamical_mass_ratio_upper_error = swinbank_low_redshift_all_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio = swinbank_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_low_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SWINBANK LOW REDSHIFT ROT

swinbank_low_redshift_rot_velocity = table_swinbank_low_redshift_rot['V']
swinbank_low_redshift_rot_velocity_upper_error = table_swinbank_low_redshift_rot['V_err']
swinbank_low_redshift_rot_velocity_lower_error = table_swinbank_low_redshift_rot['V_err']
swinbank_low_redshift_rot_mass = table_swinbank_low_redshift_rot['log(M*)']
swinbank_low_redshift_rot_radius = table_swinbank_low_redshift_rot['radius']
swinbank_low_redshift_rot_sigma = table_swinbank_low_redshift_rot['sig_in']
swinbank_low_redshift_rot_sigma_upper_error = table_swinbank_low_redshift_rot['sig_in_err']
swinbank_low_redshift_rot_sigma_lower_error = table_swinbank_low_redshift_rot['sig_in_err']

swinbank_low_redshift_rot_dictionary = return_betas(swinbank_low_redshift_rot_mass,
                                    swinbank_low_redshift_rot_velocity,
                                    swinbank_low_redshift_rot_velocity_upper_error,
                                    swinbank_low_redshift_rot_velocity_lower_error,
                                    swinbank_low_redshift_rot_sigma,
                                    swinbank_low_redshift_rot_sigma_upper_error,
                                    swinbank_low_redshift_rot_sigma_lower_error,
                                    swinbank_low_redshift_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_low_redshift_rot_dictionary['fit_mass']
swinbank_low_redshift_rot_log_velocity = swinbank_low_redshift_rot_dictionary['log_velocity']
swinbank_low_redshift_rot_log_velocity_upper_error = swinbank_low_redshift_rot_dictionary['log_velocity_upper_error']
swinbank_low_redshift_rot_log_velocity_lower_error = swinbank_low_redshift_rot_dictionary['log_velocity_lower_error']
swinbank_low_redshift_rot_v_tot = swinbank_low_redshift_rot_dictionary['v_tot']
swinbank_low_redshift_rot_v_tot_lower_error = swinbank_low_redshift_rot_dictionary['v_tot_lower_error']
swinbank_low_redshift_rot_v_tot_upper_error = swinbank_low_redshift_rot_dictionary['v_tot_upper_error']
swinbank_low_redshift_rot_log_v_tot = swinbank_low_redshift_rot_dictionary['log_v_tot']
swinbank_low_redshift_rot_log_v_tot_lower_error = swinbank_low_redshift_rot_dictionary['log_v_tot_lower_error']
swinbank_low_redshift_rot_log_v_tot_upper_error = swinbank_low_redshift_rot_dictionary['log_v_tot_upper_error']
swinbank_low_redshift_rot_log_v_tot_av_error = swinbank_low_redshift_rot_dictionary['log_v_tot_av_error']
swinbank_low_redshift_rot_ang_mom = swinbank_low_redshift_rot_dictionary['ang_mom']
swinbank_low_redshift_rot_ang_mom_lower_error = swinbank_low_redshift_rot_dictionary['ang_mom_lower_error']
swinbank_low_redshift_rot_ang_mom_upper_error = swinbank_low_redshift_rot_dictionary['ang_mom_upper_error']
swinbank_low_redshift_rot_log_ang_mom = swinbank_low_redshift_rot_dictionary['log_ang_mom']
swinbank_low_redshift_rot_log_ang_mom_lower_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_lower_error']
swinbank_low_redshift_rot_log_ang_mom_upper_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_upper_error']
swinbank_low_redshift_rot_log_ang_mom_av_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_av_error']
swinbank_low_redshift_rot_ang_mom_tot = swinbank_low_redshift_rot_dictionary['ang_mom_tot']
swinbank_low_redshift_rot_ang_mom_tot_lower_error = swinbank_low_redshift_rot_dictionary['ang_mom_tot_lower_error']
swinbank_low_redshift_rot_ang_mom_tot_upper_error = swinbank_low_redshift_rot_dictionary['ang_mom_tot_upper_error']
swinbank_low_redshift_rot_log_ang_mom_tot = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot']
swinbank_low_redshift_rot_log_ang_mom_tot_lower_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_lower_error']
swinbank_low_redshift_rot_log_ang_mom_tot_upper_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_upper_error']
swinbank_low_redshift_rot_log_ang_mom_tot_av_error = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_av_error']
swinbank_low_redshift_rot_log_velocity_mass_fit_params = swinbank_low_redshift_rot_dictionary['log_velocity_mass_fit_params']
print 'swinbank_low_redshift_rot velocity vs. mass fit params: %s' % swinbank_low_redshift_rot_log_velocity_mass_fit_params
swinbank_low_redshift_rot_log_velocity_mass_fit_line = swinbank_low_redshift_rot_dictionary['log_velocity_mass_fit_line']
swinbank_low_redshift_rot_log_velocity_mass_fit_line_lower = swinbank_low_redshift_rot_dictionary['log_velocity_mass_fit_line_lower']
swinbank_low_redshift_rot_log_velocity_mass_fit_line_upper = swinbank_low_redshift_rot_dictionary['log_velocity_mass_fit_line_upper']
swinbank_low_redshift_rot_log_v_tot_mass_fit_params = swinbank_low_redshift_rot_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_low_redshift_rot v_tot vs. mass fit params: %s' % swinbank_low_redshift_rot_log_v_tot_mass_fit_params
swinbank_low_redshift_rot_log_v_tot_mass_fit_line = swinbank_low_redshift_rot_dictionary['log_v_tot_mass_fit_line']
swinbank_low_redshift_rot_log_v_tot_mass_fit_line_lower = swinbank_low_redshift_rot_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_low_redshift_rot_log_v_tot_mass_fit_line_upper = swinbank_low_redshift_rot_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_low_redshift_rot_log_ang_mom_mass_fit_params = swinbank_low_redshift_rot_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_low_redshift_rot ang_mom vs. mass fit params: %s' % swinbank_low_redshift_rot_log_ang_mom_mass_fit_params
swinbank_low_redshift_rot_log_ang_mom_mass_fit_line = swinbank_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line']
swinbank_low_redshift_rot_log_ang_mom_mass_fit_line_lower = swinbank_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_low_redshift_rot_log_ang_mom_mass_fit_line_upper = swinbank_low_redshift_rot_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_low_redshift_rot ang_mom_tot vs. mass fit params: %s' % swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params
swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_line = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_line_lower = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_line_upper = swinbank_low_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_low_redshift_rot_median_dynamical_mass_ratio = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_ratio']
swinbank_low_redshift_rot_median_dynamical_mass_ratio_lower_error = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_low_redshift_rot_median_dynamical_mass_ratio_upper_error = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_low_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SWINBANK LOW REDSHIFT DISP

swinbank_low_redshift_disp_velocity = table_swinbank_low_redshift_disp['V']
swinbank_low_redshift_disp_velocity_upper_error = table_swinbank_low_redshift_disp['V_err']
swinbank_low_redshift_disp_velocity_lower_error = table_swinbank_low_redshift_disp['V_err']
swinbank_low_redshift_disp_mass = table_swinbank_low_redshift_disp['log(M*)']
swinbank_low_redshift_disp_radius = table_swinbank_low_redshift_disp['radius']
swinbank_low_redshift_disp_sigma = table_swinbank_low_redshift_disp['sig_in']
swinbank_low_redshift_disp_sigma_upper_error = table_swinbank_low_redshift_disp['sig_in_err']
swinbank_low_redshift_disp_sigma_lower_error = table_swinbank_low_redshift_disp['sig_in_err']

swinbank_low_redshift_disp_dictionary = return_betas(swinbank_low_redshift_disp_mass,
                                    swinbank_low_redshift_disp_velocity,
                                    swinbank_low_redshift_disp_velocity_upper_error,
                                    swinbank_low_redshift_disp_velocity_lower_error,
                                    swinbank_low_redshift_disp_sigma,
                                    swinbank_low_redshift_disp_sigma_upper_error,
                                    swinbank_low_redshift_disp_sigma_lower_error,
                                    swinbank_low_redshift_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_low_redshift_disp_dictionary['fit_mass']
swinbank_low_redshift_disp_log_velocity = swinbank_low_redshift_disp_dictionary['log_velocity']
swinbank_low_redshift_disp_log_velocity_upper_error = swinbank_low_redshift_disp_dictionary['log_velocity_upper_error']
swinbank_low_redshift_disp_log_velocity_lower_error = swinbank_low_redshift_disp_dictionary['log_velocity_lower_error']
swinbank_low_redshift_disp_v_tot = swinbank_low_redshift_disp_dictionary['v_tot']
swinbank_low_redshift_disp_v_tot_lower_error = swinbank_low_redshift_disp_dictionary['v_tot_lower_error']
swinbank_low_redshift_disp_v_tot_upper_error = swinbank_low_redshift_disp_dictionary['v_tot_upper_error']
swinbank_low_redshift_disp_log_v_tot = swinbank_low_redshift_disp_dictionary['log_v_tot']
swinbank_low_redshift_disp_log_v_tot_lower_error = swinbank_low_redshift_disp_dictionary['log_v_tot_lower_error']
swinbank_low_redshift_disp_log_v_tot_upper_error = swinbank_low_redshift_disp_dictionary['log_v_tot_upper_error']
swinbank_low_redshift_disp_log_v_tot_av_error = swinbank_low_redshift_disp_dictionary['log_v_tot_av_error']
swinbank_low_redshift_disp_ang_mom = swinbank_low_redshift_disp_dictionary['ang_mom']
swinbank_low_redshift_disp_ang_mom_lower_error = swinbank_low_redshift_disp_dictionary['ang_mom_lower_error']
swinbank_low_redshift_disp_ang_mom_upper_error = swinbank_low_redshift_disp_dictionary['ang_mom_upper_error']
swinbank_low_redshift_disp_log_ang_mom = swinbank_low_redshift_disp_dictionary['log_ang_mom']
swinbank_low_redshift_disp_log_ang_mom_lower_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_lower_error']
swinbank_low_redshift_disp_log_ang_mom_upper_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_upper_error']
swinbank_low_redshift_disp_log_ang_mom_av_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_av_error']
swinbank_low_redshift_disp_ang_mom_tot = swinbank_low_redshift_disp_dictionary['ang_mom_tot']
swinbank_low_redshift_disp_ang_mom_tot_lower_error = swinbank_low_redshift_disp_dictionary['ang_mom_tot_lower_error']
swinbank_low_redshift_disp_ang_mom_tot_upper_error = swinbank_low_redshift_disp_dictionary['ang_mom_tot_upper_error']
swinbank_low_redshift_disp_log_ang_mom_tot = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot']
swinbank_low_redshift_disp_log_ang_mom_tot_lower_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_lower_error']
swinbank_low_redshift_disp_log_ang_mom_tot_upper_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_upper_error']
swinbank_low_redshift_disp_log_ang_mom_tot_av_error = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_av_error']
swinbank_low_redshift_disp_log_velocity_mass_fit_params = swinbank_low_redshift_disp_dictionary['log_velocity_mass_fit_params']
print 'swinbank_low_redshift_disp velocity vs. mass fit params: %s' % swinbank_low_redshift_disp_log_velocity_mass_fit_params
swinbank_low_redshift_disp_log_velocity_mass_fit_line = swinbank_low_redshift_disp_dictionary['log_velocity_mass_fit_line']
swinbank_low_redshift_disp_log_velocity_mass_fit_line_lower = swinbank_low_redshift_disp_dictionary['log_velocity_mass_fit_line_lower']
swinbank_low_redshift_disp_log_velocity_mass_fit_line_upper = swinbank_low_redshift_disp_dictionary['log_velocity_mass_fit_line_upper']
swinbank_low_redshift_disp_log_v_tot_mass_fit_params = swinbank_low_redshift_disp_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_low_redshift_disp v_tot vs. mass fit params: %s' % swinbank_low_redshift_disp_log_v_tot_mass_fit_params
swinbank_low_redshift_disp_log_v_tot_mass_fit_line = swinbank_low_redshift_disp_dictionary['log_v_tot_mass_fit_line']
swinbank_low_redshift_disp_log_v_tot_mass_fit_line_lower = swinbank_low_redshift_disp_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_low_redshift_disp_log_v_tot_mass_fit_line_upper = swinbank_low_redshift_disp_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_low_redshift_disp_log_ang_mom_mass_fit_params = swinbank_low_redshift_disp_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_low_redshift_disp ang_mom vs. mass fit params: %s' % swinbank_low_redshift_disp_log_ang_mom_mass_fit_params
swinbank_low_redshift_disp_log_ang_mom_mass_fit_line = swinbank_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line']
swinbank_low_redshift_disp_log_ang_mom_mass_fit_line_lower = swinbank_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_low_redshift_disp_log_ang_mom_mass_fit_line_upper = swinbank_low_redshift_disp_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_low_redshift_disp ang_mom_tot vs. mass fit params: %s' % swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params
swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_line = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_line_lower = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_line_upper = swinbank_low_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_low_redshift_disp_median_dynamical_mass_ratio = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_ratio']
swinbank_low_redshift_disp_median_dynamical_mass_ratio_lower_error = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_low_redshift_disp_median_dynamical_mass_ratio_upper_error = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_low_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['swinbank_low_redshift',
                                0.65,
                                swinbank_low_redshift_all_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_all_log_velocity_mass_fit_params[3],
                                swinbank_low_redshift_all_log_velocity_mass_fit_params[0],
                                swinbank_low_redshift_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_all_log_velocity_mass_fit_params[0] - swinbank_low_redshift_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_all_log_velocity_mass_fit_params[2] - swinbank_low_redshift_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                swinbank_low_redshift_rot_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_rot_log_velocity_mass_fit_params[3],
                                swinbank_low_redshift_rot_log_velocity_mass_fit_params[0],
                                swinbank_low_redshift_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_rot_log_velocity_mass_fit_params[0] - swinbank_low_redshift_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_rot_log_velocity_mass_fit_params[2] - swinbank_low_redshift_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                swinbank_low_redshift_disp_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_disp_log_velocity_mass_fit_params[3],
                                swinbank_low_redshift_disp_log_velocity_mass_fit_params[0],
                                swinbank_low_redshift_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_disp_log_velocity_mass_fit_params[0] - swinbank_low_redshift_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_disp_log_velocity_mass_fit_params[2] - swinbank_low_redshift_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['swinbank_low_redshift',
                                0.65,
                                swinbank_low_redshift_all_log_v_tot_mass_fit_params[4],
                                swinbank_low_redshift_all_log_v_tot_mass_fit_params[3],
                                swinbank_low_redshift_all_log_v_tot_mass_fit_params[0],
                                swinbank_low_redshift_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_all_log_v_tot_mass_fit_params[0] - swinbank_low_redshift_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_all_log_v_tot_mass_fit_params[2] - swinbank_low_redshift_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                swinbank_low_redshift_rot_log_v_tot_mass_fit_params[4],
                                swinbank_low_redshift_rot_log_v_tot_mass_fit_params[3],
                                swinbank_low_redshift_rot_log_v_tot_mass_fit_params[0],
                                swinbank_low_redshift_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_rot_log_v_tot_mass_fit_params[0] - swinbank_low_redshift_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_rot_log_v_tot_mass_fit_params[2] - swinbank_low_redshift_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                swinbank_low_redshift_disp_log_v_tot_mass_fit_params[4],
                                swinbank_low_redshift_disp_log_v_tot_mass_fit_params[3],
                                swinbank_low_redshift_disp_log_v_tot_mass_fit_params[0],
                                swinbank_low_redshift_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_low_redshift_disp_log_v_tot_mass_fit_params[0] - swinbank_low_redshift_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_disp_log_v_tot_mass_fit_params[2] - swinbank_low_redshift_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['swinbank_low_redshift',
                                0.65,
                                swinbank_low_redshift_all_log_ang_mom_mass_fit_params[4],
                                swinbank_low_redshift_all_log_ang_mom_mass_fit_params[3],
                                swinbank_low_redshift_all_log_ang_mom_mass_fit_params[0],
                                swinbank_low_redshift_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_all_log_ang_mom_mass_fit_params[0] - swinbank_low_redshift_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_all_log_ang_mom_mass_fit_params[2] - swinbank_low_redshift_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[4],
                                swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[3],
                                swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[0],
                                swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[0] - swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[2] - swinbank_low_redshift_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[4],
                                swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[3],
                                swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[0],
                                swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[0] - swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[2] - swinbank_low_redshift_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['swinbank_low_redshift',
                                0.65,
                                swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[0] - swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[2] - swinbank_low_redshift_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[2] - swinbank_low_redshift_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[2] - swinbank_low_redshift_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['swinbank_low_redshift',
                                0.65,
                                swinbank_low_redshift_all_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_all_median_dynamical_mass_ratio,
                                swinbank_low_redshift_all_median_dynamical_mass_ratio_lower_error,
                                swinbank_low_redshift_all_median_dynamical_mass_ratio_upper_error,
                                swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio,
                                swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_low_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                swinbank_low_redshift_rot_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_rot_median_dynamical_mass_ratio,
                                swinbank_low_redshift_rot_median_dynamical_mass_ratio_lower_error,
                                swinbank_low_redshift_rot_median_dynamical_mass_ratio_upper_error,
                                swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio,
                                swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_low_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                swinbank_low_redshift_disp_log_velocity_mass_fit_params[4],
                                swinbank_low_redshift_disp_median_dynamical_mass_ratio,
                                swinbank_low_redshift_disp_median_dynamical_mass_ratio_lower_error,
                                swinbank_low_redshift_disp_median_dynamical_mass_ratio_upper_error,
                                swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio,
                                swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_low_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error])

# DYNAMICAL PROPERTIES FOR SWINBANK HIGH REDSHIFT ALL

swinbank_high_redshift_all_velocity = table_swinbank_high_redshift_all['V']
swinbank_high_redshift_all_velocity_upper_error = table_swinbank_high_redshift_all['V_err']
swinbank_high_redshift_all_velocity_lower_error = table_swinbank_high_redshift_all['V_err']
swinbank_high_redshift_all_mass = table_swinbank_high_redshift_all['log(M*)']
swinbank_high_redshift_all_radius = table_swinbank_high_redshift_all['radius']
swinbank_high_redshift_all_sigma = table_swinbank_high_redshift_all['sig_in']
swinbank_high_redshift_all_sigma_upper_error = table_swinbank_high_redshift_all['sig_in_err']
swinbank_high_redshift_all_sigma_lower_error = table_swinbank_high_redshift_all['sig_in_err']

swinbank_high_redshift_all_dictionary = return_betas(swinbank_high_redshift_all_mass,
                                    swinbank_high_redshift_all_velocity,
                                    swinbank_high_redshift_all_velocity_upper_error,
                                    swinbank_high_redshift_all_velocity_lower_error,
                                    swinbank_high_redshift_all_sigma,
                                    swinbank_high_redshift_all_sigma_upper_error,
                                    swinbank_high_redshift_all_sigma_lower_error,
                                    swinbank_high_redshift_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_high_redshift_all_dictionary['fit_mass']
swinbank_high_redshift_all_log_velocity = swinbank_high_redshift_all_dictionary['log_velocity']
swinbank_high_redshift_all_log_velocity_upper_error = swinbank_high_redshift_all_dictionary['log_velocity_upper_error']
swinbank_high_redshift_all_log_velocity_lower_error = swinbank_high_redshift_all_dictionary['log_velocity_lower_error']
swinbank_high_redshift_all_v_tot = swinbank_high_redshift_all_dictionary['v_tot']
swinbank_high_redshift_all_v_tot_lower_error = swinbank_high_redshift_all_dictionary['v_tot_lower_error']
swinbank_high_redshift_all_v_tot_upper_error = swinbank_high_redshift_all_dictionary['v_tot_upper_error']
swinbank_high_redshift_all_log_v_tot = swinbank_high_redshift_all_dictionary['log_v_tot']
swinbank_high_redshift_all_log_v_tot_lower_error = swinbank_high_redshift_all_dictionary['log_v_tot_lower_error']
swinbank_high_redshift_all_log_v_tot_upper_error = swinbank_high_redshift_all_dictionary['log_v_tot_upper_error']
swinbank_high_redshift_all_log_v_tot_av_error = swinbank_high_redshift_all_dictionary['log_v_tot_av_error']
swinbank_high_redshift_all_ang_mom = swinbank_high_redshift_all_dictionary['ang_mom']
swinbank_high_redshift_all_ang_mom_lower_error = swinbank_high_redshift_all_dictionary['ang_mom_lower_error']
swinbank_high_redshift_all_ang_mom_upper_error = swinbank_high_redshift_all_dictionary['ang_mom_upper_error']
swinbank_high_redshift_all_log_ang_mom = swinbank_high_redshift_all_dictionary['log_ang_mom']
swinbank_high_redshift_all_log_ang_mom_lower_error = swinbank_high_redshift_all_dictionary['log_ang_mom_lower_error']
swinbank_high_redshift_all_log_ang_mom_upper_error = swinbank_high_redshift_all_dictionary['log_ang_mom_upper_error']
swinbank_high_redshift_all_log_ang_mom_av_error = swinbank_high_redshift_all_dictionary['log_ang_mom_av_error']
swinbank_high_redshift_all_ang_mom_tot = swinbank_high_redshift_all_dictionary['ang_mom_tot']
swinbank_high_redshift_all_ang_mom_tot_lower_error = swinbank_high_redshift_all_dictionary['ang_mom_tot_lower_error']
swinbank_high_redshift_all_ang_mom_tot_upper_error = swinbank_high_redshift_all_dictionary['ang_mom_tot_upper_error']
swinbank_high_redshift_all_log_ang_mom_tot = swinbank_high_redshift_all_dictionary['log_ang_mom_tot']
swinbank_high_redshift_all_log_ang_mom_tot_lower_error = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_lower_error']
swinbank_high_redshift_all_log_ang_mom_tot_upper_error = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_upper_error']
swinbank_high_redshift_all_log_ang_mom_tot_av_error = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_av_error']
swinbank_high_redshift_all_log_velocity_mass_fit_params = swinbank_high_redshift_all_dictionary['log_velocity_mass_fit_params']
print 'swinbank_high_redshift_all velocity vs. mass fit params: %s' % swinbank_high_redshift_all_log_velocity_mass_fit_params
swinbank_high_redshift_all_log_velocity_mass_fit_line = swinbank_high_redshift_all_dictionary['log_velocity_mass_fit_line']
swinbank_high_redshift_all_log_velocity_mass_fit_line_lower = swinbank_high_redshift_all_dictionary['log_velocity_mass_fit_line_lower']
swinbank_high_redshift_all_log_velocity_mass_fit_line_upper = swinbank_high_redshift_all_dictionary['log_velocity_mass_fit_line_upper']
swinbank_high_redshift_all_log_v_tot_mass_fit_params = swinbank_high_redshift_all_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_high_redshift_all v_tot vs. mass fit params: %s' % swinbank_high_redshift_all_log_v_tot_mass_fit_params
swinbank_high_redshift_all_log_v_tot_mass_fit_line = swinbank_high_redshift_all_dictionary['log_v_tot_mass_fit_line']
swinbank_high_redshift_all_log_v_tot_mass_fit_line_lower = swinbank_high_redshift_all_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_high_redshift_all_log_v_tot_mass_fit_line_upper = swinbank_high_redshift_all_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_high_redshift_all_log_ang_mom_mass_fit_params = swinbank_high_redshift_all_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_high_redshift_all ang_mom vs. mass fit params: %s' % swinbank_high_redshift_all_log_ang_mom_mass_fit_params
swinbank_high_redshift_all_log_ang_mom_mass_fit_line = swinbank_high_redshift_all_dictionary['log_ang_mom_mass_fit_line']
swinbank_high_redshift_all_log_ang_mom_mass_fit_line_lower = swinbank_high_redshift_all_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_high_redshift_all_log_ang_mom_mass_fit_line_upper = swinbank_high_redshift_all_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_high_redshift_all ang_mom_tot vs. mass fit params: %s' % swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params
swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_line = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_line_lower = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_line_upper = swinbank_high_redshift_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_high_redshift_all_median_dynamical_mass_ratio = swinbank_high_redshift_all_dictionary['median_dynamical_mass_ratio']
swinbank_high_redshift_all_median_dynamical_mass_ratio_lower_error = swinbank_high_redshift_all_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_high_redshift_all_median_dynamical_mass_ratio_upper_error = swinbank_high_redshift_all_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio = swinbank_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_high_redshift_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SWINBANK HIGH REDSHIFT ROT

swinbank_high_redshift_rot_velocity = table_swinbank_high_redshift_rot['V']
swinbank_high_redshift_rot_velocity_upper_error = table_swinbank_high_redshift_rot['V_err']
swinbank_high_redshift_rot_velocity_lower_error = table_swinbank_high_redshift_rot['V_err']
swinbank_high_redshift_rot_mass = table_swinbank_high_redshift_rot['log(M*)']
swinbank_high_redshift_rot_radius = table_swinbank_high_redshift_rot['radius']
swinbank_high_redshift_rot_sigma = table_swinbank_high_redshift_rot['sig_in']
swinbank_high_redshift_rot_sigma_upper_error = table_swinbank_high_redshift_rot['sig_in_err']
swinbank_high_redshift_rot_sigma_lower_error = table_swinbank_high_redshift_rot['sig_in_err']

swinbank_high_redshift_rot_dictionary = return_betas(swinbank_high_redshift_rot_mass,
                                    swinbank_high_redshift_rot_velocity,
                                    swinbank_high_redshift_rot_velocity_upper_error,
                                    swinbank_high_redshift_rot_velocity_lower_error,
                                    swinbank_high_redshift_rot_sigma,
                                    swinbank_high_redshift_rot_sigma_upper_error,
                                    swinbank_high_redshift_rot_sigma_lower_error,
                                    swinbank_high_redshift_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_high_redshift_rot_dictionary['fit_mass']
swinbank_high_redshift_rot_log_velocity = swinbank_high_redshift_rot_dictionary['log_velocity']
swinbank_high_redshift_rot_log_velocity_upper_error = swinbank_high_redshift_rot_dictionary['log_velocity_upper_error']
swinbank_high_redshift_rot_log_velocity_lower_error = swinbank_high_redshift_rot_dictionary['log_velocity_lower_error']
swinbank_high_redshift_rot_v_tot = swinbank_high_redshift_rot_dictionary['v_tot']
swinbank_high_redshift_rot_v_tot_lower_error = swinbank_high_redshift_rot_dictionary['v_tot_lower_error']
swinbank_high_redshift_rot_v_tot_upper_error = swinbank_high_redshift_rot_dictionary['v_tot_upper_error']
swinbank_high_redshift_rot_log_v_tot = swinbank_high_redshift_rot_dictionary['log_v_tot']
swinbank_high_redshift_rot_log_v_tot_lower_error = swinbank_high_redshift_rot_dictionary['log_v_tot_lower_error']
swinbank_high_redshift_rot_log_v_tot_upper_error = swinbank_high_redshift_rot_dictionary['log_v_tot_upper_error']
swinbank_high_redshift_rot_log_v_tot_av_error = swinbank_high_redshift_rot_dictionary['log_v_tot_av_error']
swinbank_high_redshift_rot_ang_mom = swinbank_high_redshift_rot_dictionary['ang_mom']
swinbank_high_redshift_rot_ang_mom_lower_error = swinbank_high_redshift_rot_dictionary['ang_mom_lower_error']
swinbank_high_redshift_rot_ang_mom_upper_error = swinbank_high_redshift_rot_dictionary['ang_mom_upper_error']
swinbank_high_redshift_rot_log_ang_mom = swinbank_high_redshift_rot_dictionary['log_ang_mom']
swinbank_high_redshift_rot_log_ang_mom_lower_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_lower_error']
swinbank_high_redshift_rot_log_ang_mom_upper_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_upper_error']
swinbank_high_redshift_rot_log_ang_mom_av_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_av_error']
swinbank_high_redshift_rot_ang_mom_tot = swinbank_high_redshift_rot_dictionary['ang_mom_tot']
swinbank_high_redshift_rot_ang_mom_tot_lower_error = swinbank_high_redshift_rot_dictionary['ang_mom_tot_lower_error']
swinbank_high_redshift_rot_ang_mom_tot_upper_error = swinbank_high_redshift_rot_dictionary['ang_mom_tot_upper_error']
swinbank_high_redshift_rot_log_ang_mom_tot = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot']
swinbank_high_redshift_rot_log_ang_mom_tot_lower_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_lower_error']
swinbank_high_redshift_rot_log_ang_mom_tot_upper_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_upper_error']
swinbank_high_redshift_rot_log_ang_mom_tot_av_error = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_av_error']
swinbank_high_redshift_rot_log_velocity_mass_fit_params = swinbank_high_redshift_rot_dictionary['log_velocity_mass_fit_params']
print 'swinbank_high_redshift_rot velocity vs. mass fit params: %s' % swinbank_high_redshift_rot_log_velocity_mass_fit_params
swinbank_high_redshift_rot_log_velocity_mass_fit_line = swinbank_high_redshift_rot_dictionary['log_velocity_mass_fit_line']
swinbank_high_redshift_rot_log_velocity_mass_fit_line_lower = swinbank_high_redshift_rot_dictionary['log_velocity_mass_fit_line_lower']
swinbank_high_redshift_rot_log_velocity_mass_fit_line_upper = swinbank_high_redshift_rot_dictionary['log_velocity_mass_fit_line_upper']
swinbank_high_redshift_rot_log_v_tot_mass_fit_params = swinbank_high_redshift_rot_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_high_redshift_rot v_tot vs. mass fit params: %s' % swinbank_high_redshift_rot_log_v_tot_mass_fit_params
swinbank_high_redshift_rot_log_v_tot_mass_fit_line = swinbank_high_redshift_rot_dictionary['log_v_tot_mass_fit_line']
swinbank_high_redshift_rot_log_v_tot_mass_fit_line_lower = swinbank_high_redshift_rot_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_high_redshift_rot_log_v_tot_mass_fit_line_upper = swinbank_high_redshift_rot_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_high_redshift_rot_log_ang_mom_mass_fit_params = swinbank_high_redshift_rot_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_high_redshift_rot ang_mom vs. mass fit params: %s' % swinbank_high_redshift_rot_log_ang_mom_mass_fit_params
swinbank_high_redshift_rot_log_ang_mom_mass_fit_line = swinbank_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line']
swinbank_high_redshift_rot_log_ang_mom_mass_fit_line_lower = swinbank_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_high_redshift_rot_log_ang_mom_mass_fit_line_upper = swinbank_high_redshift_rot_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_high_redshift_rot ang_mom_tot vs. mass fit params: %s' % swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params
swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_line = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_line_lower = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_line_upper = swinbank_high_redshift_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_high_redshift_rot_median_dynamical_mass_ratio = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_ratio']
swinbank_high_redshift_rot_median_dynamical_mass_ratio_lower_error = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_high_redshift_rot_median_dynamical_mass_ratio_upper_error = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_high_redshift_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR SWINBANK HIGH REDSHIFT DISP

swinbank_high_redshift_disp_velocity = table_swinbank_high_redshift_disp['V']
swinbank_high_redshift_disp_velocity_upper_error = table_swinbank_high_redshift_disp['V_err']
swinbank_high_redshift_disp_velocity_lower_error = table_swinbank_high_redshift_disp['V_err']
swinbank_high_redshift_disp_mass = table_swinbank_high_redshift_disp['log(M*)']
swinbank_high_redshift_disp_radius = table_swinbank_high_redshift_disp['radius']
swinbank_high_redshift_disp_sigma = table_swinbank_high_redshift_disp['sig_in']
swinbank_high_redshift_disp_sigma_upper_error = table_swinbank_high_redshift_disp['sig_in_err']
swinbank_high_redshift_disp_sigma_lower_error = table_swinbank_high_redshift_disp['sig_in_err']

swinbank_high_redshift_disp_dictionary = return_betas(swinbank_high_redshift_disp_mass,
                                    swinbank_high_redshift_disp_velocity,
                                    swinbank_high_redshift_disp_velocity_upper_error,
                                    swinbank_high_redshift_disp_velocity_lower_error,
                                    swinbank_high_redshift_disp_sigma,
                                    swinbank_high_redshift_disp_sigma_upper_error,
                                    swinbank_high_redshift_disp_sigma_lower_error,
                                    swinbank_high_redshift_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = swinbank_high_redshift_disp_dictionary['fit_mass']
swinbank_high_redshift_disp_log_velocity = swinbank_high_redshift_disp_dictionary['log_velocity']
swinbank_high_redshift_disp_log_velocity_upper_error = swinbank_high_redshift_disp_dictionary['log_velocity_upper_error']
swinbank_high_redshift_disp_log_velocity_lower_error = swinbank_high_redshift_disp_dictionary['log_velocity_lower_error']
swinbank_high_redshift_disp_v_tot = swinbank_high_redshift_disp_dictionary['v_tot']
swinbank_high_redshift_disp_v_tot_lower_error = swinbank_high_redshift_disp_dictionary['v_tot_lower_error']
swinbank_high_redshift_disp_v_tot_upper_error = swinbank_high_redshift_disp_dictionary['v_tot_upper_error']
swinbank_high_redshift_disp_log_v_tot = swinbank_high_redshift_disp_dictionary['log_v_tot']
swinbank_high_redshift_disp_log_v_tot_lower_error = swinbank_high_redshift_disp_dictionary['log_v_tot_lower_error']
swinbank_high_redshift_disp_log_v_tot_upper_error = swinbank_high_redshift_disp_dictionary['log_v_tot_upper_error']
swinbank_high_redshift_disp_log_v_tot_av_error = swinbank_high_redshift_disp_dictionary['log_v_tot_av_error']
swinbank_high_redshift_disp_ang_mom = swinbank_high_redshift_disp_dictionary['ang_mom']
swinbank_high_redshift_disp_ang_mom_lower_error = swinbank_high_redshift_disp_dictionary['ang_mom_lower_error']
swinbank_high_redshift_disp_ang_mom_upper_error = swinbank_high_redshift_disp_dictionary['ang_mom_upper_error']
swinbank_high_redshift_disp_log_ang_mom = swinbank_high_redshift_disp_dictionary['log_ang_mom']
swinbank_high_redshift_disp_log_ang_mom_lower_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_lower_error']
swinbank_high_redshift_disp_log_ang_mom_upper_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_upper_error']
swinbank_high_redshift_disp_log_ang_mom_av_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_av_error']
swinbank_high_redshift_disp_ang_mom_tot = swinbank_high_redshift_disp_dictionary['ang_mom_tot']
swinbank_high_redshift_disp_ang_mom_tot_lower_error = swinbank_high_redshift_disp_dictionary['ang_mom_tot_lower_error']
swinbank_high_redshift_disp_ang_mom_tot_upper_error = swinbank_high_redshift_disp_dictionary['ang_mom_tot_upper_error']
swinbank_high_redshift_disp_log_ang_mom_tot = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot']
swinbank_high_redshift_disp_log_ang_mom_tot_lower_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_lower_error']
swinbank_high_redshift_disp_log_ang_mom_tot_upper_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_upper_error']
swinbank_high_redshift_disp_log_ang_mom_tot_av_error = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_av_error']
swinbank_high_redshift_disp_log_velocity_mass_fit_params = swinbank_high_redshift_disp_dictionary['log_velocity_mass_fit_params']
print 'swinbank_high_redshift_disp velocity vs. mass fit params: %s' % swinbank_high_redshift_disp_log_velocity_mass_fit_params
swinbank_high_redshift_disp_log_velocity_mass_fit_line = swinbank_high_redshift_disp_dictionary['log_velocity_mass_fit_line']
swinbank_high_redshift_disp_log_velocity_mass_fit_line_lower = swinbank_high_redshift_disp_dictionary['log_velocity_mass_fit_line_lower']
swinbank_high_redshift_disp_log_velocity_mass_fit_line_upper = swinbank_high_redshift_disp_dictionary['log_velocity_mass_fit_line_upper']
swinbank_high_redshift_disp_log_v_tot_mass_fit_params = swinbank_high_redshift_disp_dictionary['log_v_tot_mass_fit_params']
print 'swinbank_high_redshift_disp v_tot vs. mass fit params: %s' % swinbank_high_redshift_disp_log_v_tot_mass_fit_params
swinbank_high_redshift_disp_log_v_tot_mass_fit_line = swinbank_high_redshift_disp_dictionary['log_v_tot_mass_fit_line']
swinbank_high_redshift_disp_log_v_tot_mass_fit_line_lower = swinbank_high_redshift_disp_dictionary['log_v_tot_mass_fit_line_lower']
swinbank_high_redshift_disp_log_v_tot_mass_fit_line_upper = swinbank_high_redshift_disp_dictionary['log_v_tot_mass_fit_line_upper']
swinbank_high_redshift_disp_log_ang_mom_mass_fit_params = swinbank_high_redshift_disp_dictionary['log_ang_mom_mass_fit_params']
print 'swinbank_high_redshift_disp ang_mom vs. mass fit params: %s' % swinbank_high_redshift_disp_log_ang_mom_mass_fit_params
swinbank_high_redshift_disp_log_ang_mom_mass_fit_line = swinbank_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line']
swinbank_high_redshift_disp_log_ang_mom_mass_fit_line_lower = swinbank_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line_lower']
swinbank_high_redshift_disp_log_ang_mom_mass_fit_line_upper = swinbank_high_redshift_disp_dictionary['log_ang_mom_mass_fit_line_upper']
swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'swinbank_high_redshift_disp ang_mom_tot vs. mass fit params: %s' % swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params
swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_line = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line']
swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_line_lower = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_line_upper = swinbank_high_redshift_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
swinbank_high_redshift_disp_median_dynamical_mass_ratio = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_ratio']
swinbank_high_redshift_disp_median_dynamical_mass_ratio_lower_error = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_ratio_lower_error']
swinbank_high_redshift_disp_median_dynamical_mass_ratio_upper_error = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_ratio_upper_error']
swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error = swinbank_high_redshift_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['swinbank_high_redshift',
                                1.25,
                                swinbank_high_redshift_all_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_all_log_velocity_mass_fit_params[3],
                                swinbank_high_redshift_all_log_velocity_mass_fit_params[0],
                                swinbank_high_redshift_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_all_log_velocity_mass_fit_params[0] - swinbank_high_redshift_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_all_log_velocity_mass_fit_params[2] - swinbank_high_redshift_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                swinbank_high_redshift_rot_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_rot_log_velocity_mass_fit_params[3],
                                swinbank_high_redshift_rot_log_velocity_mass_fit_params[0],
                                swinbank_high_redshift_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_rot_log_velocity_mass_fit_params[0] - swinbank_high_redshift_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_rot_log_velocity_mass_fit_params[2] - swinbank_high_redshift_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                swinbank_high_redshift_disp_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_disp_log_velocity_mass_fit_params[3],
                                swinbank_high_redshift_disp_log_velocity_mass_fit_params[0],
                                swinbank_high_redshift_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_disp_log_velocity_mass_fit_params[0] - swinbank_high_redshift_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_disp_log_velocity_mass_fit_params[2] - swinbank_high_redshift_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['swinbank_high_redshift',
                                1.25,
                                swinbank_high_redshift_all_log_v_tot_mass_fit_params[4],
                                swinbank_high_redshift_all_log_v_tot_mass_fit_params[3],
                                swinbank_high_redshift_all_log_v_tot_mass_fit_params[0],
                                swinbank_high_redshift_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_all_log_v_tot_mass_fit_params[0] - swinbank_high_redshift_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_all_log_v_tot_mass_fit_params[2] - swinbank_high_redshift_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                swinbank_high_redshift_rot_log_v_tot_mass_fit_params[4],
                                swinbank_high_redshift_rot_log_v_tot_mass_fit_params[3],
                                swinbank_high_redshift_rot_log_v_tot_mass_fit_params[0],
                                swinbank_high_redshift_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_rot_log_v_tot_mass_fit_params[0] - swinbank_high_redshift_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_rot_log_v_tot_mass_fit_params[2] - swinbank_high_redshift_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                swinbank_high_redshift_disp_log_v_tot_mass_fit_params[4],
                                swinbank_high_redshift_disp_log_v_tot_mass_fit_params[3],
                                swinbank_high_redshift_disp_log_v_tot_mass_fit_params[0],
                                swinbank_high_redshift_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((swinbank_high_redshift_disp_log_v_tot_mass_fit_params[0] - swinbank_high_redshift_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_disp_log_v_tot_mass_fit_params[2] - swinbank_high_redshift_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['swinbank_high_redshift',
                                1.25,
                                swinbank_high_redshift_all_log_ang_mom_mass_fit_params[4],
                                swinbank_high_redshift_all_log_ang_mom_mass_fit_params[3],
                                swinbank_high_redshift_all_log_ang_mom_mass_fit_params[0],
                                swinbank_high_redshift_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_all_log_ang_mom_mass_fit_params[0] - swinbank_high_redshift_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_all_log_ang_mom_mass_fit_params[2] - swinbank_high_redshift_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[4],
                                swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[3],
                                swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[0],
                                swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[0] - swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[2] - swinbank_high_redshift_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[4],
                                swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[3],
                                swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[0],
                                swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[0] - swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[2] - swinbank_high_redshift_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['swinbank_high_redshift',
                                1.25,
                                swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[0] - swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[2] - swinbank_high_redshift_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0] - swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[2] - swinbank_high_redshift_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[4],
                                swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[3],
                                swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0],
                                swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0] - swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[2] - swinbank_high_redshift_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['swinbank_high_redshift',
                                1.25,
                                swinbank_high_redshift_all_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_all_median_dynamical_mass_ratio,
                                swinbank_high_redshift_all_median_dynamical_mass_ratio_lower_error,
                                swinbank_high_redshift_all_median_dynamical_mass_ratio_upper_error,
                                swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio,
                                swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_high_redshift_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                swinbank_high_redshift_rot_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_rot_median_dynamical_mass_ratio,
                                swinbank_high_redshift_rot_median_dynamical_mass_ratio_lower_error,
                                swinbank_high_redshift_rot_median_dynamical_mass_ratio_upper_error,
                                swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio,
                                swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_high_redshift_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                swinbank_high_redshift_disp_log_velocity_mass_fit_params[4],
                                swinbank_high_redshift_disp_median_dynamical_mass_ratio,
                                swinbank_high_redshift_disp_median_dynamical_mass_ratio_lower_error,
                                swinbank_high_redshift_disp_median_dynamical_mass_ratio_upper_error,
                                swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio,
                                swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                swinbank_high_redshift_disp_median_dynamical_mass_with_sigma_ratio_upper_error])

# DYNAMICAL PROPERTIES FOR DYNAMO

dynamo_all_velocity = table_dynamo_all['V2.2Rr']
dynamo_all_velocity_upper_error = table_dynamo_all['e_V2.2Rr']
dynamo_all_velocity_lower_error = table_dynamo_all['e_V2.2Rr']
dynamo_all_mass = table_dynamo_all['Mstar']
dynamo_all_radius = table_dynamo_all['radius']
dynamo_all_sigma = table_dynamo_all['sigmac']
dynamo_all_sigma_upper_error = np.repeat(2, len(dynamo_all_sigma))
dynamo_all_sigma_lower_error = np.repeat(2, len(dynamo_all_sigma))

dynamo_all_dictionary = return_betas(dynamo_all_mass,
                                    dynamo_all_velocity,
                                    dynamo_all_velocity_upper_error,
                                    dynamo_all_velocity_lower_error,
                                    dynamo_all_sigma,
                                    dynamo_all_sigma_upper_error,
                                    dynamo_all_sigma_lower_error,
                                    dynamo_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = dynamo_all_dictionary['fit_mass']
dynamo_all_log_velocity = dynamo_all_dictionary['log_velocity']
dynamo_all_log_velocity_upper_error = dynamo_all_dictionary['log_velocity_upper_error']
dynamo_all_log_velocity_lower_error = dynamo_all_dictionary['log_velocity_lower_error']
dynamo_all_v_tot = dynamo_all_dictionary['v_tot']
dynamo_all_v_tot_lower_error = dynamo_all_dictionary['v_tot_lower_error']
dynamo_all_v_tot_upper_error = dynamo_all_dictionary['v_tot_upper_error']
dynamo_all_log_v_tot = dynamo_all_dictionary['log_v_tot']
dynamo_all_log_v_tot_lower_error = dynamo_all_dictionary['log_v_tot_lower_error']
dynamo_all_log_v_tot_upper_error = dynamo_all_dictionary['log_v_tot_upper_error']
dynamo_all_log_v_tot_av_error = dynamo_all_dictionary['log_v_tot_av_error']
dynamo_all_ang_mom = dynamo_all_dictionary['ang_mom']
dynamo_all_ang_mom_lower_error = dynamo_all_dictionary['ang_mom_lower_error']
dynamo_all_ang_mom_upper_error = dynamo_all_dictionary['ang_mom_upper_error']
dynamo_all_log_ang_mom = dynamo_all_dictionary['log_ang_mom']
dynamo_all_log_ang_mom_lower_error = dynamo_all_dictionary['log_ang_mom_lower_error']
dynamo_all_log_ang_mom_upper_error = dynamo_all_dictionary['log_ang_mom_upper_error']
dynamo_all_log_ang_mom_av_error = dynamo_all_dictionary['log_ang_mom_av_error']
dynamo_all_ang_mom_tot = dynamo_all_dictionary['ang_mom_tot']
dynamo_all_ang_mom_tot_lower_error = dynamo_all_dictionary['ang_mom_tot_lower_error']
dynamo_all_ang_mom_tot_upper_error = dynamo_all_dictionary['ang_mom_tot_upper_error']
dynamo_all_log_ang_mom_tot = dynamo_all_dictionary['log_ang_mom_tot']
dynamo_all_log_ang_mom_tot_lower_error = dynamo_all_dictionary['log_ang_mom_tot_lower_error']
dynamo_all_log_ang_mom_tot_upper_error = dynamo_all_dictionary['log_ang_mom_tot_upper_error']
dynamo_all_log_ang_mom_tot_av_error = dynamo_all_dictionary['log_ang_mom_tot_av_error']
dynamo_all_log_velocity_mass_fit_params = dynamo_all_dictionary['log_velocity_mass_fit_params']
print 'dynamo_all velocity vs. mass fit params: %s' % dynamo_all_log_velocity_mass_fit_params
dynamo_all_log_velocity_mass_fit_line = dynamo_all_dictionary['log_velocity_mass_fit_line']
dynamo_all_log_velocity_mass_fit_line_lower = dynamo_all_dictionary['log_velocity_mass_fit_line_lower']
dynamo_all_log_velocity_mass_fit_line_upper = dynamo_all_dictionary['log_velocity_mass_fit_line_upper']
dynamo_all_log_v_tot_mass_fit_params = dynamo_all_dictionary['log_v_tot_mass_fit_params']
print 'dynamo_all v_tot vs. mass fit params: %s' % dynamo_all_log_v_tot_mass_fit_params
dynamo_all_log_v_tot_mass_fit_line = dynamo_all_dictionary['log_v_tot_mass_fit_line']
dynamo_all_log_v_tot_mass_fit_line_lower = dynamo_all_dictionary['log_v_tot_mass_fit_line_lower']
dynamo_all_log_v_tot_mass_fit_line_upper = dynamo_all_dictionary['log_v_tot_mass_fit_line_upper']
dynamo_all_log_ang_mom_mass_fit_params = dynamo_all_dictionary['log_ang_mom_mass_fit_params']
print 'dynamo_all ang_mom vs. mass fit params: %s' % dynamo_all_log_ang_mom_mass_fit_params
dynamo_all_log_ang_mom_mass_fit_line = dynamo_all_dictionary['log_ang_mom_mass_fit_line']
dynamo_all_log_ang_mom_mass_fit_line_lower = dynamo_all_dictionary['log_ang_mom_mass_fit_line_lower']
dynamo_all_log_ang_mom_mass_fit_line_upper = dynamo_all_dictionary['log_ang_mom_mass_fit_line_upper']
dynamo_all_log_ang_mom_tot_mass_fit_params = dynamo_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'dynamo_all ang_mom_tot vs. mass fit params: %s' % dynamo_all_log_ang_mom_tot_mass_fit_params
dynamo_all_log_ang_mom_tot_mass_fit_line = dynamo_all_dictionary['log_ang_mom_tot_mass_fit_line']
dynamo_all_log_ang_mom_tot_mass_fit_line_lower = dynamo_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
dynamo_all_log_ang_mom_tot_mass_fit_line_upper = dynamo_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
dynamo_all_median_dynamical_mass_ratio = dynamo_all_dictionary['median_dynamical_mass_ratio']
dynamo_all_median_dynamical_mass_ratio_lower_error = dynamo_all_dictionary['median_dynamical_mass_ratio_lower_error']
dynamo_all_median_dynamical_mass_ratio_upper_error = dynamo_all_dictionary['median_dynamical_mass_ratio_upper_error']
dynamo_all_median_dynamical_mass_with_sigma_ratio = dynamo_all_dictionary['median_dynamical_mass_with_sigma_ratio']
dynamo_all_median_dynamical_mass_with_sigma_ratio_lower_error = dynamo_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
dynamo_all_median_dynamical_mass_with_sigma_ratio_upper_error = dynamo_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['dynamo',
                                0.1,
                                dynamo_all_log_velocity_mass_fit_params[4],
                                dynamo_all_log_velocity_mass_fit_params[3],
                                dynamo_all_log_velocity_mass_fit_params[0],
                                dynamo_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((dynamo_all_log_velocity_mass_fit_params[0] - dynamo_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((dynamo_all_log_velocity_mass_fit_params[2] - dynamo_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
v_tot_evolution_list.append(['dynamo',
                                0.1,
                                dynamo_all_log_v_tot_mass_fit_params[4],
                                dynamo_all_log_v_tot_mass_fit_params[3],
                                dynamo_all_log_v_tot_mass_fit_params[0],
                                dynamo_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((dynamo_all_log_v_tot_mass_fit_params[0] - dynamo_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((dynamo_all_log_v_tot_mass_fit_params[2] - dynamo_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_evolution_list.append(['dynamo',
                                0.1,
                                dynamo_all_log_ang_mom_mass_fit_params[4],
                                dynamo_all_log_ang_mom_mass_fit_params[3],
                                dynamo_all_log_ang_mom_mass_fit_params[0],
                                dynamo_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((dynamo_all_log_ang_mom_mass_fit_params[0] - dynamo_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((dynamo_all_log_ang_mom_mass_fit_params[2] - dynamo_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_tot_evolution_list.append(['dynamo',
                                0.1,
                                dynamo_all_log_ang_mom_tot_mass_fit_params[4],
                                dynamo_all_log_ang_mom_tot_mass_fit_params[3],
                                dynamo_all_log_ang_mom_tot_mass_fit_params[0],
                                dynamo_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((dynamo_all_log_ang_mom_tot_mass_fit_params[0] - dynamo_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((dynamo_all_log_ang_mom_tot_mass_fit_params[2] - dynamo_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
dyn_mass_evolution_list.append(['dynamo',
                                0.1,
                                dynamo_all_log_velocity_mass_fit_params[4],
                                dynamo_all_median_dynamical_mass_ratio,
                                dynamo_all_median_dynamical_mass_ratio_lower_error,
                                dynamo_all_median_dynamical_mass_ratio_upper_error,
                                dynamo_all_median_dynamical_mass_with_sigma_ratio,
                                dynamo_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                dynamo_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])

# DYNAMICAL PROPERTIES FOR KROSS ALL

kross_all_velocity = table_kross_all['VC']
kross_all_velocity_upper_error = table_kross_all['VC_ERR_H']
kross_all_velocity_lower_error = table_kross_all['VC_ERR_L']
kross_all_mass = np.log10(table_kross_all['MASS'])
kross_all_radius = table_kross_all['R_IM']
kross_all_sigma = table_kross_all['SIGMA0']
kross_all_sigma_upper_error = table_kross_all['SIGMA0_ERR']
kross_all_sigma_lower_error = table_kross_all['SIGMA0_ERR']

kross_all_dictionary = return_betas(kross_all_mass,
                                    kross_all_velocity,
                                    kross_all_velocity_upper_error,
                                    kross_all_velocity_lower_error,
                                    kross_all_sigma,
                                    kross_all_sigma_upper_error,
                                    kross_all_sigma_lower_error,
                                    kross_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kross_all_dictionary['fit_mass']
kross_all_log_velocity = kross_all_dictionary['log_velocity']
kross_all_log_velocity_upper_error = kross_all_dictionary['log_velocity_upper_error']
kross_all_log_velocity_lower_error = kross_all_dictionary['log_velocity_lower_error']
kross_all_v_tot = kross_all_dictionary['v_tot']
kross_all_v_tot_lower_error = kross_all_dictionary['v_tot_lower_error']
kross_all_v_tot_upper_error = kross_all_dictionary['v_tot_upper_error']
kross_all_log_v_tot = kross_all_dictionary['log_v_tot']
kross_all_log_v_tot_lower_error = kross_all_dictionary['log_v_tot_lower_error']
kross_all_log_v_tot_upper_error = kross_all_dictionary['log_v_tot_upper_error']
kross_all_log_v_tot_av_error = kross_all_dictionary['log_v_tot_av_error']
kross_all_ang_mom = kross_all_dictionary['ang_mom']
kross_all_ang_mom_lower_error = kross_all_dictionary['ang_mom_lower_error']
kross_all_ang_mom_upper_error = kross_all_dictionary['ang_mom_upper_error']
kross_all_log_ang_mom = kross_all_dictionary['log_ang_mom']
kross_all_log_ang_mom_lower_error = kross_all_dictionary['log_ang_mom_lower_error']
kross_all_log_ang_mom_upper_error = kross_all_dictionary['log_ang_mom_upper_error']
kross_all_log_ang_mom_av_error = kross_all_dictionary['log_ang_mom_av_error']
kross_all_ang_mom_tot = kross_all_dictionary['ang_mom_tot']
kross_all_ang_mom_tot_lower_error = kross_all_dictionary['ang_mom_tot_lower_error']
kross_all_ang_mom_tot_upper_error = kross_all_dictionary['ang_mom_tot_upper_error']
kross_all_log_ang_mom_tot = kross_all_dictionary['log_ang_mom_tot']
kross_all_log_ang_mom_tot_lower_error = kross_all_dictionary['log_ang_mom_tot_lower_error']
kross_all_log_ang_mom_tot_upper_error = kross_all_dictionary['log_ang_mom_tot_upper_error']
kross_all_log_ang_mom_tot_av_error = kross_all_dictionary['log_ang_mom_tot_av_error']
kross_all_log_velocity_mass_fit_params = kross_all_dictionary['log_velocity_mass_fit_params']
print 'kross_all velocity vs. mass fit params: %s' % kross_all_log_velocity_mass_fit_params
kross_all_log_velocity_mass_fit_line = kross_all_dictionary['log_velocity_mass_fit_line']
kross_all_log_velocity_mass_fit_line_lower = kross_all_dictionary['log_velocity_mass_fit_line_lower']
kross_all_log_velocity_mass_fit_line_upper = kross_all_dictionary['log_velocity_mass_fit_line_upper']
kross_all_log_v_tot_mass_fit_params = kross_all_dictionary['log_v_tot_mass_fit_params']
print 'kross_all v_tot vs. mass fit params: %s' % kross_all_log_v_tot_mass_fit_params
kross_all_log_v_tot_mass_fit_line = kross_all_dictionary['log_v_tot_mass_fit_line']
kross_all_log_v_tot_mass_fit_line_lower = kross_all_dictionary['log_v_tot_mass_fit_line_lower']
kross_all_log_v_tot_mass_fit_line_upper = kross_all_dictionary['log_v_tot_mass_fit_line_upper']
kross_all_log_ang_mom_mass_fit_params = kross_all_dictionary['log_ang_mom_mass_fit_params']
print 'kross_all ang_mom vs. mass fit params: %s' % kross_all_log_ang_mom_mass_fit_params
kross_all_log_ang_mom_mass_fit_line = kross_all_dictionary['log_ang_mom_mass_fit_line']
kross_all_log_ang_mom_mass_fit_line_lower = kross_all_dictionary['log_ang_mom_mass_fit_line_lower']
kross_all_log_ang_mom_mass_fit_line_upper = kross_all_dictionary['log_ang_mom_mass_fit_line_upper']
kross_all_log_ang_mom_tot_mass_fit_params = kross_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kross_all ang_mom_tot vs. mass fit params: %s' % kross_all_log_ang_mom_tot_mass_fit_params
kross_all_log_ang_mom_tot_mass_fit_line = kross_all_dictionary['log_ang_mom_tot_mass_fit_line']
kross_all_log_ang_mom_tot_mass_fit_line_lower = kross_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kross_all_log_ang_mom_tot_mass_fit_line_upper = kross_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kross_all_median_dynamical_mass_ratio = kross_all_dictionary['median_dynamical_mass_ratio']
kross_all_median_dynamical_mass_ratio_lower_error = kross_all_dictionary['median_dynamical_mass_ratio_lower_error']
kross_all_median_dynamical_mass_ratio_upper_error = kross_all_dictionary['median_dynamical_mass_ratio_upper_error']
kross_all_median_dynamical_mass_with_sigma_ratio = kross_all_dictionary['median_dynamical_mass_with_sigma_ratio']
kross_all_median_dynamical_mass_with_sigma_ratio_lower_error = kross_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kross_all_median_dynamical_mass_with_sigma_ratio_upper_error = kross_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR KROSS ROT

kross_rot_velocity = table_kross_rot['VC']
kross_rot_velocity_upper_error = table_kross_rot['VC_ERR_H']
kross_rot_velocity_lower_error = table_kross_rot['VC_ERR_L']
kross_rot_mass = np.log10(table_kross_rot['MASS'])
kross_rot_radius = table_kross_rot['R_IM']
kross_rot_sigma = table_kross_rot['SIGMA0']
kross_rot_sigma_upper_error = table_kross_rot['SIGMA0_ERR']
kross_rot_sigma_lower_error = table_kross_rot['SIGMA0_ERR']

kross_rot_dictionary = return_betas(kross_rot_mass,
                                    kross_rot_velocity,
                                    kross_rot_velocity_upper_error,
                                    kross_rot_velocity_lower_error,
                                    kross_rot_sigma,
                                    kross_rot_sigma_upper_error,
                                    kross_rot_sigma_lower_error,
                                    kross_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kross_rot_dictionary['fit_mass']
kross_rot_log_velocity = kross_rot_dictionary['log_velocity']
kross_rot_log_velocity_upper_error = kross_rot_dictionary['log_velocity_upper_error']
kross_rot_log_velocity_lower_error = kross_rot_dictionary['log_velocity_lower_error']
kross_rot_v_tot = kross_rot_dictionary['v_tot']
kross_rot_v_tot_lower_error = kross_rot_dictionary['v_tot_lower_error']
kross_rot_v_tot_upper_error = kross_rot_dictionary['v_tot_upper_error']
kross_rot_log_v_tot = kross_rot_dictionary['log_v_tot']
kross_rot_log_v_tot_lower_error = kross_rot_dictionary['log_v_tot_lower_error']
kross_rot_log_v_tot_upper_error = kross_rot_dictionary['log_v_tot_upper_error']
kross_rot_log_v_tot_av_error = kross_rot_dictionary['log_v_tot_av_error']
kross_rot_ang_mom = kross_rot_dictionary['ang_mom']
kross_rot_ang_mom_lower_error = kross_rot_dictionary['ang_mom_lower_error']
kross_rot_ang_mom_upper_error = kross_rot_dictionary['ang_mom_upper_error']
kross_rot_log_ang_mom = kross_rot_dictionary['log_ang_mom']
kross_rot_log_ang_mom_lower_error = kross_rot_dictionary['log_ang_mom_lower_error']
kross_rot_log_ang_mom_upper_error = kross_rot_dictionary['log_ang_mom_upper_error']
kross_rot_log_ang_mom_av_error = kross_rot_dictionary['log_ang_mom_av_error']
kross_rot_ang_mom_tot = kross_rot_dictionary['ang_mom_tot']
kross_rot_ang_mom_tot_lower_error = kross_rot_dictionary['ang_mom_tot_lower_error']
kross_rot_ang_mom_tot_upper_error = kross_rot_dictionary['ang_mom_tot_upper_error']
kross_rot_log_ang_mom_tot = kross_rot_dictionary['log_ang_mom_tot']
kross_rot_log_ang_mom_tot_lower_error = kross_rot_dictionary['log_ang_mom_tot_lower_error']
kross_rot_log_ang_mom_tot_upper_error = kross_rot_dictionary['log_ang_mom_tot_upper_error']
kross_rot_log_ang_mom_tot_av_error = kross_rot_dictionary['log_ang_mom_tot_av_error']
kross_rot_log_velocity_mass_fit_params = kross_rot_dictionary['log_velocity_mass_fit_params']
print 'kross_rot velocity vs. mass fit params: %s' % kross_rot_log_velocity_mass_fit_params
kross_rot_log_velocity_mass_fit_line = kross_rot_dictionary['log_velocity_mass_fit_line']
kross_rot_log_velocity_mass_fit_line_lower = kross_rot_dictionary['log_velocity_mass_fit_line_lower']
kross_rot_log_velocity_mass_fit_line_upper = kross_rot_dictionary['log_velocity_mass_fit_line_upper']
kross_rot_log_v_tot_mass_fit_params = kross_rot_dictionary['log_v_tot_mass_fit_params']
print 'kross_rot v_tot vs. mass fit params: %s' % kross_rot_log_v_tot_mass_fit_params
kross_rot_log_v_tot_mass_fit_line = kross_rot_dictionary['log_v_tot_mass_fit_line']
kross_rot_log_v_tot_mass_fit_line_lower = kross_rot_dictionary['log_v_tot_mass_fit_line_lower']
kross_rot_log_v_tot_mass_fit_line_upper = kross_rot_dictionary['log_v_tot_mass_fit_line_upper']
kross_rot_log_ang_mom_mass_fit_params = kross_rot_dictionary['log_ang_mom_mass_fit_params']
print 'kross_rot ang_mom vs. mass fit params: %s' % kross_rot_log_ang_mom_mass_fit_params
kross_rot_log_ang_mom_mass_fit_line = kross_rot_dictionary['log_ang_mom_mass_fit_line']
kross_rot_log_ang_mom_mass_fit_line_lower = kross_rot_dictionary['log_ang_mom_mass_fit_line_lower']
kross_rot_log_ang_mom_mass_fit_line_upper = kross_rot_dictionary['log_ang_mom_mass_fit_line_upper']
kross_rot_log_ang_mom_tot_mass_fit_params = kross_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kross_rot ang_mom_tot vs. mass fit params: %s' % kross_rot_log_ang_mom_tot_mass_fit_params
kross_rot_log_ang_mom_tot_mass_fit_line = kross_rot_dictionary['log_ang_mom_tot_mass_fit_line']
kross_rot_log_ang_mom_tot_mass_fit_line_lower = kross_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kross_rot_log_ang_mom_tot_mass_fit_line_upper = kross_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kross_rot_median_dynamical_mass_ratio = kross_rot_dictionary['median_dynamical_mass_ratio']
kross_rot_median_dynamical_mass_ratio_lower_error = kross_rot_dictionary['median_dynamical_mass_ratio_lower_error']
kross_rot_median_dynamical_mass_ratio_upper_error = kross_rot_dictionary['median_dynamical_mass_ratio_upper_error']
kross_rot_median_dynamical_mass_with_sigma_ratio = kross_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
kross_rot_median_dynamical_mass_with_sigma_ratio_lower_error = kross_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kross_rot_median_dynamical_mass_with_sigma_ratio_upper_error = kross_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR KROSS DISP

kross_disp_velocity = table_kross_disp['VC']
kross_disp_velocity_upper_error = table_kross_disp['VC_ERR_H']
kross_disp_velocity_lower_error = table_kross_disp['VC_ERR_L']
kross_disp_mass = np.log10(table_kross_disp['MASS'])
kross_disp_radius = table_kross_disp['R_IM']
kross_disp_sigma = table_kross_disp['SIGMA0']
kross_disp_sigma_upper_error = table_kross_disp['SIGMA0_ERR']
kross_disp_sigma_lower_error = table_kross_disp['SIGMA0_ERR']

kross_disp_dictionary = return_betas(kross_disp_mass,
                                    kross_disp_velocity,
                                    kross_disp_velocity_upper_error,
                                    kross_disp_velocity_lower_error,
                                    kross_disp_sigma,
                                    kross_disp_sigma_upper_error,
                                    kross_disp_sigma_lower_error,
                                    kross_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kross_disp_dictionary['fit_mass']
kross_disp_log_velocity = kross_disp_dictionary['log_velocity']
kross_disp_log_velocity_upper_error = kross_disp_dictionary['log_velocity_upper_error']
kross_disp_log_velocity_lower_error = kross_disp_dictionary['log_velocity_lower_error']
kross_disp_v_tot = kross_disp_dictionary['v_tot']
kross_disp_v_tot_lower_error = kross_disp_dictionary['v_tot_lower_error']
kross_disp_v_tot_upper_error = kross_disp_dictionary['v_tot_upper_error']
kross_disp_log_v_tot = kross_disp_dictionary['log_v_tot']
kross_disp_log_v_tot_lower_error = kross_disp_dictionary['log_v_tot_lower_error']
kross_disp_log_v_tot_upper_error = kross_disp_dictionary['log_v_tot_upper_error']
kross_disp_log_v_tot_av_error = kross_disp_dictionary['log_v_tot_av_error']
kross_disp_ang_mom = kross_disp_dictionary['ang_mom']
kross_disp_ang_mom_lower_error = kross_disp_dictionary['ang_mom_lower_error']
kross_disp_ang_mom_upper_error = kross_disp_dictionary['ang_mom_upper_error']
kross_disp_log_ang_mom = kross_disp_dictionary['log_ang_mom']
kross_disp_log_ang_mom_lower_error = kross_disp_dictionary['log_ang_mom_lower_error']
kross_disp_log_ang_mom_upper_error = kross_disp_dictionary['log_ang_mom_upper_error']
kross_disp_log_ang_mom_av_error = kross_disp_dictionary['log_ang_mom_av_error']
kross_disp_ang_mom_tot = kross_disp_dictionary['ang_mom_tot']
kross_disp_ang_mom_tot_lower_error = kross_disp_dictionary['ang_mom_tot_lower_error']
kross_disp_ang_mom_tot_upper_error = kross_disp_dictionary['ang_mom_tot_upper_error']
kross_disp_log_ang_mom_tot = kross_disp_dictionary['log_ang_mom_tot']
kross_disp_log_ang_mom_tot_lower_error = kross_disp_dictionary['log_ang_mom_tot_lower_error']
kross_disp_log_ang_mom_tot_upper_error = kross_disp_dictionary['log_ang_mom_tot_upper_error']
kross_disp_log_ang_mom_tot_av_error = kross_disp_dictionary['log_ang_mom_tot_av_error']
kross_disp_log_velocity_mass_fit_params = kross_disp_dictionary['log_velocity_mass_fit_params']
print 'kross_disp velocity vs. mass fit params: %s' % kross_disp_log_velocity_mass_fit_params
kross_disp_log_velocity_mass_fit_line = kross_disp_dictionary['log_velocity_mass_fit_line']
kross_disp_log_velocity_mass_fit_line_lower = kross_disp_dictionary['log_velocity_mass_fit_line_lower']
kross_disp_log_velocity_mass_fit_line_upper = kross_disp_dictionary['log_velocity_mass_fit_line_upper']
kross_disp_log_v_tot_mass_fit_params = kross_disp_dictionary['log_v_tot_mass_fit_params']
print 'kross_disp v_tot vs. mass fit params: %s' % kross_disp_log_v_tot_mass_fit_params
kross_disp_log_v_tot_mass_fit_line = kross_disp_dictionary['log_v_tot_mass_fit_line']
kross_disp_log_v_tot_mass_fit_line_lower = kross_disp_dictionary['log_v_tot_mass_fit_line_lower']
kross_disp_log_v_tot_mass_fit_line_upper = kross_disp_dictionary['log_v_tot_mass_fit_line_upper']
kross_disp_log_ang_mom_mass_fit_params = kross_disp_dictionary['log_ang_mom_mass_fit_params']
print 'kross_disp ang_mom vs. mass fit params: %s' % kross_disp_log_ang_mom_mass_fit_params
kross_disp_log_ang_mom_mass_fit_line = kross_disp_dictionary['log_ang_mom_mass_fit_line']
kross_disp_log_ang_mom_mass_fit_line_lower = kross_disp_dictionary['log_ang_mom_mass_fit_line_lower']
kross_disp_log_ang_mom_mass_fit_line_upper = kross_disp_dictionary['log_ang_mom_mass_fit_line_upper']
kross_disp_log_ang_mom_tot_mass_fit_params = kross_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kross_disp ang_mom_tot vs. mass fit params: %s' % kross_disp_log_ang_mom_tot_mass_fit_params
kross_disp_log_ang_mom_tot_mass_fit_line = kross_disp_dictionary['log_ang_mom_tot_mass_fit_line']
kross_disp_log_ang_mom_tot_mass_fit_line_lower = kross_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kross_disp_log_ang_mom_tot_mass_fit_line_upper = kross_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kross_disp_median_dynamical_mass_ratio = kross_disp_dictionary['median_dynamical_mass_ratio']
kross_disp_median_dynamical_mass_ratio_lower_error = kross_disp_dictionary['median_dynamical_mass_ratio_lower_error']
kross_disp_median_dynamical_mass_ratio_upper_error = kross_disp_dictionary['median_dynamical_mass_ratio_upper_error']
kross_disp_median_dynamical_mass_with_sigma_ratio = kross_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
kross_disp_median_dynamical_mass_with_sigma_ratio_lower_error = kross_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kross_disp_median_dynamical_mass_with_sigma_ratio_upper_error = kross_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['kross',
                                0.9,
                                kross_all_log_velocity_mass_fit_params[4],
                                kross_all_log_velocity_mass_fit_params[3],
                                kross_all_log_velocity_mass_fit_params[0],
                                kross_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_all_log_velocity_mass_fit_params[0] - kross_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kross_all_log_velocity_mass_fit_params[2] - kross_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                kross_rot_log_velocity_mass_fit_params[4],
                                kross_rot_log_velocity_mass_fit_params[3],
                                kross_rot_log_velocity_mass_fit_params[0],
                                kross_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_rot_log_velocity_mass_fit_params[0] - kross_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kross_rot_log_velocity_mass_fit_params[2] - kross_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                kross_disp_log_velocity_mass_fit_params[4],
                                kross_disp_log_velocity_mass_fit_params[3],
                                kross_disp_log_velocity_mass_fit_params[0],
                                kross_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_disp_log_velocity_mass_fit_params[0] - kross_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kross_disp_log_velocity_mass_fit_params[2] - kross_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['kross',
                                0.9,
                                kross_all_log_v_tot_mass_fit_params[4],
                                kross_all_log_v_tot_mass_fit_params[3],
                                kross_all_log_v_tot_mass_fit_params[0],
                                kross_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_all_log_v_tot_mass_fit_params[0] - kross_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kross_all_log_v_tot_mass_fit_params[2] - kross_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                kross_rot_log_v_tot_mass_fit_params[4],
                                kross_rot_log_v_tot_mass_fit_params[3],
                                kross_rot_log_v_tot_mass_fit_params[0],
                                kross_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_rot_log_v_tot_mass_fit_params[0] - kross_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kross_rot_log_v_tot_mass_fit_params[2] - kross_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                kross_disp_log_v_tot_mass_fit_params[4],
                                kross_disp_log_v_tot_mass_fit_params[3],
                                kross_disp_log_v_tot_mass_fit_params[0],
                                kross_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kross_disp_log_v_tot_mass_fit_params[0] - kross_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kross_disp_log_v_tot_mass_fit_params[2] - kross_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['kross',
                                0.9,
                                kross_all_log_ang_mom_mass_fit_params[4],
                                kross_all_log_ang_mom_mass_fit_params[3],
                                kross_all_log_ang_mom_mass_fit_params[0],
                                kross_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_all_log_ang_mom_mass_fit_params[0] - kross_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kross_all_log_ang_mom_mass_fit_params[2] - kross_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                kross_rot_log_ang_mom_mass_fit_params[4],
                                kross_rot_log_ang_mom_mass_fit_params[3],
                                kross_rot_log_ang_mom_mass_fit_params[0],
                                kross_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_rot_log_ang_mom_mass_fit_params[0] - kross_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kross_rot_log_ang_mom_mass_fit_params[2] - kross_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                kross_disp_log_ang_mom_mass_fit_params[4],
                                kross_disp_log_ang_mom_mass_fit_params[3],
                                kross_disp_log_ang_mom_mass_fit_params[0],
                                kross_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_disp_log_ang_mom_mass_fit_params[0] - kross_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kross_disp_log_ang_mom_mass_fit_params[2] - kross_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['kross',
                                0.9,
                                kross_all_log_ang_mom_tot_mass_fit_params[4],
                                kross_all_log_ang_mom_tot_mass_fit_params[3],
                                kross_all_log_ang_mom_tot_mass_fit_params[0],
                                kross_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_all_log_ang_mom_tot_mass_fit_params[0] - kross_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kross_all_log_ang_mom_tot_mass_fit_params[2] - kross_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                kross_rot_log_ang_mom_tot_mass_fit_params[4],
                                kross_rot_log_ang_mom_tot_mass_fit_params[3],
                                kross_rot_log_ang_mom_tot_mass_fit_params[0],
                                kross_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_rot_log_ang_mom_tot_mass_fit_params[0] - kross_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kross_rot_log_ang_mom_tot_mass_fit_params[2] - kross_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                kross_disp_log_ang_mom_tot_mass_fit_params[4],
                                kross_disp_log_ang_mom_tot_mass_fit_params[3],
                                kross_disp_log_ang_mom_tot_mass_fit_params[0],
                                kross_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kross_disp_log_ang_mom_tot_mass_fit_params[0] - kross_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kross_disp_log_ang_mom_tot_mass_fit_params[2] - kross_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['kross',
                                0.9,
                                kross_all_log_velocity_mass_fit_params[4],
                                kross_all_median_dynamical_mass_ratio,
                                kross_all_median_dynamical_mass_ratio_lower_error,
                                kross_all_median_dynamical_mass_ratio_upper_error,
                                kross_all_median_dynamical_mass_with_sigma_ratio,
                                kross_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kross_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                kross_rot_log_velocity_mass_fit_params[4],
                                kross_rot_median_dynamical_mass_ratio,
                                kross_rot_median_dynamical_mass_ratio_lower_error,
                                kross_rot_median_dynamical_mass_ratio_upper_error,
                                kross_rot_median_dynamical_mass_with_sigma_ratio,
                                kross_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kross_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                kross_disp_log_velocity_mass_fit_params[4],
                                kross_disp_median_dynamical_mass_ratio,
                                kross_disp_median_dynamical_mass_ratio_lower_error,
                                kross_disp_median_dynamical_mass_ratio_upper_error,
                                kross_disp_median_dynamical_mass_with_sigma_ratio,
                                kross_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kross_disp_median_dynamical_mass_with_sigma_ratio_upper_error])

# DYNAMICAL PROPERTIES FOR MASSIV ALL

massiv_all_velocity = table_massiv_all['Vmax']
massiv_all_velocity_upper_error = table_massiv_all['vmax_error']
massiv_all_velocity_lower_error = table_massiv_all['vmax_error']
massiv_all_mass = table_massiv_all['M_star']
massiv_all_radius = table_massiv_all['Re']
massiv_all_sigma = table_massiv_all['sigma']
massiv_all_sigma_upper_error = table_massiv_all['sigma_error']
massiv_all_sigma_lower_error = table_massiv_all['sigma_error']

massiv_all_dictionary = return_betas(massiv_all_mass,
                                    massiv_all_velocity,
                                    massiv_all_velocity_upper_error,
                                    massiv_all_velocity_lower_error,
                                    massiv_all_sigma,
                                    massiv_all_sigma_upper_error,
                                    massiv_all_sigma_lower_error,
                                    massiv_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = massiv_all_dictionary['fit_mass']
massiv_all_log_velocity = massiv_all_dictionary['log_velocity']
massiv_all_log_velocity_upper_error = massiv_all_dictionary['log_velocity_upper_error']
massiv_all_log_velocity_lower_error = massiv_all_dictionary['log_velocity_lower_error']
massiv_all_v_tot = massiv_all_dictionary['v_tot']
massiv_all_v_tot_lower_error = massiv_all_dictionary['v_tot_lower_error']
massiv_all_v_tot_upper_error = massiv_all_dictionary['v_tot_upper_error']
massiv_all_log_v_tot = massiv_all_dictionary['log_v_tot']
massiv_all_log_v_tot_lower_error = massiv_all_dictionary['log_v_tot_lower_error']
massiv_all_log_v_tot_upper_error = massiv_all_dictionary['log_v_tot_upper_error']
massiv_all_log_v_tot_av_error = massiv_all_dictionary['log_v_tot_av_error']
massiv_all_ang_mom = massiv_all_dictionary['ang_mom']
massiv_all_ang_mom_lower_error = massiv_all_dictionary['ang_mom_lower_error']
massiv_all_ang_mom_upper_error = massiv_all_dictionary['ang_mom_upper_error']
massiv_all_log_ang_mom = massiv_all_dictionary['log_ang_mom']
massiv_all_log_ang_mom_lower_error = massiv_all_dictionary['log_ang_mom_lower_error']
massiv_all_log_ang_mom_upper_error = massiv_all_dictionary['log_ang_mom_upper_error']
massiv_all_log_ang_mom_av_error = massiv_all_dictionary['log_ang_mom_av_error']
massiv_all_ang_mom_tot = massiv_all_dictionary['ang_mom_tot']
massiv_all_ang_mom_tot_lower_error = massiv_all_dictionary['ang_mom_tot_lower_error']
massiv_all_ang_mom_tot_upper_error = massiv_all_dictionary['ang_mom_tot_upper_error']
massiv_all_log_ang_mom_tot = massiv_all_dictionary['log_ang_mom_tot']
massiv_all_log_ang_mom_tot_lower_error = massiv_all_dictionary['log_ang_mom_tot_lower_error']
massiv_all_log_ang_mom_tot_upper_error = massiv_all_dictionary['log_ang_mom_tot_upper_error']
massiv_all_log_ang_mom_tot_av_error = massiv_all_dictionary['log_ang_mom_tot_av_error']
massiv_all_log_velocity_mass_fit_params = massiv_all_dictionary['log_velocity_mass_fit_params']
print 'massiv_all velocity vs. mass fit params: %s' % massiv_all_log_velocity_mass_fit_params
massiv_all_log_velocity_mass_fit_line = massiv_all_dictionary['log_velocity_mass_fit_line']
massiv_all_log_velocity_mass_fit_line_lower = massiv_all_dictionary['log_velocity_mass_fit_line_lower']
massiv_all_log_velocity_mass_fit_line_upper = massiv_all_dictionary['log_velocity_mass_fit_line_upper']
massiv_all_log_v_tot_mass_fit_params = massiv_all_dictionary['log_v_tot_mass_fit_params']
print 'massiv_all v_tot vs. mass fit params: %s' % massiv_all_log_v_tot_mass_fit_params
massiv_all_log_v_tot_mass_fit_line = massiv_all_dictionary['log_v_tot_mass_fit_line']
massiv_all_log_v_tot_mass_fit_line_lower = massiv_all_dictionary['log_v_tot_mass_fit_line_lower']
massiv_all_log_v_tot_mass_fit_line_upper = massiv_all_dictionary['log_v_tot_mass_fit_line_upper']
massiv_all_log_ang_mom_mass_fit_params = massiv_all_dictionary['log_ang_mom_mass_fit_params']
print 'massiv_all ang_mom vs. mass fit params: %s' % massiv_all_log_ang_mom_mass_fit_params
massiv_all_log_ang_mom_mass_fit_line = massiv_all_dictionary['log_ang_mom_mass_fit_line']
massiv_all_log_ang_mom_mass_fit_line_lower = massiv_all_dictionary['log_ang_mom_mass_fit_line_lower']
massiv_all_log_ang_mom_mass_fit_line_upper = massiv_all_dictionary['log_ang_mom_mass_fit_line_upper']
massiv_all_log_ang_mom_tot_mass_fit_params = massiv_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'massiv_all ang_mom_tot vs. mass fit params: %s' % massiv_all_log_ang_mom_tot_mass_fit_params
massiv_all_log_ang_mom_tot_mass_fit_line = massiv_all_dictionary['log_ang_mom_tot_mass_fit_line']
massiv_all_log_ang_mom_tot_mass_fit_line_lower = massiv_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
massiv_all_log_ang_mom_tot_mass_fit_line_upper = massiv_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
massiv_all_median_dynamical_mass_ratio = massiv_all_dictionary['median_dynamical_mass_ratio']
massiv_all_median_dynamical_mass_ratio_lower_error = massiv_all_dictionary['median_dynamical_mass_ratio_lower_error']
massiv_all_median_dynamical_mass_ratio_upper_error = massiv_all_dictionary['median_dynamical_mass_ratio_upper_error']
massiv_all_median_dynamical_mass_with_sigma_ratio = massiv_all_dictionary['median_dynamical_mass_with_sigma_ratio']
massiv_all_median_dynamical_mass_with_sigma_ratio_lower_error = massiv_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
massiv_all_median_dynamical_mass_with_sigma_ratio_upper_error = massiv_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR MASSIV ROT

massiv_rot_velocity = table_massiv_rot['Vmax']
massiv_rot_velocity_upper_error = table_massiv_rot['vmax_error']
massiv_rot_velocity_lower_error = table_massiv_rot['vmax_error']
massiv_rot_mass = table_massiv_rot['M_star']
massiv_rot_radius = table_massiv_rot['Re']
massiv_rot_sigma = table_massiv_rot['sigma']
massiv_rot_sigma_upper_error = table_massiv_rot['sigma_error']
massiv_rot_sigma_lower_error = table_massiv_rot['sigma_error']

massiv_rot_dictionary = return_betas(massiv_rot_mass,
                                    massiv_rot_velocity,
                                    massiv_rot_velocity_upper_error,
                                    massiv_rot_velocity_lower_error,
                                    massiv_rot_sigma,
                                    massiv_rot_sigma_upper_error,
                                    massiv_rot_sigma_lower_error,
                                    massiv_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = massiv_rot_dictionary['fit_mass']
massiv_rot_log_velocity = massiv_rot_dictionary['log_velocity']
massiv_rot_log_velocity_upper_error = massiv_rot_dictionary['log_velocity_upper_error']
massiv_rot_log_velocity_lower_error = massiv_rot_dictionary['log_velocity_lower_error']
massiv_rot_v_tot = massiv_rot_dictionary['v_tot']
massiv_rot_v_tot_lower_error = massiv_rot_dictionary['v_tot_lower_error']
massiv_rot_v_tot_upper_error = massiv_rot_dictionary['v_tot_upper_error']
massiv_rot_log_v_tot = massiv_rot_dictionary['log_v_tot']
massiv_rot_log_v_tot_lower_error = massiv_rot_dictionary['log_v_tot_lower_error']
massiv_rot_log_v_tot_upper_error = massiv_rot_dictionary['log_v_tot_upper_error']
massiv_rot_log_v_tot_av_error = massiv_rot_dictionary['log_v_tot_av_error']
massiv_rot_ang_mom = massiv_rot_dictionary['ang_mom']
massiv_rot_ang_mom_lower_error = massiv_rot_dictionary['ang_mom_lower_error']
massiv_rot_ang_mom_upper_error = massiv_rot_dictionary['ang_mom_upper_error']
massiv_rot_log_ang_mom = massiv_rot_dictionary['log_ang_mom']
massiv_rot_log_ang_mom_lower_error = massiv_rot_dictionary['log_ang_mom_lower_error']
massiv_rot_log_ang_mom_upper_error = massiv_rot_dictionary['log_ang_mom_upper_error']
massiv_rot_log_ang_mom_av_error = massiv_rot_dictionary['log_ang_mom_av_error']
massiv_rot_ang_mom_tot = massiv_rot_dictionary['ang_mom_tot']
massiv_rot_ang_mom_tot_lower_error = massiv_rot_dictionary['ang_mom_tot_lower_error']
massiv_rot_ang_mom_tot_upper_error = massiv_rot_dictionary['ang_mom_tot_upper_error']
massiv_rot_log_ang_mom_tot = massiv_rot_dictionary['log_ang_mom_tot']
massiv_rot_log_ang_mom_tot_lower_error = massiv_rot_dictionary['log_ang_mom_tot_lower_error']
massiv_rot_log_ang_mom_tot_upper_error = massiv_rot_dictionary['log_ang_mom_tot_upper_error']
massiv_rot_log_ang_mom_tot_av_error = massiv_rot_dictionary['log_ang_mom_tot_av_error']
massiv_rot_log_velocity_mass_fit_params = massiv_rot_dictionary['log_velocity_mass_fit_params']
print 'massiv_rot velocity vs. mass fit params: %s' % massiv_rot_log_velocity_mass_fit_params
massiv_rot_log_velocity_mass_fit_line = massiv_rot_dictionary['log_velocity_mass_fit_line']
massiv_rot_log_velocity_mass_fit_line_lower = massiv_rot_dictionary['log_velocity_mass_fit_line_lower']
massiv_rot_log_velocity_mass_fit_line_upper = massiv_rot_dictionary['log_velocity_mass_fit_line_upper']
massiv_rot_log_v_tot_mass_fit_params = massiv_rot_dictionary['log_v_tot_mass_fit_params']
print 'massiv_rot v_tot vs. mass fit params: %s' % massiv_rot_log_v_tot_mass_fit_params
massiv_rot_log_v_tot_mass_fit_line = massiv_rot_dictionary['log_v_tot_mass_fit_line']
massiv_rot_log_v_tot_mass_fit_line_lower = massiv_rot_dictionary['log_v_tot_mass_fit_line_lower']
massiv_rot_log_v_tot_mass_fit_line_upper = massiv_rot_dictionary['log_v_tot_mass_fit_line_upper']
massiv_rot_log_ang_mom_mass_fit_params = massiv_rot_dictionary['log_ang_mom_mass_fit_params']
print 'massiv_rot ang_mom vs. mass fit params: %s' % massiv_rot_log_ang_mom_mass_fit_params
massiv_rot_log_ang_mom_mass_fit_line = massiv_rot_dictionary['log_ang_mom_mass_fit_line']
massiv_rot_log_ang_mom_mass_fit_line_lower = massiv_rot_dictionary['log_ang_mom_mass_fit_line_lower']
massiv_rot_log_ang_mom_mass_fit_line_upper = massiv_rot_dictionary['log_ang_mom_mass_fit_line_upper']
massiv_rot_log_ang_mom_tot_mass_fit_params = massiv_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'massiv_rot ang_mom_tot vs. mass fit params: %s' % massiv_rot_log_ang_mom_tot_mass_fit_params
massiv_rot_log_ang_mom_tot_mass_fit_line = massiv_rot_dictionary['log_ang_mom_tot_mass_fit_line']
massiv_rot_log_ang_mom_tot_mass_fit_line_lower = massiv_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
massiv_rot_log_ang_mom_tot_mass_fit_line_upper = massiv_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
massiv_rot_median_dynamical_mass_ratio = massiv_rot_dictionary['median_dynamical_mass_ratio']
massiv_rot_median_dynamical_mass_ratio_lower_error = massiv_rot_dictionary['median_dynamical_mass_ratio_lower_error']
massiv_rot_median_dynamical_mass_ratio_upper_error = massiv_rot_dictionary['median_dynamical_mass_ratio_upper_error']
massiv_rot_median_dynamical_mass_with_sigma_ratio = massiv_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
massiv_rot_median_dynamical_mass_with_sigma_ratio_lower_error = massiv_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
massiv_rot_median_dynamical_mass_with_sigma_ratio_upper_error = massiv_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR MASSIV DISP

massiv_disp_velocity = table_massiv_disp['Vmax']
massiv_disp_velocity_upper_error = table_massiv_disp['vmax_error']
massiv_disp_velocity_lower_error = table_massiv_disp['vmax_error']
massiv_disp_mass = table_massiv_disp['M_star']
massiv_disp_radius = table_massiv_disp['Re']
massiv_disp_sigma = table_massiv_disp['sigma']
massiv_disp_sigma_upper_error = table_massiv_disp['sigma_error']
massiv_disp_sigma_lower_error = table_massiv_disp['sigma_error']

massiv_disp_dictionary = return_betas(massiv_disp_mass,
                                    massiv_disp_velocity,
                                    massiv_disp_velocity_upper_error,
                                    massiv_disp_velocity_lower_error,
                                    massiv_disp_sigma,
                                    massiv_disp_sigma_upper_error,
                                    massiv_disp_sigma_lower_error,
                                    massiv_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = massiv_disp_dictionary['fit_mass']
massiv_disp_log_velocity = massiv_disp_dictionary['log_velocity']
massiv_disp_log_velocity_upper_error = massiv_disp_dictionary['log_velocity_upper_error']
massiv_disp_log_velocity_lower_error = massiv_disp_dictionary['log_velocity_lower_error']
massiv_disp_v_tot = massiv_disp_dictionary['v_tot']
massiv_disp_v_tot_lower_error = massiv_disp_dictionary['v_tot_lower_error']
massiv_disp_v_tot_upper_error = massiv_disp_dictionary['v_tot_upper_error']
massiv_disp_log_v_tot = massiv_disp_dictionary['log_v_tot']
massiv_disp_log_v_tot_lower_error = massiv_disp_dictionary['log_v_tot_lower_error']
massiv_disp_log_v_tot_upper_error = massiv_disp_dictionary['log_v_tot_upper_error']
massiv_disp_log_v_tot_av_error = massiv_disp_dictionary['log_v_tot_av_error']
massiv_disp_ang_mom = massiv_disp_dictionary['ang_mom']
massiv_disp_ang_mom_lower_error = massiv_disp_dictionary['ang_mom_lower_error']
massiv_disp_ang_mom_upper_error = massiv_disp_dictionary['ang_mom_upper_error']
massiv_disp_log_ang_mom = massiv_disp_dictionary['log_ang_mom']
massiv_disp_log_ang_mom_lower_error = massiv_disp_dictionary['log_ang_mom_lower_error']
massiv_disp_log_ang_mom_upper_error = massiv_disp_dictionary['log_ang_mom_upper_error']
massiv_disp_log_ang_mom_av_error = massiv_disp_dictionary['log_ang_mom_av_error']
massiv_disp_ang_mom_tot = massiv_disp_dictionary['ang_mom_tot']
massiv_disp_ang_mom_tot_lower_error = massiv_disp_dictionary['ang_mom_tot_lower_error']
massiv_disp_ang_mom_tot_upper_error = massiv_disp_dictionary['ang_mom_tot_upper_error']
massiv_disp_log_ang_mom_tot = massiv_disp_dictionary['log_ang_mom_tot']
massiv_disp_log_ang_mom_tot_lower_error = massiv_disp_dictionary['log_ang_mom_tot_lower_error']
massiv_disp_log_ang_mom_tot_upper_error = massiv_disp_dictionary['log_ang_mom_tot_upper_error']
massiv_disp_log_ang_mom_tot_av_error = massiv_disp_dictionary['log_ang_mom_tot_av_error']
massiv_disp_log_velocity_mass_fit_params = massiv_disp_dictionary['log_velocity_mass_fit_params']
print 'massiv_disp velocity vs. mass fit params: %s' % massiv_disp_log_velocity_mass_fit_params
massiv_disp_log_velocity_mass_fit_line = massiv_disp_dictionary['log_velocity_mass_fit_line']
massiv_disp_log_velocity_mass_fit_line_lower = massiv_disp_dictionary['log_velocity_mass_fit_line_lower']
massiv_disp_log_velocity_mass_fit_line_upper = massiv_disp_dictionary['log_velocity_mass_fit_line_upper']
massiv_disp_log_v_tot_mass_fit_params = massiv_disp_dictionary['log_v_tot_mass_fit_params']
print 'massiv_disp v_tot vs. mass fit params: %s' % massiv_disp_log_v_tot_mass_fit_params
massiv_disp_log_v_tot_mass_fit_line = massiv_disp_dictionary['log_v_tot_mass_fit_line']
massiv_disp_log_v_tot_mass_fit_line_lower = massiv_disp_dictionary['log_v_tot_mass_fit_line_lower']
massiv_disp_log_v_tot_mass_fit_line_upper = massiv_disp_dictionary['log_v_tot_mass_fit_line_upper']
massiv_disp_log_ang_mom_mass_fit_params = massiv_disp_dictionary['log_ang_mom_mass_fit_params']
print 'massiv_disp ang_mom vs. mass fit params: %s' % massiv_disp_log_ang_mom_mass_fit_params
massiv_disp_log_ang_mom_mass_fit_line = massiv_disp_dictionary['log_ang_mom_mass_fit_line']
massiv_disp_log_ang_mom_mass_fit_line_lower = massiv_disp_dictionary['log_ang_mom_mass_fit_line_lower']
massiv_disp_log_ang_mom_mass_fit_line_upper = massiv_disp_dictionary['log_ang_mom_mass_fit_line_upper']
massiv_disp_log_ang_mom_tot_mass_fit_params = massiv_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'massiv_disp ang_mom_tot vs. mass fit params: %s' % massiv_disp_log_ang_mom_tot_mass_fit_params
massiv_disp_log_ang_mom_tot_mass_fit_line = massiv_disp_dictionary['log_ang_mom_tot_mass_fit_line']
massiv_disp_log_ang_mom_tot_mass_fit_line_lower = massiv_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
massiv_disp_log_ang_mom_tot_mass_fit_line_upper = massiv_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
massiv_disp_median_dynamical_mass_ratio = massiv_disp_dictionary['median_dynamical_mass_ratio']
massiv_disp_median_dynamical_mass_ratio_lower_error = massiv_disp_dictionary['median_dynamical_mass_ratio_lower_error']
massiv_disp_median_dynamical_mass_ratio_upper_error = massiv_disp_dictionary['median_dynamical_mass_ratio_upper_error']
massiv_disp_median_dynamical_mass_with_sigma_ratio = massiv_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
massiv_disp_median_dynamical_mass_with_sigma_ratio_lower_error = massiv_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
massiv_disp_median_dynamical_mass_with_sigma_ratio_upper_error = massiv_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['massiv',
                                1.2,
                                massiv_all_log_velocity_mass_fit_params[4],
                                massiv_all_log_velocity_mass_fit_params[3],
                                massiv_all_log_velocity_mass_fit_params[0],
                                massiv_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_all_log_velocity_mass_fit_params[0] - massiv_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((massiv_all_log_velocity_mass_fit_params[2] - massiv_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                massiv_rot_log_velocity_mass_fit_params[4],
                                massiv_rot_log_velocity_mass_fit_params[3],
                                massiv_rot_log_velocity_mass_fit_params[0],
                                massiv_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_rot_log_velocity_mass_fit_params[0] - massiv_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((massiv_rot_log_velocity_mass_fit_params[2] - massiv_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                massiv_disp_log_velocity_mass_fit_params[4],
                                massiv_disp_log_velocity_mass_fit_params[3],
                                massiv_disp_log_velocity_mass_fit_params[0],
                                massiv_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_disp_log_velocity_mass_fit_params[0] - massiv_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((massiv_disp_log_velocity_mass_fit_params[2] - massiv_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['massiv',
                                1.2,
                                massiv_all_log_v_tot_mass_fit_params[4],
                                massiv_all_log_v_tot_mass_fit_params[3],
                                massiv_all_log_v_tot_mass_fit_params[0],
                                massiv_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_all_log_v_tot_mass_fit_params[0] - massiv_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((massiv_all_log_v_tot_mass_fit_params[2] - massiv_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                massiv_rot_log_v_tot_mass_fit_params[4],
                                massiv_rot_log_v_tot_mass_fit_params[3],
                                massiv_rot_log_v_tot_mass_fit_params[0],
                                massiv_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_rot_log_v_tot_mass_fit_params[0] - massiv_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((massiv_rot_log_v_tot_mass_fit_params[2] - massiv_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                massiv_disp_log_v_tot_mass_fit_params[4],
                                massiv_disp_log_v_tot_mass_fit_params[3],
                                massiv_disp_log_v_tot_mass_fit_params[0],
                                massiv_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((massiv_disp_log_v_tot_mass_fit_params[0] - massiv_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((massiv_disp_log_v_tot_mass_fit_params[2] - massiv_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['massiv',
                                1.2,
                                massiv_all_log_ang_mom_mass_fit_params[4],
                                massiv_all_log_ang_mom_mass_fit_params[3],
                                massiv_all_log_ang_mom_mass_fit_params[0],
                                massiv_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_all_log_ang_mom_mass_fit_params[0] - massiv_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((massiv_all_log_ang_mom_mass_fit_params[2] - massiv_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                massiv_rot_log_ang_mom_mass_fit_params[4],
                                massiv_rot_log_ang_mom_mass_fit_params[3],
                                massiv_rot_log_ang_mom_mass_fit_params[0],
                                massiv_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_rot_log_ang_mom_mass_fit_params[0] - massiv_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((massiv_rot_log_ang_mom_mass_fit_params[2] - massiv_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                massiv_disp_log_ang_mom_mass_fit_params[4],
                                massiv_disp_log_ang_mom_mass_fit_params[3],
                                massiv_disp_log_ang_mom_mass_fit_params[0],
                                massiv_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_disp_log_ang_mom_mass_fit_params[0] - massiv_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((massiv_disp_log_ang_mom_mass_fit_params[2] - massiv_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['massiv',
                                1.2,
                                massiv_all_log_ang_mom_tot_mass_fit_params[4],
                                massiv_all_log_ang_mom_tot_mass_fit_params[3],
                                massiv_all_log_ang_mom_tot_mass_fit_params[0],
                                massiv_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_all_log_ang_mom_tot_mass_fit_params[0] - massiv_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((massiv_all_log_ang_mom_tot_mass_fit_params[2] - massiv_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                massiv_rot_log_ang_mom_tot_mass_fit_params[4],
                                massiv_rot_log_ang_mom_tot_mass_fit_params[3],
                                massiv_rot_log_ang_mom_tot_mass_fit_params[0],
                                massiv_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_rot_log_ang_mom_tot_mass_fit_params[0] - massiv_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((massiv_rot_log_ang_mom_tot_mass_fit_params[2] - massiv_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                massiv_disp_log_ang_mom_tot_mass_fit_params[4],
                                massiv_disp_log_ang_mom_tot_mass_fit_params[3],
                                massiv_disp_log_ang_mom_tot_mass_fit_params[0],
                                massiv_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((massiv_disp_log_ang_mom_tot_mass_fit_params[0] - massiv_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((massiv_disp_log_ang_mom_tot_mass_fit_params[2] - massiv_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['massiv',
                                1.2,
                                massiv_all_log_velocity_mass_fit_params[4],
                                massiv_all_median_dynamical_mass_ratio,
                                massiv_all_median_dynamical_mass_ratio_lower_error,
                                massiv_all_median_dynamical_mass_ratio_upper_error,
                                massiv_all_median_dynamical_mass_with_sigma_ratio,
                                massiv_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                massiv_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                massiv_rot_log_velocity_mass_fit_params[4],
                                massiv_rot_median_dynamical_mass_ratio,
                                massiv_rot_median_dynamical_mass_ratio_lower_error,
                                massiv_rot_median_dynamical_mass_ratio_upper_error,
                                massiv_rot_median_dynamical_mass_with_sigma_ratio,
                                massiv_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                massiv_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                massiv_disp_log_velocity_mass_fit_params[4],
                                massiv_disp_median_dynamical_mass_ratio,
                                massiv_disp_median_dynamical_mass_ratio_lower_error,
                                massiv_disp_median_dynamical_mass_ratio_upper_error,
                                massiv_disp_median_dynamical_mass_with_sigma_ratio,
                                massiv_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                massiv_disp_median_dynamical_mass_with_sigma_ratio_upper_error])


# DYNAMICAL PROPERTIES FOR AMAZE CLEAN (There is only all and rot here)

amaze_all_velocity = table_amaze_all['vmax']
amaze_all_velocity_upper_error = table_amaze_all['vmax_upper_error']
amaze_all_velocity_lower_error = table_amaze_all['vmax_lower_error']
amaze_all_mass = table_amaze_all['Mstar']
amaze_all_radius = table_amaze_all['r_half']
amaze_all_sigma = table_amaze_all['sigma']
amaze_all_sigma_upper_error = table_amaze_all['sigma_upper_error']
amaze_all_sigma_lower_error = table_amaze_all['sigma_lower_error']

amaze_all_dictionary = return_betas(amaze_all_mass,
                                    amaze_all_velocity,
                                    amaze_all_velocity_upper_error,
                                    amaze_all_velocity_lower_error,
                                    amaze_all_sigma,
                                    amaze_all_sigma_upper_error,
                                    amaze_all_sigma_lower_error,
                                    amaze_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = amaze_all_dictionary['fit_mass']
amaze_all_log_velocity = amaze_all_dictionary['log_velocity']
amaze_all_log_velocity_upper_error = amaze_all_dictionary['log_velocity_upper_error']
amaze_all_log_velocity_lower_error = amaze_all_dictionary['log_velocity_lower_error']
amaze_all_v_tot = amaze_all_dictionary['v_tot']
amaze_all_v_tot_lower_error = amaze_all_dictionary['v_tot_lower_error']
amaze_all_v_tot_upper_error = amaze_all_dictionary['v_tot_upper_error']
amaze_all_log_v_tot = amaze_all_dictionary['log_v_tot']
amaze_all_log_v_tot_lower_error = amaze_all_dictionary['log_v_tot_lower_error']
amaze_all_log_v_tot_upper_error = amaze_all_dictionary['log_v_tot_upper_error']
amaze_all_log_v_tot_av_error = amaze_all_dictionary['log_v_tot_av_error']
amaze_all_ang_mom = amaze_all_dictionary['ang_mom']
amaze_all_ang_mom_lower_error = amaze_all_dictionary['ang_mom_lower_error']
amaze_all_ang_mom_upper_error = amaze_all_dictionary['ang_mom_upper_error']
amaze_all_log_ang_mom = amaze_all_dictionary['log_ang_mom']
amaze_all_log_ang_mom_lower_error = amaze_all_dictionary['log_ang_mom_lower_error']
amaze_all_log_ang_mom_upper_error = amaze_all_dictionary['log_ang_mom_upper_error']
amaze_all_log_ang_mom_av_error = amaze_all_dictionary['log_ang_mom_av_error']
amaze_all_ang_mom_tot = amaze_all_dictionary['ang_mom_tot']
amaze_all_ang_mom_tot_lower_error = amaze_all_dictionary['ang_mom_tot_lower_error']
amaze_all_ang_mom_tot_upper_error = amaze_all_dictionary['ang_mom_tot_upper_error']
amaze_all_log_ang_mom_tot = amaze_all_dictionary['log_ang_mom_tot']
amaze_all_log_ang_mom_tot_lower_error = amaze_all_dictionary['log_ang_mom_tot_lower_error']
amaze_all_log_ang_mom_tot_upper_error = amaze_all_dictionary['log_ang_mom_tot_upper_error']
amaze_all_log_ang_mom_tot_av_error = amaze_all_dictionary['log_ang_mom_tot_av_error']
amaze_all_log_velocity_mass_fit_params = amaze_all_dictionary['log_velocity_mass_fit_params']
print 'amaze_all velocity vs. mass fit params: %s' % amaze_all_log_velocity_mass_fit_params
amaze_all_log_velocity_mass_fit_line = amaze_all_dictionary['log_velocity_mass_fit_line']
amaze_all_log_velocity_mass_fit_line_lower = amaze_all_dictionary['log_velocity_mass_fit_line_lower']
amaze_all_log_velocity_mass_fit_line_upper = amaze_all_dictionary['log_velocity_mass_fit_line_upper']
amaze_all_log_v_tot_mass_fit_params = amaze_all_dictionary['log_v_tot_mass_fit_params']
print 'amaze_all v_tot vs. mass fit params: %s' % amaze_all_log_v_tot_mass_fit_params
amaze_all_log_v_tot_mass_fit_line = amaze_all_dictionary['log_v_tot_mass_fit_line']
amaze_all_log_v_tot_mass_fit_line_lower = amaze_all_dictionary['log_v_tot_mass_fit_line_lower']
amaze_all_log_v_tot_mass_fit_line_upper = amaze_all_dictionary['log_v_tot_mass_fit_line_upper']
amaze_all_log_ang_mom_mass_fit_params = amaze_all_dictionary['log_ang_mom_mass_fit_params']
print 'amaze_all ang_mom vs. mass fit params: %s' % amaze_all_log_ang_mom_mass_fit_params
amaze_all_log_ang_mom_mass_fit_line = amaze_all_dictionary['log_ang_mom_mass_fit_line']
amaze_all_log_ang_mom_mass_fit_line_lower = amaze_all_dictionary['log_ang_mom_mass_fit_line_lower']
amaze_all_log_ang_mom_mass_fit_line_upper = amaze_all_dictionary['log_ang_mom_mass_fit_line_upper']
amaze_all_log_ang_mom_tot_mass_fit_params = amaze_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'amaze_all ang_mom_tot vs. mass fit params: %s' % amaze_all_log_ang_mom_tot_mass_fit_params
amaze_all_log_ang_mom_tot_mass_fit_line = amaze_all_dictionary['log_ang_mom_tot_mass_fit_line']
amaze_all_log_ang_mom_tot_mass_fit_line_lower = amaze_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
amaze_all_log_ang_mom_tot_mass_fit_line_upper = amaze_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
amaze_all_median_dynamical_mass_ratio = amaze_all_dictionary['median_dynamical_mass_ratio']
amaze_all_median_dynamical_mass_ratio_lower_error = amaze_all_dictionary['median_dynamical_mass_ratio_lower_error']
amaze_all_median_dynamical_mass_ratio_upper_error = amaze_all_dictionary['median_dynamical_mass_ratio_upper_error']
amaze_all_median_dynamical_mass_with_sigma_ratio = amaze_all_dictionary['median_dynamical_mass_with_sigma_ratio']
amaze_all_median_dynamical_mass_with_sigma_ratio_lower_error = amaze_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
amaze_all_median_dynamical_mass_with_sigma_ratio_upper_error = amaze_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['amaze',
                                3.0,
                                amaze_all_log_velocity_mass_fit_params[4],
                                amaze_all_log_velocity_mass_fit_params[3],
                                amaze_all_log_velocity_mass_fit_params[0],
                                amaze_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((amaze_all_log_velocity_mass_fit_params[0] - amaze_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((amaze_all_log_velocity_mass_fit_params[2] - amaze_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
v_tot_evolution_list.append(['amaze',
                                3.0,
                                amaze_all_log_v_tot_mass_fit_params[4],
                                amaze_all_log_v_tot_mass_fit_params[3],
                                amaze_all_log_v_tot_mass_fit_params[0],
                                amaze_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((amaze_all_log_v_tot_mass_fit_params[0] - amaze_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((amaze_all_log_v_tot_mass_fit_params[2] - amaze_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_evolution_list.append(['amaze',
                                3.0,
                                amaze_all_log_ang_mom_mass_fit_params[4],
                                amaze_all_log_ang_mom_mass_fit_params[3],
                                amaze_all_log_ang_mom_mass_fit_params[0],
                                amaze_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((amaze_all_log_ang_mom_mass_fit_params[0] - amaze_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((amaze_all_log_ang_mom_mass_fit_params[2] - amaze_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
ang_mom_tot_evolution_list.append(['amaze',
                                3.0,
                                amaze_all_log_ang_mom_tot_mass_fit_params[4],
                                amaze_all_log_ang_mom_tot_mass_fit_params[3],
                                amaze_all_log_ang_mom_tot_mass_fit_params[0],
                                amaze_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((amaze_all_log_ang_mom_tot_mass_fit_params[0] - amaze_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((amaze_all_log_ang_mom_tot_mass_fit_params[2] - amaze_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])
dyn_mass_evolution_list.append(['amaze',
                                3.0,
                                amaze_all_log_velocity_mass_fit_params[4],
                                amaze_all_median_dynamical_mass_ratio,
                                amaze_all_median_dynamical_mass_ratio_lower_error,
                                amaze_all_median_dynamical_mass_ratio_upper_error,
                                amaze_all_median_dynamical_mass_with_sigma_ratio,
                                amaze_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                amaze_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0,
                                0])

# DYNAMICAL PROPERTIES FOR KDS ALL

kds_all_mass = table_kds_all['Mstar']
kds_all_indices = np.where(kds_all_mass > 9.0)[0]
kds_all_mass = table_kds_all['Mstar'][kds_all_indices]
kds_all_velocity = table_kds_all['2d_intrinsic_Vmax_3Rd'][kds_all_indices]
kds_all_velocity_upper_error = table_kds_all['2d_intrinsic_Vmax_3Rd_upper_error'][kds_all_indices]
kds_all_velocity_lower_error = table_kds_all['2d_intrinsic_Vmax_3Rd_lower_error'][kds_all_indices]
kds_all_radius = table_kds_all['R_e(Kpc)'][kds_all_indices]
kds_all_sigma = table_kds_all['mean_intrinsic_sigma'][kds_all_indices]
kds_all_sigma_upper_error = table_kds_all['Sigma_upper_error'][kds_all_indices]
kds_all_sigma_lower_error = table_kds_all['Sigma_lower_error'][kds_all_indices]

# only fit those where mass is greater than 9.0


kds_all_dictionary = return_betas(kds_all_mass,
                                    kds_all_velocity,
                                    kds_all_velocity_upper_error,
                                    kds_all_velocity_lower_error,
                                    kds_all_sigma,
                                    kds_all_sigma_upper_error,
                                    kds_all_sigma_lower_error,
                                    kds_all_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kds_all_dictionary['fit_mass']
kds_all_log_velocity = kds_all_dictionary['log_velocity']
kds_all_log_velocity_upper_error = kds_all_dictionary['log_velocity_upper_error']
kds_all_log_velocity_lower_error = kds_all_dictionary['log_velocity_lower_error']
kds_all_v_tot = kds_all_dictionary['v_tot']
kds_all_v_tot_lower_error = kds_all_dictionary['v_tot_lower_error']
kds_all_v_tot_upper_error = kds_all_dictionary['v_tot_upper_error']
kds_all_log_v_tot = kds_all_dictionary['log_v_tot']
kds_all_log_v_tot_lower_error = kds_all_dictionary['log_v_tot_lower_error']
kds_all_log_v_tot_upper_error = kds_all_dictionary['log_v_tot_upper_error']
kds_all_log_v_tot_av_error = kds_all_dictionary['log_v_tot_av_error']
kds_all_ang_mom = kds_all_dictionary['ang_mom']
kds_all_ang_mom_lower_error = kds_all_dictionary['ang_mom_lower_error']
kds_all_ang_mom_upper_error = kds_all_dictionary['ang_mom_upper_error']
kds_all_log_ang_mom = kds_all_dictionary['log_ang_mom']
kds_all_log_ang_mom_lower_error = kds_all_dictionary['log_ang_mom_lower_error']
kds_all_log_ang_mom_upper_error = kds_all_dictionary['log_ang_mom_upper_error']
kds_all_log_ang_mom_av_error = kds_all_dictionary['log_ang_mom_av_error']
kds_all_ang_mom_tot = kds_all_dictionary['ang_mom_tot']
kds_all_ang_mom_tot_lower_error = kds_all_dictionary['ang_mom_tot_lower_error']
kds_all_ang_mom_tot_upper_error = kds_all_dictionary['ang_mom_tot_upper_error']
kds_all_log_ang_mom_tot = kds_all_dictionary['log_ang_mom_tot']
kds_all_log_ang_mom_tot_lower_error = kds_all_dictionary['log_ang_mom_tot_lower_error']
kds_all_log_ang_mom_tot_upper_error = kds_all_dictionary['log_ang_mom_tot_upper_error']
kds_all_log_ang_mom_tot_av_error = kds_all_dictionary['log_ang_mom_tot_av_error']
kds_all_log_velocity_mass_fit_params = kds_all_dictionary['log_velocity_mass_fit_params']
print 'kds_all velocity vs. mass fit params: %s' % kds_all_log_velocity_mass_fit_params
kds_all_log_velocity_mass_fit_line = kds_all_dictionary['log_velocity_mass_fit_line']
kds_all_log_velocity_mass_fit_line_lower = kds_all_dictionary['log_velocity_mass_fit_line_lower']
kds_all_log_velocity_mass_fit_line_upper = kds_all_dictionary['log_velocity_mass_fit_line_upper']
kds_all_log_v_tot_mass_fit_params = kds_all_dictionary['log_v_tot_mass_fit_params']
print 'kds_all v_tot vs. mass fit params: %s' % kds_all_log_v_tot_mass_fit_params
kds_all_log_v_tot_mass_fit_line = kds_all_dictionary['log_v_tot_mass_fit_line']
kds_all_log_v_tot_mass_fit_line_lower = kds_all_dictionary['log_v_tot_mass_fit_line_lower']
kds_all_log_v_tot_mass_fit_line_upper = kds_all_dictionary['log_v_tot_mass_fit_line_upper']
kds_all_log_ang_mom_mass_fit_params = kds_all_dictionary['log_ang_mom_mass_fit_params']
print 'kds_all ang_mom vs. mass fit params: %s' % kds_all_log_ang_mom_mass_fit_params
kds_all_log_ang_mom_mass_fit_line = kds_all_dictionary['log_ang_mom_mass_fit_line']
kds_all_log_ang_mom_mass_fit_line_lower = kds_all_dictionary['log_ang_mom_mass_fit_line_lower']
kds_all_log_ang_mom_mass_fit_line_upper = kds_all_dictionary['log_ang_mom_mass_fit_line_upper']
kds_all_log_ang_mom_tot_mass_fit_params = kds_all_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kds_all ang_mom_tot vs. mass fit params: %s' % kds_all_log_ang_mom_tot_mass_fit_params
kds_all_log_ang_mom_tot_mass_fit_line = kds_all_dictionary['log_ang_mom_tot_mass_fit_line']
kds_all_log_ang_mom_tot_mass_fit_line_lower = kds_all_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kds_all_log_ang_mom_tot_mass_fit_line_upper = kds_all_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kds_all_median_dynamical_mass_ratio = kds_all_dictionary['median_dynamical_mass_ratio']
kds_all_median_dynamical_mass_ratio_lower_error = kds_all_dictionary['median_dynamical_mass_ratio_lower_error']
kds_all_median_dynamical_mass_ratio_upper_error = kds_all_dictionary['median_dynamical_mass_ratio_upper_error']
kds_all_median_dynamical_mass_with_sigma_ratio = kds_all_dictionary['median_dynamical_mass_with_sigma_ratio']
kds_all_median_dynamical_mass_with_sigma_ratio_lower_error = kds_all_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kds_all_median_dynamical_mass_with_sigma_ratio_upper_error = kds_all_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR KDS ROT

kds_rot_mass = table_kds_rot['Mstar']
kds_rot_indices = np.where(kds_rot_mass > 9.0)[0]
kds_rot_mass = table_kds_rot['Mstar'][kds_rot_indices]
kds_rot_velocity = table_kds_rot['2d_intrinsic_Vmax_3Rd'][kds_rot_indices]
kds_rot_velocity_upper_error = table_kds_rot['2d_intrinsic_Vmax_3Rd_upper_error'][kds_rot_indices]
kds_rot_velocity_lower_error = table_kds_rot['2d_intrinsic_Vmax_3Rd_lower_error'][kds_rot_indices]
kds_rot_radius = table_kds_rot['R_e(Kpc)'][kds_rot_indices]
kds_rot_sigma = table_kds_rot['mean_intrinsic_sigma'][kds_rot_indices]
kds_rot_sigma_upper_error = table_kds_rot['Sigma_upper_error'][kds_rot_indices]
kds_rot_sigma_lower_error = table_kds_rot['Sigma_lower_error'][kds_rot_indices]

kds_rot_dictionary = return_betas(kds_rot_mass,
                                    kds_rot_velocity,
                                    kds_rot_velocity_upper_error,
                                    kds_rot_velocity_lower_error,
                                    kds_rot_sigma,
                                    kds_rot_sigma_upper_error,
                                    kds_rot_sigma_lower_error,
                                    kds_rot_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kds_rot_dictionary['fit_mass']
kds_rot_log_velocity = kds_rot_dictionary['log_velocity']
kds_rot_log_velocity_upper_error = kds_rot_dictionary['log_velocity_upper_error']
kds_rot_log_velocity_lower_error = kds_rot_dictionary['log_velocity_lower_error']
kds_rot_v_tot = kds_rot_dictionary['v_tot']
kds_rot_v_tot_lower_error = kds_rot_dictionary['v_tot_lower_error']
kds_rot_v_tot_upper_error = kds_rot_dictionary['v_tot_upper_error']
kds_rot_log_v_tot = kds_rot_dictionary['log_v_tot']
kds_rot_log_v_tot_lower_error = kds_rot_dictionary['log_v_tot_lower_error']
kds_rot_log_v_tot_upper_error = kds_rot_dictionary['log_v_tot_upper_error']
kds_rot_log_v_tot_av_error = kds_rot_dictionary['log_v_tot_av_error']
kds_rot_ang_mom = kds_rot_dictionary['ang_mom']
kds_rot_ang_mom_lower_error = kds_rot_dictionary['ang_mom_lower_error']
kds_rot_ang_mom_upper_error = kds_rot_dictionary['ang_mom_upper_error']
kds_rot_log_ang_mom = kds_rot_dictionary['log_ang_mom']
kds_rot_log_ang_mom_lower_error = kds_rot_dictionary['log_ang_mom_lower_error']
kds_rot_log_ang_mom_upper_error = kds_rot_dictionary['log_ang_mom_upper_error']
kds_rot_log_ang_mom_av_error = kds_rot_dictionary['log_ang_mom_av_error']
kds_rot_ang_mom_tot = kds_rot_dictionary['ang_mom_tot']
kds_rot_ang_mom_tot_lower_error = kds_rot_dictionary['ang_mom_tot_lower_error']
kds_rot_ang_mom_tot_upper_error = kds_rot_dictionary['ang_mom_tot_upper_error']
kds_rot_log_ang_mom_tot = kds_rot_dictionary['log_ang_mom_tot']
kds_rot_log_ang_mom_tot_lower_error = kds_rot_dictionary['log_ang_mom_tot_lower_error']
kds_rot_log_ang_mom_tot_upper_error = kds_rot_dictionary['log_ang_mom_tot_upper_error']
kds_rot_log_ang_mom_tot_av_error = kds_rot_dictionary['log_ang_mom_tot_av_error']
kds_rot_log_velocity_mass_fit_params = kds_rot_dictionary['log_velocity_mass_fit_params']
print 'kds_rot velocity vs. mass fit params: %s' % kds_rot_log_velocity_mass_fit_params
kds_rot_log_velocity_mass_fit_line = kds_rot_dictionary['log_velocity_mass_fit_line']
kds_rot_log_velocity_mass_fit_line_lower = kds_rot_dictionary['log_velocity_mass_fit_line_lower']
kds_rot_log_velocity_mass_fit_line_upper = kds_rot_dictionary['log_velocity_mass_fit_line_upper']
kds_rot_log_v_tot_mass_fit_params = kds_rot_dictionary['log_v_tot_mass_fit_params']
print 'kds_rot v_tot vs. mass fit params: %s' % kds_rot_log_v_tot_mass_fit_params
kds_rot_log_v_tot_mass_fit_line = kds_rot_dictionary['log_v_tot_mass_fit_line']
kds_rot_log_v_tot_mass_fit_line_lower = kds_rot_dictionary['log_v_tot_mass_fit_line_lower']
kds_rot_log_v_tot_mass_fit_line_upper = kds_rot_dictionary['log_v_tot_mass_fit_line_upper']
kds_rot_log_ang_mom_mass_fit_params = kds_rot_dictionary['log_ang_mom_mass_fit_params']
print 'kds_rot ang_mom vs. mass fit params: %s' % kds_rot_log_ang_mom_mass_fit_params
kds_rot_log_ang_mom_mass_fit_line = kds_rot_dictionary['log_ang_mom_mass_fit_line']
kds_rot_log_ang_mom_mass_fit_line_lower = kds_rot_dictionary['log_ang_mom_mass_fit_line_lower']
kds_rot_log_ang_mom_mass_fit_line_upper = kds_rot_dictionary['log_ang_mom_mass_fit_line_upper']
kds_rot_log_ang_mom_tot_mass_fit_params = kds_rot_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kds_rot ang_mom_tot vs. mass fit params: %s' % kds_rot_log_ang_mom_tot_mass_fit_params
kds_rot_log_ang_mom_tot_mass_fit_line = kds_rot_dictionary['log_ang_mom_tot_mass_fit_line']
kds_rot_log_ang_mom_tot_mass_fit_line_lower = kds_rot_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kds_rot_log_ang_mom_tot_mass_fit_line_upper = kds_rot_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kds_rot_median_dynamical_mass_ratio = kds_rot_dictionary['median_dynamical_mass_ratio']
kds_rot_median_dynamical_mass_ratio_lower_error = kds_rot_dictionary['median_dynamical_mass_ratio_lower_error']
kds_rot_median_dynamical_mass_ratio_upper_error = kds_rot_dictionary['median_dynamical_mass_ratio_upper_error']
kds_rot_median_dynamical_mass_with_sigma_ratio = kds_rot_dictionary['median_dynamical_mass_with_sigma_ratio']
kds_rot_median_dynamical_mass_with_sigma_ratio_lower_error = kds_rot_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kds_rot_median_dynamical_mass_with_sigma_ratio_upper_error = kds_rot_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']

# DYNAMICAL PROPERTIES FOR KDS DISP

kds_disp_mass = table_kds_disp['Mstar']
kds_disp_indices = np.where(kds_disp_mass > 9.0)[0]
kds_disp_mass = table_kds_disp['Mstar'][kds_disp_indices]
kds_disp_velocity = table_kds_disp['2d_intrinsic_Vmax_3Rd'][kds_disp_indices]
kds_disp_velocity_upper_error = table_kds_disp['2d_intrinsic_Vmax_3Rd_upper_error'][kds_disp_indices]
kds_disp_velocity_lower_error = table_kds_disp['2d_intrinsic_Vmax_3Rd_lower_error'][kds_disp_indices]
kds_disp_radius = table_kds_disp['R_e(Kpc)'][kds_disp_indices]
kds_disp_sigma = table_kds_disp['mean_intrinsic_sigma'][kds_disp_indices]
kds_disp_sigma_upper_error = table_kds_disp['Sigma_upper_error'][kds_disp_indices]
kds_disp_sigma_lower_error = table_kds_disp['Sigma_lower_error'][kds_disp_indices]

kds_disp_dictionary = return_betas(kds_disp_mass,
                                   kds_disp_velocity,
                                   kds_disp_velocity_upper_error,
                                   kds_disp_velocity_lower_error,
                                   kds_disp_sigma,
                                   kds_disp_sigma_upper_error,
                                   kds_disp_sigma_lower_error,
                                   kds_disp_radius)
# assign the dictionary outputs to keywords in standardised form
fit_mass = kds_disp_dictionary['fit_mass']
kds_disp_log_velocity = kds_disp_dictionary['log_velocity']
kds_disp_log_velocity_upper_error = kds_disp_dictionary['log_velocity_upper_error']
kds_disp_log_velocity_lower_error = kds_disp_dictionary['log_velocity_lower_error']
kds_disp_v_tot = kds_disp_dictionary['v_tot']
kds_disp_v_tot_lower_error = kds_disp_dictionary['v_tot_lower_error']
kds_disp_v_tot_upper_error = kds_disp_dictionary['v_tot_upper_error']
kds_disp_log_v_tot = kds_disp_dictionary['log_v_tot']
kds_disp_log_v_tot_lower_error = kds_disp_dictionary['log_v_tot_lower_error']
kds_disp_log_v_tot_upper_error = kds_disp_dictionary['log_v_tot_upper_error']
kds_disp_log_v_tot_av_error = kds_disp_dictionary['log_v_tot_av_error']
kds_disp_ang_mom = kds_disp_dictionary['ang_mom']
kds_disp_ang_mom_lower_error = kds_disp_dictionary['ang_mom_lower_error']
kds_disp_ang_mom_upper_error = kds_disp_dictionary['ang_mom_upper_error']
kds_disp_log_ang_mom = kds_disp_dictionary['log_ang_mom']
kds_disp_log_ang_mom_lower_error = kds_disp_dictionary['log_ang_mom_lower_error']
kds_disp_log_ang_mom_upper_error = kds_disp_dictionary['log_ang_mom_upper_error']
kds_disp_log_ang_mom_av_error = kds_disp_dictionary['log_ang_mom_av_error']
kds_disp_ang_mom_tot = kds_disp_dictionary['ang_mom_tot']
kds_disp_ang_mom_tot_lower_error = kds_disp_dictionary['ang_mom_tot_lower_error']
kds_disp_ang_mom_tot_upper_error = kds_disp_dictionary['ang_mom_tot_upper_error']
kds_disp_log_ang_mom_tot = kds_disp_dictionary['log_ang_mom_tot']
kds_disp_log_ang_mom_tot_lower_error = kds_disp_dictionary['log_ang_mom_tot_lower_error']
kds_disp_log_ang_mom_tot_upper_error = kds_disp_dictionary['log_ang_mom_tot_upper_error']
kds_disp_log_ang_mom_tot_av_error = kds_disp_dictionary['log_ang_mom_tot_av_error']
kds_disp_log_velocity_mass_fit_params = kds_disp_dictionary['log_velocity_mass_fit_params']
print 'kds_disp velocity vs. mass fit params: %s' % kds_disp_log_velocity_mass_fit_params
kds_disp_log_velocity_mass_fit_line = kds_disp_dictionary['log_velocity_mass_fit_line']
kds_disp_log_velocity_mass_fit_line_lower = kds_disp_dictionary['log_velocity_mass_fit_line_lower']
kds_disp_log_velocity_mass_fit_line_upper = kds_disp_dictionary['log_velocity_mass_fit_line_upper']
kds_disp_log_v_tot_mass_fit_params = kds_disp_dictionary['log_v_tot_mass_fit_params']
print 'kds_disp v_tot vs. mass fit params: %s' % kds_disp_log_v_tot_mass_fit_params
kds_disp_log_v_tot_mass_fit_line = kds_disp_dictionary['log_v_tot_mass_fit_line']
kds_disp_log_v_tot_mass_fit_line_lower = kds_disp_dictionary['log_v_tot_mass_fit_line_lower']
kds_disp_log_v_tot_mass_fit_line_upper = kds_disp_dictionary['log_v_tot_mass_fit_line_upper']
kds_disp_log_ang_mom_mass_fit_params = kds_disp_dictionary['log_ang_mom_mass_fit_params']
print 'kds_disp ang_mom vs. mass fit params: %s' % kds_disp_log_ang_mom_mass_fit_params
kds_disp_log_ang_mom_mass_fit_line = kds_disp_dictionary['log_ang_mom_mass_fit_line']
kds_disp_log_ang_mom_mass_fit_line_lower = kds_disp_dictionary['log_ang_mom_mass_fit_line_lower']
kds_disp_log_ang_mom_mass_fit_line_upper = kds_disp_dictionary['log_ang_mom_mass_fit_line_upper']
kds_disp_log_ang_mom_tot_mass_fit_params = kds_disp_dictionary['log_ang_mom_tot_mass_fit_params']
print 'kds_disp ang_mom_tot vs. mass fit params: %s' % kds_disp_log_ang_mom_tot_mass_fit_params
kds_disp_log_ang_mom_tot_mass_fit_line = kds_disp_dictionary['log_ang_mom_tot_mass_fit_line']
kds_disp_log_ang_mom_tot_mass_fit_line_lower = kds_disp_dictionary['log_ang_mom_tot_mass_fit_line_lower']
kds_disp_log_ang_mom_tot_mass_fit_line_upper = kds_disp_dictionary['log_ang_mom_tot_mass_fit_line_upper']
kds_disp_median_dynamical_mass_ratio = kds_disp_dictionary['median_dynamical_mass_ratio']
kds_disp_median_dynamical_mass_ratio_lower_error = kds_disp_dictionary['median_dynamical_mass_ratio_lower_error']
kds_disp_median_dynamical_mass_ratio_upper_error = kds_disp_dictionary['median_dynamical_mass_ratio_upper_error']
kds_disp_median_dynamical_mass_with_sigma_ratio = kds_disp_dictionary['median_dynamical_mass_with_sigma_ratio']
kds_disp_median_dynamical_mass_with_sigma_ratio_lower_error = kds_disp_dictionary['median_dynamical_mass_with_sigma_ratio_lower_error']
kds_disp_median_dynamical_mass_with_sigma_ratio_upper_error = kds_disp_dictionary['median_dynamical_mass_with_sigma_ratio_upper_error']
# populate the list for velocity and append this to the evolution list
velocity_evolution_list.append(['kds',
                                3.5,
                                kds_all_log_velocity_mass_fit_params[4],
                                kds_all_log_velocity_mass_fit_params[3],
                                kds_all_log_velocity_mass_fit_params[0],
                                kds_all_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_all_log_velocity_mass_fit_params[0] - kds_all_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kds_all_log_velocity_mass_fit_params[2] - kds_all_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                kds_rot_log_velocity_mass_fit_params[4],
                                kds_rot_log_velocity_mass_fit_params[3],
                                kds_rot_log_velocity_mass_fit_params[0],
                                kds_rot_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_rot_log_velocity_mass_fit_params[0] - kds_rot_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kds_rot_log_velocity_mass_fit_params[2] - kds_rot_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),
                                kds_disp_log_velocity_mass_fit_params[4],
                                kds_disp_log_velocity_mass_fit_params[3],
                                kds_disp_log_velocity_mass_fit_params[0],
                                kds_disp_log_velocity_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_disp_log_velocity_mass_fit_params[0] - kds_disp_log_velocity_mass_fit_params[1])**2 + velocity_beta_zero_error**2),
                                np.sqrt((kds_disp_log_velocity_mass_fit_params[2] - kds_disp_log_velocity_mass_fit_params[0])**2 + velocity_beta_zero_error**2),])
v_tot_evolution_list.append(['kds',
                                3.5,
                                kds_all_log_v_tot_mass_fit_params[4],
                                kds_all_log_v_tot_mass_fit_params[3],
                                kds_all_log_v_tot_mass_fit_params[0],
                                kds_all_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_all_log_v_tot_mass_fit_params[0] - kds_all_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kds_all_log_v_tot_mass_fit_params[2] - kds_all_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                kds_rot_log_v_tot_mass_fit_params[4],
                                kds_rot_log_v_tot_mass_fit_params[3],
                                kds_rot_log_v_tot_mass_fit_params[0],
                                kds_rot_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_rot_log_v_tot_mass_fit_params[0] - kds_rot_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kds_rot_log_v_tot_mass_fit_params[2] - kds_rot_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),
                                kds_disp_log_v_tot_mass_fit_params[4],
                                kds_disp_log_v_tot_mass_fit_params[3],
                                kds_disp_log_v_tot_mass_fit_params[0],
                                kds_disp_log_v_tot_mass_fit_params[0] - velocity_beta_zero,
                                np.sqrt((kds_disp_log_v_tot_mass_fit_params[0] - kds_disp_log_v_tot_mass_fit_params[1])**2 + v_tot_beta_zero_error**2),
                                np.sqrt((kds_disp_log_v_tot_mass_fit_params[2] - kds_disp_log_v_tot_mass_fit_params[0])**2 + v_tot_beta_zero_error**2),])
ang_mom_evolution_list.append(['kds',
                                3.5,
                                kds_all_log_ang_mom_mass_fit_params[4],
                                kds_all_log_ang_mom_mass_fit_params[3],
                                kds_all_log_ang_mom_mass_fit_params[0],
                                kds_all_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_all_log_ang_mom_mass_fit_params[0] - kds_all_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kds_all_log_ang_mom_mass_fit_params[2] - kds_all_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                kds_rot_log_ang_mom_mass_fit_params[4],
                                kds_rot_log_ang_mom_mass_fit_params[3],
                                kds_rot_log_ang_mom_mass_fit_params[0],
                                kds_rot_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_rot_log_ang_mom_mass_fit_params[0] - kds_rot_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kds_rot_log_ang_mom_mass_fit_params[2] - kds_rot_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),
                                kds_disp_log_ang_mom_mass_fit_params[4],
                                kds_disp_log_ang_mom_mass_fit_params[3],
                                kds_disp_log_ang_mom_mass_fit_params[0],
                                kds_disp_log_ang_mom_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_disp_log_ang_mom_mass_fit_params[0] - kds_disp_log_ang_mom_mass_fit_params[1])**2 + ang_mom_beta_zero_error**2),
                                np.sqrt((kds_disp_log_ang_mom_mass_fit_params[2] - kds_disp_log_ang_mom_mass_fit_params[0])**2 + ang_mom_beta_zero_error**2),])
ang_mom_tot_evolution_list.append(['kds',
                                3.5,
                                kds_all_log_ang_mom_tot_mass_fit_params[4],
                                kds_all_log_ang_mom_tot_mass_fit_params[3],
                                kds_all_log_ang_mom_tot_mass_fit_params[0],
                                kds_all_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_all_log_ang_mom_tot_mass_fit_params[0] - kds_all_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kds_all_log_ang_mom_tot_mass_fit_params[2] - kds_all_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                kds_rot_log_ang_mom_tot_mass_fit_params[4],
                                kds_rot_log_ang_mom_tot_mass_fit_params[3],
                                kds_rot_log_ang_mom_tot_mass_fit_params[0],
                                kds_rot_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_rot_log_ang_mom_tot_mass_fit_params[0] - kds_rot_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kds_rot_log_ang_mom_tot_mass_fit_params[2] - kds_rot_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),
                                kds_disp_log_ang_mom_tot_mass_fit_params[4],
                                kds_disp_log_ang_mom_tot_mass_fit_params[3],
                                kds_disp_log_ang_mom_tot_mass_fit_params[0],
                                kds_disp_log_ang_mom_tot_mass_fit_params[0] - ang_mom_beta_zero,
                                np.sqrt((kds_disp_log_ang_mom_tot_mass_fit_params[0] - kds_disp_log_ang_mom_tot_mass_fit_params[1])**2 + ang_mom_tot_beta_zero_error**2),
                                np.sqrt((kds_disp_log_ang_mom_tot_mass_fit_params[2] - kds_disp_log_ang_mom_tot_mass_fit_params[0])**2 + ang_mom_tot_beta_zero_error**2),])
dyn_mass_evolution_list.append(['kds',
                                3.5,
                                kds_all_log_velocity_mass_fit_params[4],
                                kds_all_median_dynamical_mass_ratio,
                                kds_all_median_dynamical_mass_ratio_lower_error,
                                kds_all_median_dynamical_mass_ratio_upper_error,
                                kds_all_median_dynamical_mass_with_sigma_ratio,
                                kds_all_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kds_all_median_dynamical_mass_with_sigma_ratio_upper_error,
                                kds_rot_log_velocity_mass_fit_params[4],
                                kds_rot_median_dynamical_mass_ratio,
                                kds_rot_median_dynamical_mass_ratio_lower_error,
                                kds_rot_median_dynamical_mass_ratio_upper_error,
                                kds_rot_median_dynamical_mass_with_sigma_ratio,
                                kds_rot_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kds_rot_median_dynamical_mass_with_sigma_ratio_upper_error,
                                kds_disp_log_velocity_mass_fit_params[4],
                                kds_disp_median_dynamical_mass_ratio,
                                kds_disp_median_dynamical_mass_ratio_lower_error,
                                kds_disp_median_dynamical_mass_ratio_upper_error,
                                kds_disp_median_dynamical_mass_with_sigma_ratio,
                                kds_disp_median_dynamical_mass_with_sigma_ratio_lower_error,
                                kds_disp_median_dynamical_mass_with_sigma_ratio_upper_error])

# write out to the files (overwrite them if they exist)
# first the velocity evolution file
velocity_file_name = '/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/v_' + str(global_beta) + '_evolution.txt'
if os.path.isfile(velocity_file_name):

    os.system('rm %s' % velocity_file_name)

with open(velocity_file_name, 'a') as f:

    for entry in velocity_evolution_list:

        for column in entry:

            f.write('%s\t' % column)
        f.write('\n')
# v tot evolution file
v_tot_file_name = '/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/v_tot_' + str(global_beta) + '_evolution.txt'
if os.path.isfile(v_tot_file_name):

    os.system('rm %s' % v_tot_file_name)

with open(v_tot_file_name, 'a') as f:

    for entry in v_tot_evolution_list:

        for column in entry:

            f.write('%s\t' % column)
        f.write('\n')
# ang mom evolution file
ang_mom_file_name = '/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/ang_mom_' + str(global_beta) + '_evolution.txt'
if os.path.isfile(ang_mom_file_name):

    os.system('rm %s' % ang_mom_file_name)

with open(ang_mom_file_name, 'a') as f:

    for entry in ang_mom_evolution_list:

        for column in entry:

            f.write('%s\t' % column)
        f.write('\n')
# ang mom tot evolution file
ang_mom_tot_file_name = '/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/ang_mom_tot_' + str(global_beta) + '_evolution.txt'
if os.path.isfile(ang_mom_tot_file_name):

    os.system('rm %s' % ang_mom_tot_file_name)

with open(ang_mom_tot_file_name, 'a') as f:

    for entry in ang_mom_tot_evolution_list:

        for column in entry:

            f.write('%s\t' % column)
        f.write('\n')
# dyn mass evolution file
dyn_mass_file_name = '/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/dyn_mass_' + str(global_beta) + '_evolution.txt'
if os.path.isfile(dyn_mass_file_name):

    os.system('rm %s' % dyn_mass_file_name)

with open(dyn_mass_file_name, 'a') as f:

    for entry in dyn_mass_evolution_list:

        for column in entry:

            f.write('%s\t' % column)
        f.write('\n')


fig, ax = plt.subplots(1, 1, figsize=(8,8))

ax.set_xlabel(r'log(M$\boldsymbol{_{\star}}$/M$\boldsymbol{_{\odot}}$)',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax.set_ylabel(r'log(V$\boldsymbol{_{C}}$)',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax.tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax.tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

ax.plot(fit_mass,
        kross_rot_log_velocity_mass_fit_line,
        color='forestgreen',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 star-forming sample (Harrison+17)}')

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

eb1 = ax.errorbar(kds_rot_mass,
                  kds_rot_log_velocity,
                  yerr=[kds_rot_log_velocity_lower_error, kds_rot_log_velocity_upper_error],
                  ecolor='black',
                  marker='o',
                  markersize=10,
                  color='firebrick',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1 (this study)}')

eb1 = ax.errorbar(kross_all_mass,
                  kross_all_log_velocity,
                  yerr=[kross_all_log_velocity_lower_error, kross_all_log_velocity_upper_error],
                  ecolor='black',
                  marker='^',
                  markersize=10,
                  color='grey',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{dynamo}')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

ax.plot(fit_mass,
        kds_rot_log_velocity_mass_fit_line,
        color='firebrick',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1}')
ax.plot(fit_mass,
        kross_all_log_velocity_mass_fit_line,
        color='grey',
        lw=4,
        label=r'\textbf{dynamo fit}')
ax.fill_between(fit_mass,
                kds_rot_log_velocity_mass_fit_line_lower,
                kds_rot_log_velocity_mass_fit_line_upper,
                facecolor='firebrick',
                edgecolor='firebrick',
                alpha=0.25)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

eb2 = ax.errorbar(kds_disp_mass,
                  kds_disp_log_velocity,
                  yerr=[kds_disp_log_velocity_lower_error, kds_disp_log_velocity_upper_error],
                  ecolor='black',
                  marker='D',
                  markersize=10,
                  color='navajowhite',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1 (this study)}')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)
ax.plot(fit_mass,
        kds_disp_log_velocity_mass_fit_line,
        color='navajowhite',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1}')
ax.fill_between(fit_mass,
                kds_disp_log_velocity_mass_fit_line_lower,
                kds_disp_log_velocity_mass_fit_line_upper,
                facecolor='dimgrey',
                edgecolor='dimgrey',
                alpha=0.5)


[i.set_linewidth(4.0) for i in ax.spines.itervalues()]


# Tiley z ~ 1 line
lines = {'linestyle': '--'}
plt.rc('lines', **lines)


ax.plot(reyes_tf_mass,
        log_reyes_tf_velocity,
        label=r'\textbf{z$\boldsymbol{<}$0.1 disks (Reyes+11)}',
        color='black',
        lw=4.0)

cm = plt.cm.get_cmap('Greens')
a,b,d,im = ax.hist2d(kross_all_mass,
                     kross_all_log_velocity,
                     bins=20,
                     range=[[8.0,11.0],
                            [1.0, 2.5]],
                     norm=LogNorm(),
                     cmap=cm)
cbaxes = fig.add_axes([0.20, 0.91, 0.74, 0.04]) 
cb = plt.colorbar(im, cax = cbaxes, ticks=[1.0, 14, 40],orientation='horizontal')
cb.ax.set_xticklabels(['1', '15', '40'])
cb.ax.set_xlabel(r'\textbf{z $\boldsymbol{\simeq}$ 0.9 star-forming galaxies (Harrison+17)}',
                 fontsize='20',
                 weight='bold',
                 labelpad=-42)
for l in cb.ax.xaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(18)

ax.set_xlim(8.0,11.2)
ax.set_ylim(0.8,2.7)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

ax.legend(loc='lower right',
          prop={'size':16,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
# also plot systematic uncertainty on Mstar
ax.plot([8.5,8.9],[1.1,1.1],color='black',lw=4)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig.tight_layout()
ax.minorticks_on()

plt.show()

fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/tully_fisher_standard_log.png')

plt.close('all')

fig, ax = plt.subplots(1, 1, figsize=(8,8))

ax.set_xlabel(r'log(M$\boldsymbol{_{\star}}$/M$\boldsymbol{_{\odot}}$)',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax.set_ylabel(r'log$\left(\sqrt{V\boldsymbol{_{C}^{2} + 2.6\sigma_{int}^{2}}}\right)$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax.tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax.tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

ax.plot(fit_mass,
        kross_all_log_v_tot_mass_fit_line,
        color='forestgreen',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 star-forming galaxies (Harrison+17)}')

ax.plot(fit_mass,
        kross_all_log_v_tot_mass_fit_line,
        color='grey',
        lw=4,
        label=r'\textbf{test}')

ax.plot(fit_mass,
        kross_all_log_velocity_mass_fit_line,
        color='red',
        lw=4,
        label=r'\textbf{test}')

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

eb0 = ax.errorbar(kross_all_mass,
                  kross_all_log_v_tot,
                  yerr=[kross_all_log_v_tot_lower_error, kross_all_log_v_tot_upper_error],
                  ecolor='black',
                  marker='D',
                  markersize=10,
                  color='grey',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{Test}')

eb0 = ax.errorbar(kross_all_mass,
                  kross_all_log_velocity,
                  yerr=[kross_all_log_velocity_lower_error, kross_all_log_velocity_upper_error],
                  ecolor='black',
                  marker='D',
                  markersize=10,
                  color='red',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{Test}')

eb1 = ax.errorbar(kds_rot_mass,
                  kds_rot_log_v_tot,
                  yerr=[kds_rot_log_v_tot_lower_error, kds_rot_log_v_tot_upper_error],
                  ecolor='black',
                  marker='o',
                  markersize=10,
                  color='firebrick',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1 (this study)}')

eb2 = ax.errorbar(kds_disp_mass,
                  kds_disp_log_v_tot,
                  yerr=[kds_disp_log_v_tot_lower_error, kds_disp_log_v_tot_upper_error],
                  ecolor='black',
                  marker='D',
                  markersize=10,
                  color='navajowhite',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1 (this study)}')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

ax.plot(fit_mass,
        kds_all_log_v_tot_mass_fit_line,
        color='firebrick',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to full sample')
ax.fill_between(fit_mass,
                kds_all_log_v_tot_mass_fit_line_lower,
                kds_all_log_v_tot_mass_fit_line_upper,
                facecolor='firebrick',
                edgecolor='firebrick',
                alpha=0.25)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)


[i.set_linewidth(4.0) for i in ax.spines.itervalues()]


cm = plt.cm.get_cmap('Greens')
a,b,d,im = ax.hist2d(kross_all_mass,
                     kross_all_log_v_tot,
                     bins=20,
                     range=[[8.0,11.0],
                            [1.0, 2.5]],
                     norm=LogNorm(),
                     cmap=cm)
cbaxes = fig.add_axes([0.22, 0.36, 0.72, 0.04]) 
cb = plt.colorbar(im, cax = cbaxes, ticks=[1.0, 19, 40],orientation='horizontal')
cb.ax.set_xticklabels(['1', '20', '20'])
cb.ax.set_xlabel(r'\textbf{z $\boldsymbol{\simeq}$ 0.9 star-forming galaxies (Harrison+17)}',
                 fontsize='20',
                 weight='bold',
                 labelpad=-42)
for l in cb.ax.xaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(18)


# Reyes z ~ 0 line
lines = {'linestyle': '--'}
plt.rc('lines', **lines)

ax.plot(reyes_tf_mass,
        log_reyes_tf_velocity,
        label=r'\textbf{z$\boldsymbol{<}$0.1 disks (Reyes+11)}',
        color='black',
        lw=4.0)

ax.set_xlim(8.0,11.2)
ax.set_ylim(0.8,2.7)


ax.legend(loc='lower right',
          prop={'size':16,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
# also plot systematic uncertainty on Mstar
ax.plot([8.5,8.9],[1.1,1.1],color='black',lw=4)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig.tight_layout()
ax.minorticks_on()

plt.show()

fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/tully_fisher_with_sigma.png')

plt.close('all')

############ANG MOM PLOTS################

fig, ax = plt.subplots(1, 1, figsize=(8,8))

ax.set_xlabel(r'log(M$\boldsymbol{_{\star}}$/M$\boldsymbol{_{\odot}}$)',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax.set_ylabel(r'log(j$\boldsymbol{_{s}}$)',
              fontsize=28,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax.tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax.tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

cm = plt.cm.get_cmap('Greens')
a,b,d,im = ax.hist2d(kross_all_mass,
                     kross_all_log_ang_mom,
                     bins=20,
                     range=[[8.5,11.0],
                            [1.0, 3.0]],
                     norm=LogNorm(),
                     cmap=cm)
cbaxes = fig.add_axes([0.23, 0.2, 0.73, 0.04]) 
cb = plt.colorbar(im, cax = cbaxes, ticks=[1.0, 10, 40],orientation='horizontal')
cb.ax.set_xticklabels(['1', '10', '20'])
cb.ax.set_xlabel(r'\textbf{z $\boldsymbol{\simeq}$ 0.9 star-forming galaxies (Harrison+17)}',
                 fontsize='20',
                 weight='bold',
                 labelpad=-42)
for l in cb.ax.xaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(18)

ax.plot(reyes_tf_mass,
        log_fall_ang_mom,
        label=r'\textbf{z$\boldsymbol{=}$0 Spirals (Romanowsky \& Fall 12)}',
        color='black',
        lw=4)

#ax.plot(reyes_tf_mass,
#        ang_mom_ev,
#        label=r'\textbf{j$\boldsymbol{_{s}\propto (1 + z)^{-1}}$ evolution}',
#        color='red',
#        lw=4)

ax.plot(fit_mass,
        kross_all_log_ang_mom_mass_fit_line,
        color='forestgreen',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 all (Harrison+17)}')

lines = {'linestyle': '--'}
plt.rc('lines', **lines)

ax.plot(fit_mass,
        kross_rot_log_ang_mom_mass_fit_line,
        color='darkgreen',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1 (Harrison+17)}')

ax.plot(fit_mass,
        kross_disp_log_ang_mom_mass_fit_line,
        color='olive',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1 (Harrison+17)}')

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

#eb3 = ax.errorbar(amaze_all_mass,
#                  amaze_all_log_ang_mom,
#                  xerr=[amaze_all_mass_upper_error,amaze_all_mass_lower_error],
#                  yerr=[amaze_all_log_ang_mom_lower_error, amaze_all_log_ang_mom_upper_error],
#                  ecolor='black',
#                  marker='H',
#                  markersize=10,
#                  color='grey',
#                  capsize=2,
#                  elinewidth=2,
#                  label=r'\textbf{z$\boldsymbol{\sim}$3 AMAZE}')

eb1 = ax.errorbar(kds_rot_mass,
                  kds_rot_log_ang_mom,
                  yerr=[kds_rot_log_ang_mom_lower_error, kds_rot_log_ang_mom_upper_error],
                  ecolor='black',
                  marker='o',
                  markersize=10,
                  color='firebrick',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1 (this study)}')

eb1 = ax.errorbar(cresci_all_mass,
                  cresci_all_log_ang_mom,
                  yerr=[cresci_all_log_ang_mom_lower_error, cresci_all_log_ang_mom_upper_error],
                  ecolor='black',
                  marker='^',
                  markersize=10,
                  color='grey',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{massiv}')


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
ax.plot(fit_mass,
        kds_rot_log_ang_mom_mass_fit_line,
        color='firebrick',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1}')
ax.plot(fit_mass,
        cresci_all_log_ang_mom_mass_fit_line,
        color='grey',
        lw=4,
        label=r'\textbf{massiv fit line}')
ax.fill_between(fit_mass,
                kds_rot_log_ang_mom_mass_fit_line_lower,
                kds_rot_log_ang_mom_mass_fit_line_upper,
                facecolor='firebrick',
                edgecolor='firebrick',
                alpha=0.25)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

eb2 = ax.errorbar(kds_disp_mass,
                  kds_disp_log_ang_mom,
                  yerr=[kds_disp_log_ang_mom_lower_error, kds_disp_log_ang_mom_upper_error],
                  ecolor='black',
                  marker='D',
                  markersize=10,
                  color='navajowhite',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{z$\boldsymbol{\sim}$3.5 KDS V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1 (this study)}')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)
ax.plot(fit_mass,
        kds_disp_log_ang_mom_mass_fit_line,
        color='navajowhite',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1}')
ax.fill_between(fit_mass,
                kds_disp_log_ang_mom_mass_fit_line_lower,
                kds_disp_log_ang_mom_mass_fit_line_upper,
                facecolor='dimgrey',
                edgecolor='dimgrey',
                alpha=0.5)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)


[i.set_linewidth(4.0) for i in ax.spines.itervalues()]

#ax.scatter(log_vmax_kross_all,
#                  mass_kross_all,
#                  marker='.',
#                  color='black',
#                  label=r'KROSS: Harrison 2016')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)


lower_track = kross_all_log_ang_mom_mass_fit_line - 0.3
upper_track = kross_all_log_ang_mom_mass_fit_line + 0.3

#ax.fill_between(fit_mass,
#                lower_track,
#                upper_track,
#                facecolor='forestgreen',
#                edgecolor='forestgreen',
#                alpha=0.25)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)


#eb3[-1][0].set_linestyle('-.')
#eb3[-1][1].set_linestyle('-.')


# Tiley z ~ 1 line
lines = {'linestyle': '-'}
plt.rc('lines', **lines)

#ax.plot(tiley_v,
#        tiley_mass,
#        label=r'Tiley 2016 z $\sim$ 1',
#        color='orangered',
#        lw=2.0)



ax.set_xlim(8.0,11.2)
ax.set_ylim(0.5,3.8)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

ax.legend(loc='upper left',
          prop={'size':16,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
# also plot systematic uncertainty on Mstar
ax.plot([10.5,10.9],[1.4,1.4],color='black',lw=4)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig.tight_layout()
ax.minorticks_on()

plt.show()

fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/angular_momentum.png')

plt.close('all')

fig, ax = plt.subplots(1, 1, figsize=(8,8))

ax.set_xlabel(r'log(M$\boldsymbol{_{\star}}$/M$\boldsymbol{_{\odot}}$)',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax.set_ylabel(r'log(j$\boldsymbol{_{\textrm{tot}}}$)',
              fontsize=28,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax.tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax.tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

lines = {'linestyle': '-'}
plt.rc('lines', **lines)

ax.plot(reyes_tf_mass,
        log_fall_ang_mom,
        label=r'\textbf{z$\boldsymbol{=}$0 Spirals (Romanowsky \& Fall 12)}',
        color='black',
        lw=4)

#lines = {'linestyle': '--'}
#plt.rc('lines', **lines)
#ax.plot(reyes_tf_mass,
#        ang_mom_ev,
#        label=r'\textbf{predicted $\boldsymbol{(1 + z)^{-1}}$ evolution}',
#        color='black',
#        lw=4)

#ax.plot(reyes_tf_mass,
#        ang_mom_ev,
#        label=r'\textbf{j$\boldsymbol{_{s}\propto (1 + z)^{-1}}$ evolution}',
#        color='red',
#        lw=4)

ax.plot(fit_mass,
        kross_all_log_ang_mom_tot_mass_fit_line,
        color='forestgreen',
        lw=4,
        label=r'\textbf{z$\boldsymbol{\sim}$0.9 all (Harrison+17)}')

lines = {'linestyle': '--'}
plt.rc('lines', **lines)

#ax.plot(fit_mass,
#        kross_rot_log_ang_mom_mass_fit_line,
#        color='darkgreen',
#        lw=4,
#        label=r'\textbf{z$\boldsymbol{\sim}$0.9 V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{>}$ 1 (Harrison+17)}')
#ax.plot(fit_mass,
#        kross_disp_log_ang_mom_mass_fit_line,
#        color='olive',
#        lw=4,
#        label=r'\textbf{z$\boldsymbol{\sim}$0.9 V$\boldsymbol{_{C}}$/$\boldsymbol{\sigma_{int}}$ $\boldsymbol{<}$ 1 (Harrison+17)}')

cm = plt.cm.get_cmap('Greens')
a,b,d,im = ax.hist2d(kross_all_mass,
                     kross_all_log_ang_mom_tot,
                     bins=20,
                     range=[[8.5,11.0],
                            [1.0, 3.0]],
                     norm=LogNorm(),
                     cmap=cm)
cbaxes = fig.add_axes([0.23, 0.2, 0.73, 0.04]) 
cb = plt.colorbar(im, cax = cbaxes, ticks=[1.0, 10, 40],orientation='horizontal')
cb.ax.set_xticklabels(['1', '10', '20'])
cb.ax.set_xlabel(r'\textbf{z $\boldsymbol{\simeq}$ 0.9 star-forming galaxies (Harrison+17)}',
                 fontsize='20',
                 weight='bold',
                 labelpad=-42)
for l in cb.ax.xaxis.get_ticklabels():
    l.set_weight("bold")
    l.set_fontsize(18)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

#eb3 = ax.errorbar(amaze_all_mass,
#                  amaze_all_log_ang_mom,
#                  xerr=[amaze_all_mass_upper_error,amaze_all_mass_lower_error],
#                  yerr=[amaze_all_log_ang_mom_lower_error, amaze_all_log_ang_mom_upper_error],
#                  ecolor='black',
#                  marker='H',
#                  markersize=10,
#                  color='grey',
#                  capsize=2,
#                  elinewidth=2,
#                  label=r'\textbf{z$\boldsymbol{\sim}$3 AMAZE}')

eb1 = ax.errorbar(kds_all_mass,
                  kds_all_log_ang_mom_tot,
                  yerr=[kds_all_log_ang_mom_tot_lower_error, kds_all_log_ang_mom_tot_upper_error],
                  ecolor='black',
                  marker='o',
                  markersize=10,
                  color='firebrick',
                  capsize=2,
                  elinewidth=2,
                  label=r'\textbf{KDS full sample (this study)}')


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
ax.plot(fit_mass,
        kds_all_log_ang_mom_tot_mass_fit_line,
        color='firebrick',
        lw=4,
        label=r'\textbf{KDS fixed slope fit to full sample}')
ax.fill_between(fit_mass,
                kds_all_log_ang_mom_tot_mass_fit_line_lower,
                kds_all_log_ang_mom_tot_mass_fit_line_upper,
                facecolor='firebrick',
                edgecolor='firebrick',
                alpha=0.25)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

[i.set_linewidth(4.0) for i in ax.spines.itervalues()]

#ax.scatter(log_vmax_kross_all,
#                  mass_kross_all,
#                  marker='.',
#                  color='black',
#                  label=r'KROSS: Harrison 2016')

lines = {'linestyle': '-'}
plt.rc('lines', **lines)


lower_track = kross_all_log_ang_mom_mass_fit_line - 0.3
upper_track = kross_all_log_ang_mom_mass_fit_line + 0.3

#ax.fill_between(fit_mass,
#                lower_track,
#                upper_track,
#                facecolor='forestgreen',
#                edgecolor='forestgreen',
#                alpha=0.25)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)


#eb3[-1][0].set_linestyle('-.')
#eb3[-1][1].set_linestyle('-.')


# Tiley z ~ 1 line
lines = {'linestyle': '-'}
plt.rc('lines', **lines)

#ax.plot(tiley_v,
#        tiley_mass,
#        label=r'Tiley 2016 z $\sim$ 1',
#        color='orangered',
#        lw=2.0)



ax.set_xlim(8.0,11.2)
ax.set_ylim(0.5,3.8)

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

ax.legend(loc='upper left',
          prop={'size':16,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)


lines = {'linestyle': '-'}
plt.rc('lines', **lines)
# also plot systematic uncertainty on Mstar
ax.plot([10.5,10.9],[1.6,1.6],color='black',lw=4)
lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig.tight_layout()
ax.minorticks_on()

plt.show()

fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/angular_momentum_with_sigma.png')

plt.close('all')