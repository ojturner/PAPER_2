import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from scipy.stats import binned_statistic
from matplotlib import rc
import matplotlib.ticker as ticker
from lmfit import Model, Parameters
from matplotlib.ticker import ScalarFormatter
from pylab import MaxNLocator
from pylab import *
from matplotlib import pyplot


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', weight='bold')
rc('text', usetex=True)
rc('axes', linewidth=2)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

# read in the evolution tables
table_velocity_evolution = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/v_2.5_evolution.txt')
table_v_tot_evolution = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/v_tot_2.5_evolution.txt')
table_ang_mom_evolution = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/ang_mom_2.5_evolution.txt')
table_ang_mom_tot_evolution = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/ang_mom_tot_2.5_evolution.txt')
table_dynamical_mass_evolution = ascii.read('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/REDSHIFT_EVOLUTION_TABLES/dyn_mass_2.5_evolution.txt')

# make the velocity evolution figure

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'Velocity vs. Mass$\boldsymbol{\Delta \beta}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-0.6,0.6)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-0.6,0.6)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[0].set_ylim(-0.6,0.6)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

dynamo_redshift = table_velocity_evolution[5][1]
dynamo_velocity_number_all = table_velocity_evolution[5][2]
dynamo_velocity_chi_all = table_velocity_evolution[5][3]
dynamo_velocity_beta_all = table_velocity_evolution[5][4]
dynamo_velocity_delta_beta_all = table_velocity_evolution[5][5]
dynamo_velocity_delta_beta_lower_error_all = table_velocity_evolution[5][6]
dynamo_velocity_delta_beta_upper_error_all = table_velocity_evolution[5][7]
dynamo_velocity_number_rot = table_velocity_evolution[5][2]
dynamo_velocity_chi_rot = table_velocity_evolution[5][3]
dynamo_velocity_beta_rot = table_velocity_evolution[5][4]
dynamo_velocity_delta_beta_rot = table_velocity_evolution[5][5]
dynamo_velocity_delta_beta_lower_error_rot = table_velocity_evolution[5][6]
dynamo_velocity_delta_beta_upper_error_rot = table_velocity_evolution[5][7]
dynamo_velocity_number_disp = table_velocity_evolution[5][14]
dynamo_velocity_chi_disp = table_velocity_evolution[5][15]
dynamo_velocity_beta_disp = table_velocity_evolution[5][16]
dynamo_velocity_delta_beta_disp = table_velocity_evolution[5][17]
dynamo_velocity_delta_beta_lower_error_disp = table_velocity_evolution[5][18]
dynamo_velocity_delta_beta_upper_error_disp = table_velocity_evolution[5][19]

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_velocity_delta_beta_all,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_velocity_delta_beta_lower_error_all,dynamo_velocity_delta_beta_upper_error_all]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_velocity_delta_beta_rot,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_velocity_delta_beta_lower_error_rot,dynamo_velocity_delta_beta_upper_error_rot]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_low_redshift_redshift = table_velocity_evolution[3][1]
swinbank_low_redshift_velocity_number_all = table_velocity_evolution[3][2]
swinbank_low_redshift_velocity_chi_all = table_velocity_evolution[3][3]
swinbank_low_redshift_velocity_beta_all = table_velocity_evolution[3][4]
swinbank_low_redshift_velocity_delta_beta_all = table_velocity_evolution[3][5]
swinbank_low_redshift_velocity_delta_beta_lower_error_all = table_velocity_evolution[3][6]
swinbank_low_redshift_velocity_delta_beta_upper_error_all = table_velocity_evolution[3][7]
swinbank_low_redshift_velocity_number_rot = table_velocity_evolution[3][8]
swinbank_low_redshift_velocity_chi_rot = table_velocity_evolution[3][9]
swinbank_low_redshift_velocity_beta_rot = table_velocity_evolution[3][10]
swinbank_low_redshift_velocity_delta_beta_rot = table_velocity_evolution[3][11]
swinbank_low_redshift_velocity_delta_beta_lower_error_rot = table_velocity_evolution[3][12]
swinbank_low_redshift_velocity_delta_beta_upper_error_rot = table_velocity_evolution[3][13]
swinbank_low_redshift_velocity_number_disp = table_velocity_evolution[3][14]
swinbank_low_redshift_velocity_chi_disp = table_velocity_evolution[3][15]
swinbank_low_redshift_velocity_beta_disp = table_velocity_evolution[3][16]
swinbank_low_redshift_velocity_delta_beta_disp = table_velocity_evolution[3][17]
swinbank_low_redshift_velocity_delta_beta_lower_error_disp = table_velocity_evolution[3][18]
swinbank_low_redshift_velocity_delta_beta_upper_error_disp = table_velocity_evolution[3][19]

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_velocity_delta_beta_all,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_velocity_delta_beta_lower_error_all,swinbank_low_redshift_velocity_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_velocity_delta_beta_rot,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_velocity_delta_beta_lower_error_rot,swinbank_low_redshift_velocity_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_velocity_delta_beta_disp,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_velocity_delta_beta_lower_error_disp,swinbank_low_redshift_velocity_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

kross_redshift = table_velocity_evolution[6][1]
kross_velocity_number_all = table_velocity_evolution[6][2]
kross_velocity_chi_all = table_velocity_evolution[6][3]
kross_velocity_beta_all = table_velocity_evolution[6][4]
kross_velocity_delta_beta_all = table_velocity_evolution[6][5]
kross_velocity_delta_beta_lower_error_all = table_velocity_evolution[6][6]
kross_velocity_delta_beta_upper_error_all = table_velocity_evolution[6][7]
kross_velocity_number_rot = table_velocity_evolution[6][8]
kross_velocity_chi_rot = table_velocity_evolution[6][9]
kross_velocity_beta_rot = table_velocity_evolution[6][10]
kross_velocity_delta_beta_rot = table_velocity_evolution[6][11]
kross_velocity_delta_beta_lower_error_rot = table_velocity_evolution[6][12]
kross_velocity_delta_beta_upper_error_rot = table_velocity_evolution[6][13]
kross_velocity_number_disp = table_velocity_evolution[6][14]
kross_velocity_chi_disp = table_velocity_evolution[6][15]
kross_velocity_beta_disp = table_velocity_evolution[6][16]
kross_velocity_delta_beta_disp = table_velocity_evolution[6][17]
kross_velocity_delta_beta_lower_error_disp = table_velocity_evolution[6][18]
kross_velocity_delta_beta_upper_error_disp = table_velocity_evolution[6][19]

p3_all = ax[0].errorbar(kross_redshift,
            kross_velocity_delta_beta_all,
            ecolor='limegreen',
            yerr=np.array([[kross_velocity_delta_beta_lower_error_all,kross_velocity_delta_beta_upper_error_all]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_rot = ax[1].errorbar(kross_redshift,
            kross_velocity_delta_beta_rot,
            ecolor='limegreen',
            yerr=np.array([[kross_velocity_delta_beta_lower_error_rot,kross_velocity_delta_beta_upper_error_rot]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_disp = ax[2].errorbar(kross_redshift,
            kross_velocity_delta_beta_disp,
            ecolor='limegreen',
            yerr=np.array([[kross_velocity_delta_beta_lower_error_disp,kross_velocity_delta_beta_upper_error_disp]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

massiv_redshift = table_velocity_evolution[7][1]
massiv_velocity_number_all = table_velocity_evolution[7][2]
massiv_velocity_chi_all = table_velocity_evolution[7][3]
massiv_velocity_beta_all = table_velocity_evolution[7][4]
massiv_velocity_delta_beta_all = table_velocity_evolution[7][5]
massiv_velocity_delta_beta_lower_error_all = table_velocity_evolution[7][6]
massiv_velocity_delta_beta_upper_error_all = table_velocity_evolution[7][7]
massiv_velocity_number_rot = table_velocity_evolution[7][8]
massiv_velocity_chi_rot = table_velocity_evolution[7][9]
massiv_velocity_beta_rot = table_velocity_evolution[7][10]
massiv_velocity_delta_beta_rot = table_velocity_evolution[7][11]
massiv_velocity_delta_beta_lower_error_rot = table_velocity_evolution[7][12]
massiv_velocity_delta_beta_upper_error_rot = table_velocity_evolution[7][13]
massiv_velocity_number_disp = table_velocity_evolution[7][14]
massiv_velocity_chi_disp = table_velocity_evolution[7][15]
massiv_velocity_beta_disp = table_velocity_evolution[7][16]
massiv_velocity_delta_beta_disp = table_velocity_evolution[7][17]
massiv_velocity_delta_beta_lower_error_disp = table_velocity_evolution[7][18]
massiv_velocity_delta_beta_upper_error_disp = table_velocity_evolution[7][19]

p4_all = ax[0].errorbar(massiv_redshift,
            massiv_velocity_delta_beta_all,
            ecolor='olive',
            yerr=np.array([[massiv_velocity_delta_beta_lower_error_all,massiv_velocity_delta_beta_upper_error_all]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_rot = ax[1].errorbar(massiv_redshift,
            massiv_velocity_delta_beta_rot,
            ecolor='olive',
            yerr=np.array([[massiv_velocity_delta_beta_lower_error_rot,massiv_velocity_delta_beta_upper_error_rot]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_disp = ax[2].errorbar(massiv_redshift,
            massiv_velocity_delta_beta_disp,
            ecolor='olive',
            yerr=np.array([[massiv_velocity_delta_beta_lower_error_disp,massiv_velocity_delta_beta_upper_error_disp]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_high_redshift_redshift = table_velocity_evolution[4][1]
swinbank_high_redshift_velocity_number_all = table_velocity_evolution[4][2]
swinbank_high_redshift_velocity_chi_all = table_velocity_evolution[4][3]
swinbank_high_redshift_velocity_beta_all = table_velocity_evolution[4][4]
swinbank_high_redshift_velocity_delta_beta_all = table_velocity_evolution[4][5]
swinbank_high_redshift_velocity_delta_beta_lower_error_all = table_velocity_evolution[4][6]
swinbank_high_redshift_velocity_delta_beta_upper_error_all = table_velocity_evolution[4][7]
swinbank_high_redshift_velocity_number_rot = table_velocity_evolution[4][8]
swinbank_high_redshift_velocity_chi_rot = table_velocity_evolution[4][9]
swinbank_high_redshift_velocity_beta_rot = table_velocity_evolution[4][10]
swinbank_high_redshift_velocity_delta_beta_rot = table_velocity_evolution[4][11]
swinbank_high_redshift_velocity_delta_beta_lower_error_rot = table_velocity_evolution[4][12]
swinbank_high_redshift_velocity_delta_beta_upper_error_rot = table_velocity_evolution[4][13]
swinbank_high_redshift_velocity_number_disp = table_velocity_evolution[4][14]
swinbank_high_redshift_velocity_chi_disp = table_velocity_evolution[4][15]
swinbank_high_redshift_velocity_beta_disp = table_velocity_evolution[4][16]
swinbank_high_redshift_velocity_delta_beta_disp = table_velocity_evolution[4][17]
swinbank_high_redshift_velocity_delta_beta_lower_error_disp = table_velocity_evolution[4][18]
swinbank_high_redshift_velocity_delta_beta_upper_error_disp = table_velocity_evolution[4][19]

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_velocity_delta_beta_all,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_velocity_delta_beta_lower_error_all,swinbank_high_redshift_velocity_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_velocity_delta_beta_rot,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_velocity_delta_beta_lower_error_rot,swinbank_high_redshift_velocity_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_velocity_delta_beta_disp,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_velocity_delta_beta_lower_error_disp,swinbank_high_redshift_velocity_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)


sigma_low_redshift_redshift = table_velocity_evolution[1][1]
sigma_low_redshift_velocity_number_all = table_velocity_evolution[1][2]
sigma_low_redshift_velocity_chi_all = table_velocity_evolution[1][3]
sigma_low_redshift_velocity_beta_all = table_velocity_evolution[1][4]
sigma_low_redshift_velocity_delta_beta_all = table_velocity_evolution[1][5]
sigma_low_redshift_velocity_delta_beta_lower_error_all = table_velocity_evolution[1][6]
sigma_low_redshift_velocity_delta_beta_upper_error_all = table_velocity_evolution[1][7]
sigma_low_redshift_velocity_number_rot = table_velocity_evolution[1][8]
sigma_low_redshift_velocity_chi_rot = table_velocity_evolution[1][9]
sigma_low_redshift_velocity_beta_rot = table_velocity_evolution[1][10]
sigma_low_redshift_velocity_delta_beta_rot = table_velocity_evolution[1][11]
sigma_low_redshift_velocity_delta_beta_lower_error_rot = table_velocity_evolution[1][12]
sigma_low_redshift_velocity_delta_beta_upper_error_rot = table_velocity_evolution[1][13]
sigma_low_redshift_velocity_number_disp = table_velocity_evolution[1][14]
sigma_low_redshift_velocity_chi_disp = table_velocity_evolution[1][15]
sigma_low_redshift_velocity_beta_disp = table_velocity_evolution[1][16]
sigma_low_redshift_velocity_delta_beta_disp = table_velocity_evolution[1][17]
sigma_low_redshift_velocity_delta_beta_lower_error_disp = table_velocity_evolution[1][18]
sigma_low_redshift_velocity_delta_beta_upper_error_disp = table_velocity_evolution[1][19]

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_velocity_delta_beta_all,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_velocity_delta_beta_lower_error_all,sigma_low_redshift_velocity_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_velocity_delta_beta_rot,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_velocity_delta_beta_lower_error_rot,sigma_low_redshift_velocity_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_velocity_delta_beta_disp,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_velocity_delta_beta_lower_error_disp,sigma_low_redshift_velocity_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

cresci_redshift = table_velocity_evolution[0][1]
cresci_velocity_number_all = table_velocity_evolution[0][2]
cresci_velocity_chi_all = table_velocity_evolution[0][3]
cresci_velocity_beta_all = table_velocity_evolution[0][4]
cresci_velocity_delta_beta_all = table_velocity_evolution[0][5]
cresci_velocity_delta_beta_lower_error_all = table_velocity_evolution[0][6]
cresci_velocity_delta_beta_upper_error_all = table_velocity_evolution[0][7]
cresci_velocity_number_rot = table_velocity_evolution[0][2]
cresci_velocity_chi_rot = table_velocity_evolution[0][3]
cresci_velocity_beta_rot = table_velocity_evolution[0][4]
cresci_velocity_delta_beta_rot = table_velocity_evolution[0][5]
cresci_velocity_delta_beta_lower_error_rot = table_velocity_evolution[0][6]
cresci_velocity_delta_beta_upper_error_rot = table_velocity_evolution[0][7]
cresci_velocity_number_disp = table_velocity_evolution[0][14]
cresci_velocity_chi_disp = table_velocity_evolution[0][15]
cresci_velocity_beta_disp = table_velocity_evolution[0][16]
cresci_velocity_delta_beta_disp = table_velocity_evolution[0][17]
cresci_velocity_delta_beta_lower_error_disp = table_velocity_evolution[0][18]
cresci_velocity_delta_beta_upper_error_disp = table_velocity_evolution[0][19]

ax[0].errorbar(cresci_redshift,
            cresci_velocity_delta_beta_all,
            ecolor='darkorange',
            yerr=np.array([[cresci_velocity_delta_beta_lower_error_all,cresci_velocity_delta_beta_upper_error_all]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].errorbar(cresci_redshift,
            cresci_velocity_delta_beta_rot,
            ecolor='darkorange',
            yerr=np.array([[cresci_velocity_delta_beta_lower_error_rot,cresci_velocity_delta_beta_upper_error_rot]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

sigma_high_redshift_redshift = table_velocity_evolution[2][1]
sigma_high_redshift_velocity_number_all = table_velocity_evolution[2][2]
sigma_high_redshift_velocity_chi_all = table_velocity_evolution[2][3]
sigma_high_redshift_velocity_beta_all = table_velocity_evolution[2][4]
sigma_high_redshift_velocity_delta_beta_all = table_velocity_evolution[2][5]
sigma_high_redshift_velocity_delta_beta_lower_error_all = table_velocity_evolution[2][6]
sigma_high_redshift_velocity_delta_beta_upper_error_all = table_velocity_evolution[2][7]
sigma_high_redshift_velocity_number_rot = table_velocity_evolution[2][8]
sigma_high_redshift_velocity_chi_rot = table_velocity_evolution[2][9]
sigma_high_redshift_velocity_beta_rot = table_velocity_evolution[2][10]
sigma_high_redshift_velocity_delta_beta_rot = table_velocity_evolution[2][11]
sigma_high_redshift_velocity_delta_beta_lower_error_rot = table_velocity_evolution[2][12]
sigma_high_redshift_velocity_delta_beta_upper_error_rot = table_velocity_evolution[2][13]
sigma_high_redshift_velocity_number_disp = table_velocity_evolution[2][14]
sigma_high_redshift_velocity_chi_disp = table_velocity_evolution[2][15]
sigma_high_redshift_velocity_beta_disp = table_velocity_evolution[2][16]
sigma_high_redshift_velocity_delta_beta_disp = table_velocity_evolution[2][17]
sigma_high_redshift_velocity_delta_beta_lower_error_disp = table_velocity_evolution[2][18]
sigma_high_redshift_velocity_delta_beta_upper_error_disp = table_velocity_evolution[2][19]

ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_velocity_delta_beta_all,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_velocity_delta_beta_lower_error_all,sigma_high_redshift_velocity_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_velocity_delta_beta_rot,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_velocity_delta_beta_lower_error_rot,sigma_high_redshift_velocity_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_velocity_delta_beta_disp,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_velocity_delta_beta_lower_error_disp,sigma_high_redshift_velocity_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

amaze_redshift = table_velocity_evolution[8][1]
amaze_velocity_number_all = table_velocity_evolution[8][2]
amaze_velocity_chi_all = table_velocity_evolution[8][3]
amaze_velocity_beta_all = table_velocity_evolution[8][4]
amaze_velocity_delta_beta_all = table_velocity_evolution[8][5]
amaze_velocity_delta_beta_lower_error_all = table_velocity_evolution[8][6]
amaze_velocity_delta_beta_upper_error_all = table_velocity_evolution[8][7]
amaze_velocity_number_rot = table_velocity_evolution[8][2]
amaze_velocity_chi_rot = table_velocity_evolution[8][3]
amaze_velocity_beta_rot = table_velocity_evolution[8][4]
amaze_velocity_delta_beta_rot = table_velocity_evolution[8][5]
amaze_velocity_delta_beta_lower_error_rot = table_velocity_evolution[8][6]
amaze_velocity_delta_beta_upper_error_rot = table_velocity_evolution[8][7]
amaze_velocity_number_disp = table_velocity_evolution[8][14]
amaze_velocity_chi_disp = table_velocity_evolution[8][15]
amaze_velocity_beta_disp = table_velocity_evolution[8][16]
amaze_velocity_delta_beta_disp = table_velocity_evolution[8][17]
amaze_velocity_delta_beta_lower_error_disp = table_velocity_evolution[8][18]
amaze_velocity_delta_beta_upper_error_disp = table_velocity_evolution[8][19]

ax[0].errorbar(amaze_redshift,
            amaze_velocity_delta_beta_all,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_velocity_delta_beta_lower_error_all,amaze_velocity_delta_beta_upper_error_all]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].errorbar(amaze_redshift,
            amaze_velocity_delta_beta_rot,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_velocity_delta_beta_lower_error_rot,amaze_velocity_delta_beta_upper_error_rot]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

kds_redshift = table_velocity_evolution[9][1]
kds_velocity_number_all = table_velocity_evolution[9][2]
kds_velocity_chi_all = table_velocity_evolution[9][3]
kds_velocity_beta_all = table_velocity_evolution[9][4]
kds_velocity_delta_beta_all = table_velocity_evolution[9][5]
kds_velocity_delta_beta_lower_error_all = table_velocity_evolution[9][6]
kds_velocity_delta_beta_upper_error_all = table_velocity_evolution[9][7]
kds_velocity_number_rot = table_velocity_evolution[9][8]
kds_velocity_chi_rot = table_velocity_evolution[9][9]
kds_velocity_beta_rot = table_velocity_evolution[9][10]
kds_velocity_delta_beta_rot = table_velocity_evolution[9][11]
kds_velocity_delta_beta_lower_error_rot = table_velocity_evolution[9][12]
kds_velocity_delta_beta_upper_error_rot = table_velocity_evolution[9][13]
kds_velocity_number_disp = table_velocity_evolution[9][14]
kds_velocity_chi_disp = table_velocity_evolution[9][15]
kds_velocity_beta_disp = table_velocity_evolution[9][16]
kds_velocity_delta_beta_disp = table_velocity_evolution[9][17]
kds_velocity_delta_beta_lower_error_disp = table_velocity_evolution[9][18]
kds_velocity_delta_beta_upper_error_disp = table_velocity_evolution[9][19]

ax[0].errorbar(kds_redshift,
            kds_velocity_delta_beta_all,
            ecolor='firebrick',
            yerr=np.array([[kds_velocity_delta_beta_lower_error_all,kds_velocity_delta_beta_upper_error_all]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].errorbar(kds_redshift,
            kds_velocity_delta_beta_rot,
            ecolor='firebrick',
            yerr=np.array([[kds_velocity_delta_beta_lower_error_rot,kds_velocity_delta_beta_upper_error_rot]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].errorbar(kds_redshift,
            kds_velocity_delta_beta_disp,
            ecolor='firebrick',
            yerr=np.array([[kds_velocity_delta_beta_lower_error_disp,kds_velocity_delta_beta_upper_error_disp]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('upper left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='upper right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/velocity_evolution.png')
plt.close('all')

# make the v_tot evolution figure

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'V$\boldsymbol{_{tot}}$ vs. Mass$\boldsymbol{\Delta \beta}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-0.6,0.6)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-0.6,0.6)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[2].set_ylim(-0.6,0.6)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

dynamo_redshift = table_v_tot_evolution[5][1]
dynamo_v_tot_number_all = table_v_tot_evolution[5][2]
dynamo_v_tot_chi_all = table_v_tot_evolution[5][3]
dynamo_v_tot_beta_all = table_v_tot_evolution[5][4]
dynamo_v_tot_delta_beta_all = table_v_tot_evolution[5][5]
dynamo_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[5][6]
dynamo_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[5][7]
dynamo_v_tot_number_rot = table_v_tot_evolution[5][2]
dynamo_v_tot_chi_rot = table_v_tot_evolution[5][3]
dynamo_v_tot_beta_rot = table_v_tot_evolution[5][4]
dynamo_v_tot_delta_beta_rot = table_v_tot_evolution[5][5]
dynamo_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[5][6]
dynamo_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[5][7]
dynamo_v_tot_number_disp = table_v_tot_evolution[5][14]
dynamo_v_tot_chi_disp = table_v_tot_evolution[5][15]
dynamo_v_tot_beta_disp = table_v_tot_evolution[5][16]
dynamo_v_tot_delta_beta_disp = table_v_tot_evolution[5][17]
dynamo_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[5][18]
dynamo_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[5][19]

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_v_tot_delta_beta_all,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_v_tot_delta_beta_lower_error_all,dynamo_v_tot_delta_beta_upper_error_all]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_v_tot_delta_beta_rot,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_v_tot_delta_beta_lower_error_rot,dynamo_v_tot_delta_beta_upper_error_rot]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_low_redshift_redshift = table_v_tot_evolution[3][1]
swinbank_low_redshift_v_tot_number_all = table_v_tot_evolution[3][2]
swinbank_low_redshift_v_tot_chi_all = table_v_tot_evolution[3][3]
swinbank_low_redshift_v_tot_beta_all = table_v_tot_evolution[3][4]
swinbank_low_redshift_v_tot_delta_beta_all = table_v_tot_evolution[3][5]
swinbank_low_redshift_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[3][6]
swinbank_low_redshift_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[3][7]
swinbank_low_redshift_v_tot_number_rot = table_v_tot_evolution[3][8]
swinbank_low_redshift_v_tot_chi_rot = table_v_tot_evolution[3][9]
swinbank_low_redshift_v_tot_beta_rot = table_v_tot_evolution[3][10]
swinbank_low_redshift_v_tot_delta_beta_rot = table_v_tot_evolution[3][11]
swinbank_low_redshift_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[3][12]
swinbank_low_redshift_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[3][13]
swinbank_low_redshift_v_tot_number_disp = table_v_tot_evolution[3][14]
swinbank_low_redshift_v_tot_chi_disp = table_v_tot_evolution[3][15]
swinbank_low_redshift_v_tot_beta_disp = table_v_tot_evolution[3][16]
swinbank_low_redshift_v_tot_delta_beta_disp = table_v_tot_evolution[3][17]
swinbank_low_redshift_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[3][18]
swinbank_low_redshift_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[3][19]

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_v_tot_delta_beta_all,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_v_tot_delta_beta_lower_error_all,swinbank_low_redshift_v_tot_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_v_tot_delta_beta_rot,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_v_tot_delta_beta_lower_error_rot,swinbank_low_redshift_v_tot_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_v_tot_delta_beta_disp,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_v_tot_delta_beta_lower_error_disp,swinbank_low_redshift_v_tot_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

kross_redshift = table_v_tot_evolution[6][1]
kross_v_tot_number_all = table_v_tot_evolution[6][2]
kross_v_tot_chi_all = table_v_tot_evolution[6][3]
kross_v_tot_beta_all = table_v_tot_evolution[6][4]
kross_v_tot_delta_beta_all = table_v_tot_evolution[6][5]
kross_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[6][6]
kross_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[6][7]
kross_v_tot_number_rot = table_v_tot_evolution[6][8]
kross_v_tot_chi_rot = table_v_tot_evolution[6][9]
kross_v_tot_beta_rot = table_v_tot_evolution[6][10]
kross_v_tot_delta_beta_rot = table_v_tot_evolution[6][11]
kross_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[6][12]
kross_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[6][13]
kross_v_tot_number_disp = table_v_tot_evolution[6][14]
kross_v_tot_chi_disp = table_v_tot_evolution[6][15]
kross_v_tot_beta_disp = table_v_tot_evolution[6][16]
kross_v_tot_delta_beta_disp = table_v_tot_evolution[6][17]
kross_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[6][18]
kross_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[6][19]

p3_all = ax[0].errorbar(kross_redshift,
            kross_v_tot_delta_beta_all,
            ecolor='limegreen',
            yerr=np.array([[kross_v_tot_delta_beta_lower_error_all,kross_v_tot_delta_beta_upper_error_all]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_rot = ax[1].errorbar(kross_redshift,
            kross_v_tot_delta_beta_rot,
            ecolor='limegreen',
            yerr=np.array([[kross_v_tot_delta_beta_lower_error_rot,kross_v_tot_delta_beta_upper_error_rot]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_disp = ax[2].errorbar(kross_redshift,
            kross_v_tot_delta_beta_disp,
            ecolor='limegreen',
            yerr=np.array([[kross_v_tot_delta_beta_lower_error_disp,kross_v_tot_delta_beta_upper_error_disp]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

massiv_redshift = table_v_tot_evolution[7][1]
massiv_v_tot_number_all = table_v_tot_evolution[7][2]
massiv_v_tot_chi_all = table_v_tot_evolution[7][3]
massiv_v_tot_beta_all = table_v_tot_evolution[7][4]
massiv_v_tot_delta_beta_all = table_v_tot_evolution[7][5]
massiv_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[7][6]
massiv_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[7][7]
massiv_v_tot_number_rot = table_v_tot_evolution[7][8]
massiv_v_tot_chi_rot = table_v_tot_evolution[7][9]
massiv_v_tot_beta_rot = table_v_tot_evolution[7][10]
massiv_v_tot_delta_beta_rot = table_v_tot_evolution[7][11]
massiv_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[7][12]
massiv_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[7][13]
massiv_v_tot_number_disp = table_v_tot_evolution[7][14]
massiv_v_tot_chi_disp = table_v_tot_evolution[7][15]
massiv_v_tot_beta_disp = table_v_tot_evolution[7][16]
massiv_v_tot_delta_beta_disp = table_v_tot_evolution[7][17]
massiv_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[7][18]
massiv_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[7][19]

p4_all = ax[0].errorbar(massiv_redshift,
            massiv_v_tot_delta_beta_all,
            ecolor='olive',
            yerr=np.array([[massiv_v_tot_delta_beta_lower_error_all,massiv_v_tot_delta_beta_upper_error_all]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_rot = ax[1].errorbar(massiv_redshift,
            massiv_v_tot_delta_beta_rot,
            ecolor='olive',
            yerr=np.array([[massiv_v_tot_delta_beta_lower_error_rot,massiv_v_tot_delta_beta_upper_error_rot]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_disp = ax[2].errorbar(massiv_redshift,
            massiv_v_tot_delta_beta_disp,
            ecolor='olive',
            yerr=np.array([[massiv_v_tot_delta_beta_lower_error_disp,massiv_v_tot_delta_beta_upper_error_disp]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_high_redshift_redshift = table_v_tot_evolution[4][1]
swinbank_high_redshift_v_tot_number_all = table_v_tot_evolution[4][2]
swinbank_high_redshift_v_tot_chi_all = table_v_tot_evolution[4][3]
swinbank_high_redshift_v_tot_beta_all = table_v_tot_evolution[4][4]
swinbank_high_redshift_v_tot_delta_beta_all = table_v_tot_evolution[4][5]
swinbank_high_redshift_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[4][6]
swinbank_high_redshift_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[4][7]
swinbank_high_redshift_v_tot_number_rot = table_v_tot_evolution[4][8]
swinbank_high_redshift_v_tot_chi_rot = table_v_tot_evolution[4][9]
swinbank_high_redshift_v_tot_beta_rot = table_v_tot_evolution[4][10]
swinbank_high_redshift_v_tot_delta_beta_rot = table_v_tot_evolution[4][11]
swinbank_high_redshift_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[4][12]
swinbank_high_redshift_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[4][13]
swinbank_high_redshift_v_tot_number_disp = table_v_tot_evolution[4][14]
swinbank_high_redshift_v_tot_chi_disp = table_v_tot_evolution[4][15]
swinbank_high_redshift_v_tot_beta_disp = table_v_tot_evolution[4][16]
swinbank_high_redshift_v_tot_delta_beta_disp = table_v_tot_evolution[4][17]
swinbank_high_redshift_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[4][18]
swinbank_high_redshift_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[4][19]

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_v_tot_delta_beta_all,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_v_tot_delta_beta_lower_error_all,swinbank_high_redshift_v_tot_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_v_tot_delta_beta_rot,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_v_tot_delta_beta_lower_error_rot,swinbank_high_redshift_v_tot_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_v_tot_delta_beta_disp,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_v_tot_delta_beta_lower_error_disp,swinbank_high_redshift_v_tot_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)


sigma_low_redshift_redshift = table_v_tot_evolution[1][1]
sigma_low_redshift_v_tot_number_all = table_v_tot_evolution[1][2]
sigma_low_redshift_v_tot_chi_all = table_v_tot_evolution[1][3]
sigma_low_redshift_v_tot_beta_all = table_v_tot_evolution[1][4]
sigma_low_redshift_v_tot_delta_beta_all = table_v_tot_evolution[1][5]
sigma_low_redshift_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[1][6]
sigma_low_redshift_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[1][7]
sigma_low_redshift_v_tot_number_rot = table_v_tot_evolution[1][8]
sigma_low_redshift_v_tot_chi_rot = table_v_tot_evolution[1][9]
sigma_low_redshift_v_tot_beta_rot = table_v_tot_evolution[1][10]
sigma_low_redshift_v_tot_delta_beta_rot = table_v_tot_evolution[1][11]
sigma_low_redshift_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[1][12]
sigma_low_redshift_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[1][13]
sigma_low_redshift_v_tot_number_disp = table_v_tot_evolution[1][14]
sigma_low_redshift_v_tot_chi_disp = table_v_tot_evolution[1][15]
sigma_low_redshift_v_tot_beta_disp = table_v_tot_evolution[1][16]
sigma_low_redshift_v_tot_delta_beta_disp = table_v_tot_evolution[1][17]
sigma_low_redshift_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[1][18]
sigma_low_redshift_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[1][19]

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_v_tot_delta_beta_all,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_v_tot_delta_beta_lower_error_all,sigma_low_redshift_v_tot_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_v_tot_delta_beta_rot,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_v_tot_delta_beta_lower_error_rot,sigma_low_redshift_v_tot_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_v_tot_delta_beta_disp,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_v_tot_delta_beta_lower_error_disp,sigma_low_redshift_v_tot_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

cresci_redshift = table_v_tot_evolution[0][1]
cresci_v_tot_number_all = table_v_tot_evolution[0][2]
cresci_v_tot_chi_all = table_v_tot_evolution[0][3]
cresci_v_tot_beta_all = table_v_tot_evolution[0][4]
cresci_v_tot_delta_beta_all = table_v_tot_evolution[0][5]
cresci_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[0][6]
cresci_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[0][7]
cresci_v_tot_number_rot = table_v_tot_evolution[0][2]
cresci_v_tot_chi_rot = table_v_tot_evolution[0][3]
cresci_v_tot_beta_rot = table_v_tot_evolution[0][4]
cresci_v_tot_delta_beta_rot = table_v_tot_evolution[0][5]
cresci_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[0][6]
cresci_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[0][7]
cresci_v_tot_number_disp = table_v_tot_evolution[0][14]
cresci_v_tot_chi_disp = table_v_tot_evolution[0][15]
cresci_v_tot_beta_disp = table_v_tot_evolution[0][16]
cresci_v_tot_delta_beta_disp = table_v_tot_evolution[0][17]
cresci_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[0][18]
cresci_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[0][19]

ax[0].errorbar(cresci_redshift,
            cresci_v_tot_delta_beta_all,
            ecolor='darkorange',
            yerr=np.array([[cresci_v_tot_delta_beta_lower_error_all,cresci_v_tot_delta_beta_upper_error_all]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].errorbar(cresci_redshift,
            cresci_v_tot_delta_beta_rot,
            ecolor='darkorange',
            yerr=np.array([[cresci_v_tot_delta_beta_lower_error_rot,cresci_v_tot_delta_beta_upper_error_rot]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

sigma_high_redshift_redshift = table_v_tot_evolution[2][1]
sigma_high_redshift_v_tot_number_all = table_v_tot_evolution[2][2]
sigma_high_redshift_v_tot_chi_all = table_v_tot_evolution[2][3]
sigma_high_redshift_v_tot_beta_all = table_v_tot_evolution[2][4]
sigma_high_redshift_v_tot_delta_beta_all = table_v_tot_evolution[2][5]
sigma_high_redshift_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[2][6]
sigma_high_redshift_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[2][7]
sigma_high_redshift_v_tot_number_rot = table_v_tot_evolution[2][8]
sigma_high_redshift_v_tot_chi_rot = table_v_tot_evolution[2][9]
sigma_high_redshift_v_tot_beta_rot = table_v_tot_evolution[2][10]
sigma_high_redshift_v_tot_delta_beta_rot = table_v_tot_evolution[2][11]
sigma_high_redshift_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[2][12]
sigma_high_redshift_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[2][13]
sigma_high_redshift_v_tot_number_disp = table_v_tot_evolution[2][14]
sigma_high_redshift_v_tot_chi_disp = table_v_tot_evolution[2][15]
sigma_high_redshift_v_tot_beta_disp = table_v_tot_evolution[2][16]
sigma_high_redshift_v_tot_delta_beta_disp = table_v_tot_evolution[2][17]
sigma_high_redshift_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[2][18]
sigma_high_redshift_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[2][19]

ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_v_tot_delta_beta_all,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_v_tot_delta_beta_lower_error_all,sigma_high_redshift_v_tot_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_v_tot_delta_beta_rot,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_v_tot_delta_beta_lower_error_rot,sigma_high_redshift_v_tot_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_v_tot_delta_beta_disp,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_v_tot_delta_beta_lower_error_disp,sigma_high_redshift_v_tot_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

amaze_redshift = table_v_tot_evolution[8][1]
amaze_v_tot_number_all = table_v_tot_evolution[8][2]
amaze_v_tot_chi_all = table_v_tot_evolution[8][3]
amaze_v_tot_beta_all = table_v_tot_evolution[8][4]
amaze_v_tot_delta_beta_all = table_v_tot_evolution[8][5]
amaze_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[8][6]
amaze_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[8][7]
amaze_v_tot_number_rot = table_v_tot_evolution[8][2]
amaze_v_tot_chi_rot = table_v_tot_evolution[8][3]
amaze_v_tot_beta_rot = table_v_tot_evolution[8][4]
amaze_v_tot_delta_beta_rot = table_v_tot_evolution[8][5]
amaze_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[8][6]
amaze_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[8][7]
amaze_v_tot_number_disp = table_v_tot_evolution[8][14]
amaze_v_tot_chi_disp = table_v_tot_evolution[8][15]
amaze_v_tot_beta_disp = table_v_tot_evolution[8][16]
amaze_v_tot_delta_beta_disp = table_v_tot_evolution[8][17]
amaze_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[8][18]
amaze_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[8][19]

ax[0].errorbar(amaze_redshift,
            amaze_v_tot_delta_beta_all,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_v_tot_delta_beta_lower_error_all,amaze_v_tot_delta_beta_upper_error_all]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].errorbar(amaze_redshift,
            amaze_v_tot_delta_beta_rot,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_v_tot_delta_beta_lower_error_rot,amaze_v_tot_delta_beta_upper_error_rot]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

kds_redshift = table_v_tot_evolution[9][1]
kds_v_tot_number_all = table_v_tot_evolution[9][2]
kds_v_tot_chi_all = table_v_tot_evolution[9][3]
kds_v_tot_beta_all = table_v_tot_evolution[9][4]
kds_v_tot_delta_beta_all = table_v_tot_evolution[9][5]
kds_v_tot_delta_beta_lower_error_all = table_v_tot_evolution[9][6]
kds_v_tot_delta_beta_upper_error_all = table_v_tot_evolution[9][7]
kds_v_tot_number_rot = table_v_tot_evolution[9][8]
kds_v_tot_chi_rot = table_v_tot_evolution[9][9]
kds_v_tot_beta_rot = table_v_tot_evolution[9][10]
kds_v_tot_delta_beta_rot = table_v_tot_evolution[9][11]
kds_v_tot_delta_beta_lower_error_rot = table_v_tot_evolution[9][12]
kds_v_tot_delta_beta_upper_error_rot = table_v_tot_evolution[9][13]
kds_v_tot_number_disp = table_v_tot_evolution[9][14]
kds_v_tot_chi_disp = table_v_tot_evolution[9][15]
kds_v_tot_beta_disp = table_v_tot_evolution[9][16]
kds_v_tot_delta_beta_disp = table_v_tot_evolution[9][17]
kds_v_tot_delta_beta_lower_error_disp = table_v_tot_evolution[9][18]
kds_v_tot_delta_beta_upper_error_disp = table_v_tot_evolution[9][19]

ax[0].errorbar(kds_redshift,
            kds_v_tot_delta_beta_all,
            ecolor='firebrick',
            yerr=np.array([[kds_v_tot_delta_beta_lower_error_all,kds_v_tot_delta_beta_upper_error_all]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].errorbar(kds_redshift,
            kds_v_tot_delta_beta_rot,
            ecolor='firebrick',
            yerr=np.array([[kds_v_tot_delta_beta_lower_error_rot,kds_v_tot_delta_beta_upper_error_rot]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].errorbar(kds_redshift,
            kds_v_tot_delta_beta_disp,
            ecolor='firebrick',
            yerr=np.array([[kds_v_tot_delta_beta_lower_error_disp,kds_v_tot_delta_beta_upper_error_disp]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('upper left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='upper right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/v_tot_evolution.png')
plt.close('all')

# make the ang_mom evolution figure

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'j$\boldsymbol{_{gas}}$ vs. Mass$\boldsymbol{\Delta \beta}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-1.0,0.6)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-1.0,0.6)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[2].set_ylim(-1.0,0.6)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

dynamo_redshift = table_ang_mom_evolution[5][1]
dynamo_ang_mom_number_all = table_ang_mom_evolution[5][2]
dynamo_ang_mom_chi_all = table_ang_mom_evolution[5][3]
dynamo_ang_mom_beta_all = table_ang_mom_evolution[5][4]
dynamo_ang_mom_delta_beta_all = table_ang_mom_evolution[5][5]
dynamo_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[5][6]
dynamo_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[5][7]
dynamo_ang_mom_number_rot = table_ang_mom_evolution[5][2]
dynamo_ang_mom_chi_rot = table_ang_mom_evolution[5][3]
dynamo_ang_mom_beta_rot = table_ang_mom_evolution[5][4]
dynamo_ang_mom_delta_beta_rot = table_ang_mom_evolution[5][5]
dynamo_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[5][6]
dynamo_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[5][7]
dynamo_ang_mom_number_disp = table_ang_mom_evolution[5][14]
dynamo_ang_mom_chi_disp = table_ang_mom_evolution[5][15]
dynamo_ang_mom_beta_disp = table_ang_mom_evolution[5][16]
dynamo_ang_mom_delta_beta_disp = table_ang_mom_evolution[5][17]
dynamo_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[5][18]
dynamo_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[5][19]

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_ang_mom_delta_beta_all,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_ang_mom_delta_beta_lower_error_all,dynamo_ang_mom_delta_beta_upper_error_all]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_ang_mom_delta_beta_rot,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_ang_mom_delta_beta_lower_error_rot,dynamo_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_low_redshift_redshift = table_ang_mom_evolution[3][1]
swinbank_low_redshift_ang_mom_number_all = table_ang_mom_evolution[3][2]
swinbank_low_redshift_ang_mom_chi_all = table_ang_mom_evolution[3][3]
swinbank_low_redshift_ang_mom_beta_all = table_ang_mom_evolution[3][4]
swinbank_low_redshift_ang_mom_delta_beta_all = table_ang_mom_evolution[3][5]
swinbank_low_redshift_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[3][6]
swinbank_low_redshift_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[3][7]
swinbank_low_redshift_ang_mom_number_rot = table_ang_mom_evolution[3][8]
swinbank_low_redshift_ang_mom_chi_rot = table_ang_mom_evolution[3][9]
swinbank_low_redshift_ang_mom_beta_rot = table_ang_mom_evolution[3][10]
swinbank_low_redshift_ang_mom_delta_beta_rot = table_ang_mom_evolution[3][11]
swinbank_low_redshift_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[3][12]
swinbank_low_redshift_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[3][13]
swinbank_low_redshift_ang_mom_number_disp = table_ang_mom_evolution[3][14]
swinbank_low_redshift_ang_mom_chi_disp = table_ang_mom_evolution[3][15]
swinbank_low_redshift_ang_mom_beta_disp = table_ang_mom_evolution[3][16]
swinbank_low_redshift_ang_mom_delta_beta_disp = table_ang_mom_evolution[3][17]
swinbank_low_redshift_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[3][18]
swinbank_low_redshift_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[3][19]

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_delta_beta_all,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_delta_beta_lower_error_all,swinbank_low_redshift_ang_mom_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_delta_beta_rot,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_delta_beta_lower_error_rot,swinbank_low_redshift_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_delta_beta_disp,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_delta_beta_lower_error_disp,swinbank_low_redshift_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

kross_redshift = table_ang_mom_evolution[6][1]
kross_ang_mom_number_all = table_ang_mom_evolution[6][2]
kross_ang_mom_chi_all = table_ang_mom_evolution[6][3]
kross_ang_mom_beta_all = table_ang_mom_evolution[6][4]
kross_ang_mom_delta_beta_all = table_ang_mom_evolution[6][5]
kross_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[6][6]
kross_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[6][7]
kross_ang_mom_number_rot = table_ang_mom_evolution[6][8]
kross_ang_mom_chi_rot = table_ang_mom_evolution[6][9]
kross_ang_mom_beta_rot = table_ang_mom_evolution[6][10]
kross_ang_mom_delta_beta_rot = table_ang_mom_evolution[6][11]
kross_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[6][12]
kross_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[6][13]
kross_ang_mom_number_disp = table_ang_mom_evolution[6][14]
kross_ang_mom_chi_disp = table_ang_mom_evolution[6][15]
kross_ang_mom_beta_disp = table_ang_mom_evolution[6][16]
kross_ang_mom_delta_beta_disp = table_ang_mom_evolution[6][17]
kross_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[6][18]
kross_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[6][19]

p3_all = ax[0].errorbar(kross_redshift,
            kross_ang_mom_delta_beta_all,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_delta_beta_lower_error_all,kross_ang_mom_delta_beta_upper_error_all]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_rot = ax[1].errorbar(kross_redshift,
            kross_ang_mom_delta_beta_rot,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_delta_beta_lower_error_rot,kross_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_disp = ax[2].errorbar(kross_redshift,
            kross_ang_mom_delta_beta_disp,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_delta_beta_lower_error_disp,kross_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

massiv_redshift = table_ang_mom_evolution[7][1]
massiv_ang_mom_number_all = table_ang_mom_evolution[7][2]
massiv_ang_mom_chi_all = table_ang_mom_evolution[7][3]
massiv_ang_mom_beta_all = table_ang_mom_evolution[7][4]
massiv_ang_mom_delta_beta_all = table_ang_mom_evolution[7][5]
massiv_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[7][6]
massiv_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[7][7]
massiv_ang_mom_number_rot = table_ang_mom_evolution[7][8]
massiv_ang_mom_chi_rot = table_ang_mom_evolution[7][9]
massiv_ang_mom_beta_rot = table_ang_mom_evolution[7][10]
massiv_ang_mom_delta_beta_rot = table_ang_mom_evolution[7][11]
massiv_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[7][12]
massiv_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[7][13]
massiv_ang_mom_number_disp = table_ang_mom_evolution[7][14]
massiv_ang_mom_chi_disp = table_ang_mom_evolution[7][15]
massiv_ang_mom_beta_disp = table_ang_mom_evolution[7][16]
massiv_ang_mom_delta_beta_disp = table_ang_mom_evolution[7][17]
massiv_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[7][18]
massiv_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[7][19]

p4_all = ax[0].errorbar(massiv_redshift,
            massiv_ang_mom_delta_beta_all,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_delta_beta_lower_error_all,massiv_ang_mom_delta_beta_upper_error_all]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_rot = ax[1].errorbar(massiv_redshift,
            massiv_ang_mom_delta_beta_rot,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_delta_beta_lower_error_rot,massiv_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_disp = ax[2].errorbar(massiv_redshift,
            massiv_ang_mom_delta_beta_disp,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_delta_beta_lower_error_disp,massiv_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_high_redshift_redshift = table_ang_mom_evolution[4][1]
swinbank_high_redshift_ang_mom_number_all = table_ang_mom_evolution[4][2]
swinbank_high_redshift_ang_mom_chi_all = table_ang_mom_evolution[4][3]
swinbank_high_redshift_ang_mom_beta_all = table_ang_mom_evolution[4][4]
swinbank_high_redshift_ang_mom_delta_beta_all = table_ang_mom_evolution[4][5]
swinbank_high_redshift_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[4][6]
swinbank_high_redshift_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[4][7]
swinbank_high_redshift_ang_mom_number_rot = table_ang_mom_evolution[4][8]
swinbank_high_redshift_ang_mom_chi_rot = table_ang_mom_evolution[4][9]
swinbank_high_redshift_ang_mom_beta_rot = table_ang_mom_evolution[4][10]
swinbank_high_redshift_ang_mom_delta_beta_rot = table_ang_mom_evolution[4][11]
swinbank_high_redshift_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[4][12]
swinbank_high_redshift_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[4][13]
swinbank_high_redshift_ang_mom_number_disp = table_ang_mom_evolution[4][14]
swinbank_high_redshift_ang_mom_chi_disp = table_ang_mom_evolution[4][15]
swinbank_high_redshift_ang_mom_beta_disp = table_ang_mom_evolution[4][16]
swinbank_high_redshift_ang_mom_delta_beta_disp = table_ang_mom_evolution[4][17]
swinbank_high_redshift_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[4][18]
swinbank_high_redshift_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[4][19]

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_delta_beta_all,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_delta_beta_lower_error_all,swinbank_high_redshift_ang_mom_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_delta_beta_rot,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_delta_beta_lower_error_rot,swinbank_high_redshift_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_delta_beta_disp,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_delta_beta_lower_error_disp,swinbank_high_redshift_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)


sigma_low_redshift_redshift = table_ang_mom_evolution[1][1]
sigma_low_redshift_ang_mom_number_all = table_ang_mom_evolution[1][2]
sigma_low_redshift_ang_mom_chi_all = table_ang_mom_evolution[1][3]
sigma_low_redshift_ang_mom_beta_all = table_ang_mom_evolution[1][4]
sigma_low_redshift_ang_mom_delta_beta_all = table_ang_mom_evolution[1][5]
sigma_low_redshift_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[1][6]
sigma_low_redshift_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[1][7]
sigma_low_redshift_ang_mom_number_rot = table_ang_mom_evolution[1][8]
sigma_low_redshift_ang_mom_chi_rot = table_ang_mom_evolution[1][9]
sigma_low_redshift_ang_mom_beta_rot = table_ang_mom_evolution[1][10]
sigma_low_redshift_ang_mom_delta_beta_rot = table_ang_mom_evolution[1][11]
sigma_low_redshift_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[1][12]
sigma_low_redshift_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[1][13]
sigma_low_redshift_ang_mom_number_disp = table_ang_mom_evolution[1][14]
sigma_low_redshift_ang_mom_chi_disp = table_ang_mom_evolution[1][15]
sigma_low_redshift_ang_mom_beta_disp = table_ang_mom_evolution[1][16]
sigma_low_redshift_ang_mom_delta_beta_disp = table_ang_mom_evolution[1][17]
sigma_low_redshift_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[1][18]
sigma_low_redshift_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[1][19]

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_delta_beta_all,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_delta_beta_lower_error_all,sigma_low_redshift_ang_mom_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_delta_beta_rot,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_delta_beta_lower_error_rot,sigma_low_redshift_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_delta_beta_disp,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_delta_beta_lower_error_disp,sigma_low_redshift_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

cresci_redshift = table_ang_mom_evolution[0][1]
cresci_ang_mom_number_all = table_ang_mom_evolution[0][2]
cresci_ang_mom_chi_all = table_ang_mom_evolution[0][3]
cresci_ang_mom_beta_all = table_ang_mom_evolution[0][4]
cresci_ang_mom_delta_beta_all = table_ang_mom_evolution[0][5]
cresci_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[0][6]
cresci_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[0][7]
cresci_ang_mom_number_rot = table_ang_mom_evolution[0][2]
cresci_ang_mom_chi_rot = table_ang_mom_evolution[0][3]
cresci_ang_mom_beta_rot = table_ang_mom_evolution[0][4]
cresci_ang_mom_delta_beta_rot = table_ang_mom_evolution[0][5]
cresci_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[0][6]
cresci_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[0][7]
cresci_ang_mom_number_disp = table_ang_mom_evolution[0][14]
cresci_ang_mom_chi_disp = table_ang_mom_evolution[0][15]
cresci_ang_mom_beta_disp = table_ang_mom_evolution[0][16]
cresci_ang_mom_delta_beta_disp = table_ang_mom_evolution[0][17]
cresci_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[0][18]
cresci_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[0][19]

ax[0].errorbar(cresci_redshift,
            cresci_ang_mom_delta_beta_all,
            ecolor='darkorange',
            yerr=np.array([[cresci_ang_mom_delta_beta_lower_error_all,cresci_ang_mom_delta_beta_upper_error_all]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].errorbar(cresci_redshift,
            cresci_ang_mom_delta_beta_rot,
            ecolor='darkorange',
            yerr=np.array([[cresci_ang_mom_delta_beta_lower_error_rot,cresci_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

sigma_high_redshift_redshift = table_ang_mom_evolution[2][1]
sigma_high_redshift_ang_mom_number_all = table_ang_mom_evolution[2][2]
sigma_high_redshift_ang_mom_chi_all = table_ang_mom_evolution[2][3]
sigma_high_redshift_ang_mom_beta_all = table_ang_mom_evolution[2][4]
sigma_high_redshift_ang_mom_delta_beta_all = table_ang_mom_evolution[2][5]
sigma_high_redshift_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[2][6]
sigma_high_redshift_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[2][7]
sigma_high_redshift_ang_mom_number_rot = table_ang_mom_evolution[2][8]
sigma_high_redshift_ang_mom_chi_rot = table_ang_mom_evolution[2][9]
sigma_high_redshift_ang_mom_beta_rot = table_ang_mom_evolution[2][10]
sigma_high_redshift_ang_mom_delta_beta_rot = table_ang_mom_evolution[2][11]
sigma_high_redshift_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[2][12]
sigma_high_redshift_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[2][13]
sigma_high_redshift_ang_mom_number_disp = table_ang_mom_evolution[2][14]
sigma_high_redshift_ang_mom_chi_disp = table_ang_mom_evolution[2][15]
sigma_high_redshift_ang_mom_beta_disp = table_ang_mom_evolution[2][16]
sigma_high_redshift_ang_mom_delta_beta_disp = table_ang_mom_evolution[2][17]
sigma_high_redshift_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[2][18]
sigma_high_redshift_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[2][19]

ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_delta_beta_all,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_delta_beta_lower_error_all,sigma_high_redshift_ang_mom_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_delta_beta_rot,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_delta_beta_lower_error_rot,sigma_high_redshift_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_delta_beta_disp,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_delta_beta_lower_error_disp,sigma_high_redshift_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

amaze_redshift = table_ang_mom_evolution[8][1]
amaze_ang_mom_number_all = table_ang_mom_evolution[8][2]
amaze_ang_mom_chi_all = table_ang_mom_evolution[8][3]
amaze_ang_mom_beta_all = table_ang_mom_evolution[8][4]
amaze_ang_mom_delta_beta_all = table_ang_mom_evolution[8][5]
amaze_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[8][6]
amaze_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[8][7]
amaze_ang_mom_number_rot = table_ang_mom_evolution[8][2]
amaze_ang_mom_chi_rot = table_ang_mom_evolution[8][3]
amaze_ang_mom_beta_rot = table_ang_mom_evolution[8][4]
amaze_ang_mom_delta_beta_rot = table_ang_mom_evolution[8][5]
amaze_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[8][6]
amaze_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[8][7]
amaze_ang_mom_number_disp = table_ang_mom_evolution[8][14]
amaze_ang_mom_chi_disp = table_ang_mom_evolution[8][15]
amaze_ang_mom_beta_disp = table_ang_mom_evolution[8][16]
amaze_ang_mom_delta_beta_disp = table_ang_mom_evolution[8][17]
amaze_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[8][18]
amaze_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[8][19]

ax[0].errorbar(amaze_redshift,
            amaze_ang_mom_delta_beta_all,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_ang_mom_delta_beta_lower_error_all,amaze_ang_mom_delta_beta_upper_error_all]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].errorbar(amaze_redshift,
            amaze_ang_mom_delta_beta_rot,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_ang_mom_delta_beta_lower_error_rot,amaze_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

kds_redshift = table_ang_mom_evolution[9][1]
kds_ang_mom_number_all = table_ang_mom_evolution[9][2]
kds_ang_mom_chi_all = table_ang_mom_evolution[9][3]
kds_ang_mom_beta_all = table_ang_mom_evolution[9][4]
kds_ang_mom_delta_beta_all = table_ang_mom_evolution[9][5]
kds_ang_mom_delta_beta_lower_error_all = table_ang_mom_evolution[9][6]
kds_ang_mom_delta_beta_upper_error_all = table_ang_mom_evolution[9][7]
kds_ang_mom_number_rot = table_ang_mom_evolution[9][8]
kds_ang_mom_chi_rot = table_ang_mom_evolution[9][9]
kds_ang_mom_beta_rot = table_ang_mom_evolution[9][10]
kds_ang_mom_delta_beta_rot = table_ang_mom_evolution[9][11]
kds_ang_mom_delta_beta_lower_error_rot = table_ang_mom_evolution[9][12]
kds_ang_mom_delta_beta_upper_error_rot = table_ang_mom_evolution[9][13]
kds_ang_mom_number_disp = table_ang_mom_evolution[9][14]
kds_ang_mom_chi_disp = table_ang_mom_evolution[9][15]
kds_ang_mom_beta_disp = table_ang_mom_evolution[9][16]
kds_ang_mom_delta_beta_disp = table_ang_mom_evolution[9][17]
kds_ang_mom_delta_beta_lower_error_disp = table_ang_mom_evolution[9][18]
kds_ang_mom_delta_beta_upper_error_disp = table_ang_mom_evolution[9][19]

ax[0].errorbar(kds_redshift,
            kds_ang_mom_delta_beta_all,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_delta_beta_lower_error_all,kds_ang_mom_delta_beta_upper_error_all]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].errorbar(kds_redshift,
            kds_ang_mom_delta_beta_rot,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_delta_beta_lower_error_rot,kds_ang_mom_delta_beta_upper_error_rot]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].errorbar(kds_redshift,
            kds_ang_mom_delta_beta_disp,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_delta_beta_lower_error_disp,kds_ang_mom_delta_beta_upper_error_disp]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('upper left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='upper right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/ang_mom_evolution.png')
plt.close('all')

# make the velocity evolution figure

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'j$\boldsymbol{_{tot}}$ $\boldsymbol{\Delta \beta}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-1.0,0.6)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-1.0,0.6)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[2].set_ylim(-1.0,0.6)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

dynamo_redshift = table_ang_mom_tot_evolution[5][1]
dynamo_ang_mom_tot_number_all = table_ang_mom_tot_evolution[5][2]
dynamo_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[5][3]
dynamo_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[5][4]
dynamo_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[5][5]
dynamo_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[5][6]
dynamo_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[5][7]
dynamo_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[5][2]
dynamo_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[5][3]
dynamo_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[5][4]
dynamo_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[5][5]
dynamo_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[5][6]
dynamo_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[5][7]
dynamo_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[5][14]
dynamo_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[5][15]
dynamo_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[5][16]
dynamo_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[5][17]
dynamo_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[5][18]
dynamo_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[5][19]

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_ang_mom_tot_delta_beta_all,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_ang_mom_tot_delta_beta_lower_error_all,dynamo_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_ang_mom_tot_delta_beta_rot,
            ecolor='cornflowerblue',
            yerr=np.array([[dynamo_ang_mom_tot_delta_beta_lower_error_rot,dynamo_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_low_redshift_redshift = table_ang_mom_tot_evolution[3][1]
swinbank_low_redshift_ang_mom_tot_number_all = table_ang_mom_tot_evolution[3][2]
swinbank_low_redshift_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[3][3]
swinbank_low_redshift_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[3][4]
swinbank_low_redshift_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[3][5]
swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[3][6]
swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[3][7]
swinbank_low_redshift_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[3][8]
swinbank_low_redshift_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[3][9]
swinbank_low_redshift_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[3][10]
swinbank_low_redshift_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[3][11]
swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[3][12]
swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[3][13]
swinbank_low_redshift_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[3][14]
swinbank_low_redshift_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[3][15]
swinbank_low_redshift_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[3][16]
swinbank_low_redshift_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[3][17]
swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[3][18]
swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[3][19]

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_tot_delta_beta_all,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_all,swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_tot_delta_beta_rot,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_rot,swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_ang_mom_tot_delta_beta_disp,
            ecolor='navy',
            yerr=np.array([[swinbank_low_redshift_ang_mom_tot_delta_beta_lower_error_disp,swinbank_low_redshift_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

kross_redshift = table_ang_mom_tot_evolution[6][1]
kross_ang_mom_tot_number_all = table_ang_mom_tot_evolution[6][2]
kross_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[6][3]
kross_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[6][4]
kross_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[6][5]
kross_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[6][6]
kross_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[6][7]
kross_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[6][8]
kross_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[6][9]
kross_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[6][10]
kross_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[6][11]
kross_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[6][12]
kross_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[6][13]
kross_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[6][14]
kross_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[6][15]
kross_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[6][16]
kross_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[6][17]
kross_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[6][18]
kross_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[6][19]

p3_all = ax[0].errorbar(kross_redshift,
            kross_ang_mom_tot_delta_beta_all,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_tot_delta_beta_lower_error_all,kross_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_rot = ax[1].errorbar(kross_redshift,
            kross_ang_mom_tot_delta_beta_rot,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_tot_delta_beta_lower_error_rot,kross_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p2_disp = ax[2].errorbar(kross_redshift,
            kross_ang_mom_tot_delta_beta_disp,
            ecolor='limegreen',
            yerr=np.array([[kross_ang_mom_tot_delta_beta_lower_error_disp,kross_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='>',
            markersize=10,
            markerfacecolor='limegreen',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

massiv_redshift = table_ang_mom_tot_evolution[7][1]
massiv_ang_mom_tot_number_all = table_ang_mom_tot_evolution[7][2]
massiv_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[7][3]
massiv_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[7][4]
massiv_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[7][5]
massiv_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[7][6]
massiv_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[7][7]
massiv_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[7][8]
massiv_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[7][9]
massiv_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[7][10]
massiv_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[7][11]
massiv_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[7][12]
massiv_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[7][13]
massiv_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[7][14]
massiv_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[7][15]
massiv_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[7][16]
massiv_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[7][17]
massiv_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[7][18]
massiv_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[7][19]

p4_all = ax[0].errorbar(massiv_redshift,
            massiv_ang_mom_tot_delta_beta_all,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_tot_delta_beta_lower_error_all,massiv_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_rot = ax[1].errorbar(massiv_redshift,
            massiv_ang_mom_tot_delta_beta_rot,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_tot_delta_beta_lower_error_rot,massiv_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p3_disp = ax[2].errorbar(massiv_redshift,
            massiv_ang_mom_tot_delta_beta_disp,
            ecolor='olive',
            yerr=np.array([[massiv_ang_mom_tot_delta_beta_lower_error_disp,massiv_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='v',
            markersize=10,
            markerfacecolor='olive',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

swinbank_high_redshift_redshift = table_ang_mom_tot_evolution[4][1]
swinbank_high_redshift_ang_mom_tot_number_all = table_ang_mom_tot_evolution[4][2]
swinbank_high_redshift_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[4][3]
swinbank_high_redshift_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[4][4]
swinbank_high_redshift_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[4][5]
swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[4][6]
swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[4][7]
swinbank_high_redshift_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[4][8]
swinbank_high_redshift_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[4][9]
swinbank_high_redshift_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[4][10]
swinbank_high_redshift_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[4][11]
swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[4][12]
swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[4][13]
swinbank_high_redshift_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[4][14]
swinbank_high_redshift_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[4][15]
swinbank_high_redshift_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[4][16]
swinbank_high_redshift_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[4][17]
swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[4][18]
swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[4][19]

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_tot_delta_beta_all,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_all,swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_tot_delta_beta_rot,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_rot,swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_ang_mom_tot_delta_beta_disp,
            ecolor='darkgreen',
            yerr=np.array([[swinbank_high_redshift_ang_mom_tot_delta_beta_lower_error_disp,swinbank_high_redshift_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)


sigma_low_redshift_redshift = table_ang_mom_tot_evolution[1][1]
sigma_low_redshift_ang_mom_tot_number_all = table_ang_mom_tot_evolution[1][2]
sigma_low_redshift_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[1][3]
sigma_low_redshift_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[1][4]
sigma_low_redshift_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[1][5]
sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[1][6]
sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[1][7]
sigma_low_redshift_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[1][8]
sigma_low_redshift_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[1][9]
sigma_low_redshift_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[1][10]
sigma_low_redshift_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[1][11]
sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[1][12]
sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[1][13]
sigma_low_redshift_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[1][14]
sigma_low_redshift_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[1][15]
sigma_low_redshift_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[1][16]
sigma_low_redshift_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[1][17]
sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[1][18]
sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[1][19]

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_tot_delta_beta_all,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_all,sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_tot_delta_beta_rot,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_rot,sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_ang_mom_tot_delta_beta_disp,
            ecolor='teal',
            yerr=np.array([[sigma_low_redshift_ang_mom_tot_delta_beta_lower_error_disp,sigma_low_redshift_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

cresci_redshift = table_ang_mom_tot_evolution[0][1]
cresci_ang_mom_tot_number_all = table_ang_mom_tot_evolution[0][2]
cresci_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[0][3]
cresci_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[0][4]
cresci_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[0][5]
cresci_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[0][6]
cresci_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[0][7]
cresci_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[0][2]
cresci_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[0][3]
cresci_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[0][4]
cresci_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[0][5]
cresci_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[0][6]
cresci_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[0][7]
cresci_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[0][14]
cresci_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[0][15]
cresci_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[0][16]
cresci_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[0][17]
cresci_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[0][18]
cresci_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[0][19]

ax[0].errorbar(cresci_redshift,
            cresci_ang_mom_tot_delta_beta_all,
            ecolor='darkorange',
            yerr=np.array([[cresci_ang_mom_tot_delta_beta_lower_error_all,cresci_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].errorbar(cresci_redshift,
            cresci_ang_mom_tot_delta_beta_rot,
            ecolor='darkorange',
            yerr=np.array([[cresci_ang_mom_tot_delta_beta_lower_error_rot,cresci_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

sigma_high_redshift_redshift = table_ang_mom_tot_evolution[2][1]
sigma_high_redshift_ang_mom_tot_number_all = table_ang_mom_tot_evolution[2][2]
sigma_high_redshift_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[2][3]
sigma_high_redshift_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[2][4]
sigma_high_redshift_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[2][5]
sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[2][6]
sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[2][7]
sigma_high_redshift_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[2][8]
sigma_high_redshift_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[2][9]
sigma_high_redshift_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[2][10]
sigma_high_redshift_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[2][11]
sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[2][12]
sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[2][13]
sigma_high_redshift_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[2][14]
sigma_high_redshift_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[2][15]
sigma_high_redshift_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[2][16]
sigma_high_redshift_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[2][17]
sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[2][18]
sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[2][19]

ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_tot_delta_beta_all,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_all,sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_tot_delta_beta_rot,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_rot,sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_ang_mom_tot_delta_beta_disp,
            ecolor='red',
            yerr=np.array([[sigma_high_redshift_ang_mom_tot_delta_beta_lower_error_disp,sigma_high_redshift_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

amaze_redshift = table_ang_mom_tot_evolution[8][1]
amaze_ang_mom_tot_number_all = table_ang_mom_tot_evolution[8][2]
amaze_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[8][3]
amaze_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[8][4]
amaze_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[8][5]
amaze_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[8][6]
amaze_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[8][7]
amaze_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[8][2]
amaze_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[8][3]
amaze_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[8][4]
amaze_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[8][5]
amaze_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[8][6]
amaze_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[8][7]
amaze_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[8][14]
amaze_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[8][15]
amaze_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[8][16]
amaze_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[8][17]
amaze_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[8][18]
amaze_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[8][19]

ax[0].errorbar(amaze_redshift,
            amaze_ang_mom_tot_delta_beta_all,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_ang_mom_tot_delta_beta_lower_error_all,amaze_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].errorbar(amaze_redshift,
            amaze_ang_mom_tot_delta_beta_rot,
            ecolor='saddlebrown',
            yerr=np.array([[amaze_ang_mom_tot_delta_beta_lower_error_rot,amaze_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

kds_redshift = table_ang_mom_tot_evolution[9][1]
kds_ang_mom_tot_number_all = table_ang_mom_tot_evolution[9][2]
kds_ang_mom_tot_chi_all = table_ang_mom_tot_evolution[9][3]
kds_ang_mom_tot_beta_all = table_ang_mom_tot_evolution[9][4]
kds_ang_mom_tot_delta_beta_all = table_ang_mom_tot_evolution[9][5]
kds_ang_mom_tot_delta_beta_lower_error_all = table_ang_mom_tot_evolution[9][6]
kds_ang_mom_tot_delta_beta_upper_error_all = table_ang_mom_tot_evolution[9][7]
kds_ang_mom_tot_number_rot = table_ang_mom_tot_evolution[9][8]
kds_ang_mom_tot_chi_rot = table_ang_mom_tot_evolution[9][9]
kds_ang_mom_tot_beta_rot = table_ang_mom_tot_evolution[9][10]
kds_ang_mom_tot_delta_beta_rot = table_ang_mom_tot_evolution[9][11]
kds_ang_mom_tot_delta_beta_lower_error_rot = table_ang_mom_tot_evolution[9][12]
kds_ang_mom_tot_delta_beta_upper_error_rot = table_ang_mom_tot_evolution[9][13]
kds_ang_mom_tot_number_disp = table_ang_mom_tot_evolution[9][14]
kds_ang_mom_tot_chi_disp = table_ang_mom_tot_evolution[9][15]
kds_ang_mom_tot_beta_disp = table_ang_mom_tot_evolution[9][16]
kds_ang_mom_tot_delta_beta_disp = table_ang_mom_tot_evolution[9][17]
kds_ang_mom_tot_delta_beta_lower_error_disp = table_ang_mom_tot_evolution[9][18]
kds_ang_mom_tot_delta_beta_upper_error_disp = table_ang_mom_tot_evolution[9][19]

ax[0].errorbar(kds_redshift,
            kds_ang_mom_tot_delta_beta_all,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_tot_delta_beta_lower_error_all,kds_ang_mom_tot_delta_beta_upper_error_all]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].errorbar(kds_redshift,
            kds_ang_mom_tot_delta_beta_rot,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_tot_delta_beta_lower_error_rot,kds_ang_mom_tot_delta_beta_upper_error_rot]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].errorbar(kds_redshift,
            kds_ang_mom_tot_delta_beta_disp,
            ecolor='firebrick',
            yerr=np.array([[kds_ang_mom_tot_delta_beta_lower_error_disp,kds_ang_mom_tot_delta_beta_upper_error_disp]]).T,
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('upper left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='upper right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/ang_mom_tot_evolution.png')
plt.close('all')

# DYNAMICAL MASS WITHOUT SIGMA PLOT

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'M$\boldsymbol{_{vir}}$/M$\boldsymbol{_{\star}}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-4.0,8.0)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-4.0,8.0)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[2].set_ylim(-4.0,8.0)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

dynamo_redshift = table_dynamical_mass_evolution[5][1]
dynamo_dyn_ratio_number_all = table_dynamical_mass_evolution[5][2]
dynamo_dyn_ratio_all = table_dynamical_mass_evolution[5][3]
dynamo_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[5][4]
dynamo_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[5][5]
dynamo_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[5][6]
dynamo_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[5][7]
dynamo_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[5][8]
dynamo_dyn_ratio_number_rot = table_dynamical_mass_evolution[5][2]
dynamo_dyn_ratio_rot = table_dynamical_mass_evolution[5][3]
dynamo_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[5][4]
dynamo_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[5][5]
dynamo_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[5][6]
dynamo_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[5][7]
dynamo_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[5][8]

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_dyn_ratio_all,
            ecolor='cornflowerblue',
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([dynamo_redshift-0.05,
      dynamo_redshift+0.05,
      dynamo_redshift+0.05,
      dynamo_redshift-0.05],
     [dynamo_dyn_ratio_all - dynamo_dyn_ratio_lower_error_all,
      dynamo_dyn_ratio_all - dynamo_dyn_ratio_lower_error_all,
      dynamo_dyn_ratio_all + dynamo_dyn_ratio_upper_error_all,
      dynamo_dyn_ratio_all + dynamo_dyn_ratio_upper_error_all],
     'cornflowerblue',
      edgecolor='cornflowerblue',
      lw=2,
      alpha=0.25)


p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_dyn_ratio_rot,
            ecolor='cornflowerblue',
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([dynamo_redshift-0.05,
      dynamo_redshift+0.05,
      dynamo_redshift+0.05,
      dynamo_redshift-0.05],
     [dynamo_dyn_ratio_rot - dynamo_dyn_ratio_lower_error_rot,
      dynamo_dyn_ratio_rot - dynamo_dyn_ratio_lower_error_rot,
      dynamo_dyn_ratio_rot + dynamo_dyn_ratio_upper_error_rot,
      dynamo_dyn_ratio_rot + dynamo_dyn_ratio_upper_error_rot],
     'cornflowerblue',
      edgecolor='cornflowerblue',
      lw=2,
      alpha=0.25)

swinbank_low_redshift_redshift = table_dynamical_mass_evolution[3][1]
swinbank_low_redshift_dyn_ratio_number_all = table_dynamical_mass_evolution[3][2]
swinbank_low_redshift_dyn_ratio_all = table_dynamical_mass_evolution[3][3]
swinbank_low_redshift_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[3][4]
swinbank_low_redshift_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[3][5]
swinbank_low_redshift_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[3][6]
swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[3][7]
swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[3][8]
swinbank_low_redshift_dyn_ratio_number_rot = table_dynamical_mass_evolution[3][9]
swinbank_low_redshift_dyn_ratio_rot = table_dynamical_mass_evolution[3][10]
swinbank_low_redshift_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[3][11]
swinbank_low_redshift_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[3][12]
swinbank_low_redshift_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[3][13]
swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[3][14]
swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[3][15]
swinbank_low_redshift_dyn_ratio_number_disp = table_dynamical_mass_evolution[3][16]
swinbank_low_redshift_dyn_ratio_disp = table_dynamical_mass_evolution[3][17]
swinbank_low_redshift_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[3][18]
swinbank_low_redshift_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[3][19]
swinbank_low_redshift_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[3][20]
swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[3][21]
swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[3][22]

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_all,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_all - swinbank_low_redshift_dyn_ratio_lower_error_all,
      swinbank_low_redshift_dyn_ratio_all - swinbank_low_redshift_dyn_ratio_lower_error_all,
      swinbank_low_redshift_dyn_ratio_all + swinbank_low_redshift_dyn_ratio_upper_error_all,
      swinbank_low_redshift_dyn_ratio_all + swinbank_low_redshift_dyn_ratio_upper_error_all],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_rot,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_rot - swinbank_low_redshift_dyn_ratio_lower_error_rot,
      swinbank_low_redshift_dyn_ratio_rot - swinbank_low_redshift_dyn_ratio_lower_error_rot,
      swinbank_low_redshift_dyn_ratio_rot + swinbank_low_redshift_dyn_ratio_upper_error_rot,
      swinbank_low_redshift_dyn_ratio_rot + swinbank_low_redshift_dyn_ratio_upper_error_rot],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_disp,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_disp - swinbank_low_redshift_dyn_ratio_lower_error_disp,
      swinbank_low_redshift_dyn_ratio_disp - swinbank_low_redshift_dyn_ratio_lower_error_disp,
      swinbank_low_redshift_dyn_ratio_disp + swinbank_low_redshift_dyn_ratio_upper_error_disp,
      swinbank_low_redshift_dyn_ratio_disp + swinbank_low_redshift_dyn_ratio_upper_error_disp],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

kross_redshift = table_dynamical_mass_evolution[6][1]
kross_dyn_ratio_number_all = table_dynamical_mass_evolution[6][2]
kross_dyn_ratio_all = table_dynamical_mass_evolution[6][3]
kross_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[6][4]
kross_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[6][5]
kross_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[6][6]
kross_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[6][7]
kross_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[6][8]
kross_dyn_ratio_number_rot = table_dynamical_mass_evolution[6][9]
kross_dyn_ratio_rot = table_dynamical_mass_evolution[6][10]
kross_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[6][11]
kross_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[6][12]
kross_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[6][13]
kross_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[6][14]
kross_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[6][15]
kross_dyn_ratio_number_disp = table_dynamical_mass_evolution[6][16]
kross_dyn_ratio_disp = table_dynamical_mass_evolution[6][17]
kross_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[6][18]
kross_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[6][19]
kross_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[6][20]
kross_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[6][21]
kross_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[6][22]

p2_all = ax[0].errorbar(kross_redshift,
            kross_dyn_ratio_all,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_all - kross_dyn_ratio_lower_error_all,
      kross_dyn_ratio_all - kross_dyn_ratio_lower_error_all,
      kross_dyn_ratio_all + kross_dyn_ratio_upper_error_all,
      kross_dyn_ratio_all + kross_dyn_ratio_upper_error_all],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(kross_redshift,
            kross_dyn_ratio_rot,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_rot - kross_dyn_ratio_lower_error_rot,
      kross_dyn_ratio_rot - kross_dyn_ratio_lower_error_rot,
      kross_dyn_ratio_rot + kross_dyn_ratio_upper_error_rot,
      kross_dyn_ratio_rot + kross_dyn_ratio_upper_error_rot],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(kross_redshift,
            kross_dyn_ratio_disp,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_disp - kross_dyn_ratio_lower_error_disp,
      kross_dyn_ratio_disp - kross_dyn_ratio_lower_error_disp,
      kross_dyn_ratio_disp + kross_dyn_ratio_upper_error_disp,
      kross_dyn_ratio_disp + kross_dyn_ratio_upper_error_disp],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

massiv_redshift = table_dynamical_mass_evolution[7][1]
massiv_dyn_ratio_number_all = table_dynamical_mass_evolution[7][2]
massiv_dyn_ratio_all = table_dynamical_mass_evolution[7][3]
massiv_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[7][4]
massiv_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[7][5]
massiv_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[7][6]
massiv_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[7][7]
massiv_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[7][8]
massiv_dyn_ratio_number_rot = table_dynamical_mass_evolution[7][9]
massiv_dyn_ratio_rot = table_dynamical_mass_evolution[7][10]
massiv_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[7][11]
massiv_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[7][12]
massiv_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[7][13]
massiv_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[7][14]
massiv_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[7][15]
massiv_dyn_ratio_number_disp = table_dynamical_mass_evolution[7][16]
massiv_dyn_ratio_disp = table_dynamical_mass_evolution[7][17]
massiv_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[7][18]
massiv_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[7][19]
massiv_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[7][20]
massiv_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[7][21]
massiv_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[7][22]

p2_all = ax[0].errorbar(massiv_redshift,
            massiv_dyn_ratio_all,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_all - massiv_dyn_ratio_lower_error_all,
      massiv_dyn_ratio_all - massiv_dyn_ratio_lower_error_all,
      massiv_dyn_ratio_all + massiv_dyn_ratio_upper_error_all,
      massiv_dyn_ratio_all + massiv_dyn_ratio_upper_error_all],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(massiv_redshift,
            massiv_dyn_ratio_rot,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_rot - massiv_dyn_ratio_lower_error_rot,
      massiv_dyn_ratio_rot - massiv_dyn_ratio_lower_error_rot,
      massiv_dyn_ratio_rot + massiv_dyn_ratio_upper_error_rot,
      massiv_dyn_ratio_rot + massiv_dyn_ratio_upper_error_rot],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(massiv_redshift,
            massiv_dyn_ratio_disp,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_disp - massiv_dyn_ratio_lower_error_disp,
      massiv_dyn_ratio_disp - massiv_dyn_ratio_lower_error_disp,
      massiv_dyn_ratio_disp + massiv_dyn_ratio_upper_error_disp,
      massiv_dyn_ratio_disp + massiv_dyn_ratio_upper_error_disp],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

swinbank_high_redshift_redshift = table_dynamical_mass_evolution[4][1]
swinbank_high_redshift_dyn_ratio_number_all = table_dynamical_mass_evolution[4][2]
swinbank_high_redshift_dyn_ratio_all = table_dynamical_mass_evolution[4][3]
swinbank_high_redshift_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[4][4]
swinbank_high_redshift_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[4][5]
swinbank_high_redshift_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[4][6]
swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[4][7]
swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[4][8]
swinbank_high_redshift_dyn_ratio_number_rot = table_dynamical_mass_evolution[4][9]
swinbank_high_redshift_dyn_ratio_rot = table_dynamical_mass_evolution[4][10]
swinbank_high_redshift_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[4][11]
swinbank_high_redshift_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[4][12]
swinbank_high_redshift_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[4][13]
swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[4][14]
swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[4][15]
swinbank_high_redshift_dyn_ratio_number_disp = table_dynamical_mass_evolution[4][16]
swinbank_high_redshift_dyn_ratio_disp = table_dynamical_mass_evolution[4][17]
swinbank_high_redshift_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[4][18]
swinbank_high_redshift_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[4][19]
swinbank_high_redshift_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[4][20]
swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[4][21]
swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[4][22]

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_all,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_all - swinbank_high_redshift_dyn_ratio_lower_error_all,
      swinbank_high_redshift_dyn_ratio_all - swinbank_high_redshift_dyn_ratio_lower_error_all,
      swinbank_high_redshift_dyn_ratio_all + swinbank_high_redshift_dyn_ratio_upper_error_all,
      swinbank_high_redshift_dyn_ratio_all + swinbank_high_redshift_dyn_ratio_upper_error_all],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_rot,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_rot - swinbank_high_redshift_dyn_ratio_lower_error_rot,
      swinbank_high_redshift_dyn_ratio_rot - swinbank_high_redshift_dyn_ratio_lower_error_rot,
      swinbank_high_redshift_dyn_ratio_rot + swinbank_high_redshift_dyn_ratio_upper_error_rot,
      swinbank_high_redshift_dyn_ratio_rot + swinbank_high_redshift_dyn_ratio_upper_error_rot],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_disp,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_disp - swinbank_high_redshift_dyn_ratio_lower_error_disp,
      swinbank_high_redshift_dyn_ratio_disp - swinbank_high_redshift_dyn_ratio_lower_error_disp,
      swinbank_high_redshift_dyn_ratio_disp + swinbank_high_redshift_dyn_ratio_upper_error_disp,
      swinbank_high_redshift_dyn_ratio_disp + swinbank_high_redshift_dyn_ratio_upper_error_disp],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

sigma_low_redshift_redshift = table_dynamical_mass_evolution[1][1]
sigma_low_redshift_dyn_ratio_number_all = table_dynamical_mass_evolution[1][2]
sigma_low_redshift_dyn_ratio_all = table_dynamical_mass_evolution[1][3]
sigma_low_redshift_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[1][4]
sigma_low_redshift_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[1][5]
sigma_low_redshift_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[1][6]
sigma_low_redshift_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[1][7]
sigma_low_redshift_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[1][8]
sigma_low_redshift_dyn_ratio_number_rot = table_dynamical_mass_evolution[1][9]
sigma_low_redshift_dyn_ratio_rot = table_dynamical_mass_evolution[1][10]
sigma_low_redshift_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[1][11]
sigma_low_redshift_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[1][12]
sigma_low_redshift_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[1][13]
sigma_low_redshift_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[1][14]
sigma_low_redshift_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[1][15]
sigma_low_redshift_dyn_ratio_number_disp = table_dynamical_mass_evolution[1][16]
sigma_low_redshift_dyn_ratio_disp = table_dynamical_mass_evolution[1][17]
sigma_low_redshift_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[1][18]
sigma_low_redshift_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[1][19]
sigma_low_redshift_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[1][20]
sigma_low_redshift_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[1][21]
sigma_low_redshift_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[1][22]

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_all,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[0].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_all - sigma_low_redshift_dyn_ratio_lower_error_all,
      sigma_low_redshift_dyn_ratio_all - sigma_low_redshift_dyn_ratio_lower_error_all,
      sigma_low_redshift_dyn_ratio_all + sigma_low_redshift_dyn_ratio_upper_error_all,
      sigma_low_redshift_dyn_ratio_all + sigma_low_redshift_dyn_ratio_upper_error_all],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_rot,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_rot - sigma_low_redshift_dyn_ratio_lower_error_rot,
      sigma_low_redshift_dyn_ratio_rot - sigma_low_redshift_dyn_ratio_lower_error_rot,
      sigma_low_redshift_dyn_ratio_rot + sigma_low_redshift_dyn_ratio_upper_error_rot,
      sigma_low_redshift_dyn_ratio_rot + sigma_low_redshift_dyn_ratio_upper_error_rot],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_disp,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_disp - sigma_low_redshift_dyn_ratio_lower_error_disp,
      sigma_low_redshift_dyn_ratio_disp - sigma_low_redshift_dyn_ratio_lower_error_disp,
      sigma_low_redshift_dyn_ratio_disp + sigma_low_redshift_dyn_ratio_upper_error_disp,
      sigma_low_redshift_dyn_ratio_disp + sigma_low_redshift_dyn_ratio_upper_error_disp],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

cresci_redshift = table_dynamical_mass_evolution[0][1]
cresci_dyn_ratio_number_all = table_dynamical_mass_evolution[0][2]
cresci_dyn_ratio_all = table_dynamical_mass_evolution[0][3]
cresci_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[0][4]
cresci_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[0][5]
cresci_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[0][6]
cresci_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[0][7]
cresci_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[0][8]
cresci_dyn_ratio_number_rot = table_dynamical_mass_evolution[0][2]
cresci_dyn_ratio_rot = table_dynamical_mass_evolution[0][3]
cresci_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[0][4]
cresci_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[0][5]
cresci_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[0][6]
cresci_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[0][7]
cresci_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[0][8]

ax[0].errorbar(cresci_redshift,
            cresci_dyn_ratio_all,
            ecolor='darkorange',
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[0].fill([cresci_redshift-0.05,
      cresci_redshift+0.05,
      cresci_redshift+0.05,
      cresci_redshift-0.05],
     [cresci_dyn_ratio_all - cresci_dyn_ratio_lower_error_all,
      cresci_dyn_ratio_all - cresci_dyn_ratio_lower_error_all,
      cresci_dyn_ratio_all + cresci_dyn_ratio_upper_error_all,
      cresci_dyn_ratio_all + cresci_dyn_ratio_upper_error_all],
     'darkorange',
      edgecolor='darkorange',
      lw=2,
      alpha=0.25)

ax[1].errorbar(cresci_redshift,
            cresci_dyn_ratio_rot,
            ecolor='darkorange',
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].fill([cresci_redshift-0.05,
      cresci_redshift+0.05,
      cresci_redshift+0.05,
      cresci_redshift-0.05],
     [cresci_dyn_ratio_rot - cresci_dyn_ratio_lower_error_rot,
      cresci_dyn_ratio_rot - cresci_dyn_ratio_lower_error_rot,
      cresci_dyn_ratio_rot + cresci_dyn_ratio_upper_error_rot,
      cresci_dyn_ratio_rot + cresci_dyn_ratio_upper_error_rot],
     'darkorange',
      edgecolor='darkorange',
      lw=2,
      alpha=0.25)

sigma_high_redshift_redshift = table_dynamical_mass_evolution[2][1]
sigma_high_redshift_dyn_ratio_number_all = table_dynamical_mass_evolution[2][2]
sigma_high_redshift_dyn_ratio_all = table_dynamical_mass_evolution[2][3]
sigma_high_redshift_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[2][4]
sigma_high_redshift_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[2][5]
sigma_high_redshift_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[2][6]
sigma_high_redshift_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[2][7]
sigma_high_redshift_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[2][8]
sigma_high_redshift_dyn_ratio_number_rot = table_dynamical_mass_evolution[2][9]
sigma_high_redshift_dyn_ratio_rot = table_dynamical_mass_evolution[2][10]
sigma_high_redshift_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[2][11]
sigma_high_redshift_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[2][12]
sigma_high_redshift_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[2][13]
sigma_high_redshift_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[2][14]
sigma_high_redshift_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[2][15]
sigma_high_redshift_dyn_ratio_number_disp = table_dynamical_mass_evolution[2][16]
sigma_high_redshift_dyn_ratio_disp = table_dynamical_mass_evolution[2][17]
sigma_high_redshift_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[2][18]
sigma_high_redshift_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[2][19]
sigma_high_redshift_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[2][20]
sigma_high_redshift_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[2][21]
sigma_high_redshift_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[2][22]

ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_all,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[0].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_all - sigma_high_redshift_dyn_ratio_lower_error_all,
      sigma_high_redshift_dyn_ratio_all - sigma_high_redshift_dyn_ratio_lower_error_all,
      sigma_high_redshift_dyn_ratio_all + sigma_high_redshift_dyn_ratio_upper_error_all,
      sigma_high_redshift_dyn_ratio_all + sigma_high_redshift_dyn_ratio_upper_error_all],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_rot,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_rot - sigma_high_redshift_dyn_ratio_lower_error_rot,
      sigma_high_redshift_dyn_ratio_rot - sigma_high_redshift_dyn_ratio_lower_error_rot,
      sigma_high_redshift_dyn_ratio_rot + sigma_high_redshift_dyn_ratio_upper_error_rot,
      sigma_high_redshift_dyn_ratio_rot + sigma_high_redshift_dyn_ratio_upper_error_rot],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_disp,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_disp - sigma_high_redshift_dyn_ratio_lower_error_disp,
      sigma_high_redshift_dyn_ratio_disp - sigma_high_redshift_dyn_ratio_lower_error_disp,
      sigma_high_redshift_dyn_ratio_disp + sigma_high_redshift_dyn_ratio_upper_error_disp,
      sigma_high_redshift_dyn_ratio_disp + sigma_high_redshift_dyn_ratio_upper_error_disp],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

amaze_redshift = table_dynamical_mass_evolution[8][1]
amaze_dyn_ratio_number_all = table_dynamical_mass_evolution[8][2]
amaze_dyn_ratio_all = table_dynamical_mass_evolution[8][3]
amaze_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[8][4]
amaze_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[8][5]
amaze_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[8][6]
amaze_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[8][7]
amaze_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[8][8]
amaze_dyn_ratio_number_rot = table_dynamical_mass_evolution[8][2]
amaze_dyn_ratio_rot = table_dynamical_mass_evolution[8][3]
amaze_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[8][4]
amaze_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[8][5]
amaze_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[8][6]
amaze_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[8][7]
amaze_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[8][8]

ax[0].errorbar(amaze_redshift,
            amaze_dyn_ratio_all,
            ecolor='saddlebrown',
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[0].fill([amaze_redshift-0.05,
      amaze_redshift+0.05,
      amaze_redshift+0.05,
      amaze_redshift-0.05],
     [amaze_dyn_ratio_all - amaze_dyn_ratio_lower_error_all,
      amaze_dyn_ratio_all - amaze_dyn_ratio_lower_error_all,
      amaze_dyn_ratio_all + amaze_dyn_ratio_upper_error_all,
      amaze_dyn_ratio_all + amaze_dyn_ratio_upper_error_all],
     'saddlebrown',
      edgecolor='saddlebrown',
      lw=2,
      alpha=0.25)

ax[1].errorbar(amaze_redshift,
            amaze_dyn_ratio_rot,
            ecolor='saddlebrown',
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].fill([amaze_redshift-0.05,
      amaze_redshift+0.05,
      amaze_redshift+0.05,
      amaze_redshift-0.05],
     [amaze_dyn_ratio_rot - amaze_dyn_ratio_lower_error_rot,
      amaze_dyn_ratio_rot - amaze_dyn_ratio_lower_error_rot,
      amaze_dyn_ratio_rot + amaze_dyn_ratio_upper_error_rot,
      amaze_dyn_ratio_rot + amaze_dyn_ratio_upper_error_rot],
     'saddlebrown',
      edgecolor='saddlebrown',
      lw=2,
      alpha=0.25)

kds_redshift = table_dynamical_mass_evolution[9][1]
kds_dyn_ratio_number_all = table_dynamical_mass_evolution[9][2]
kds_dyn_ratio_all = table_dynamical_mass_evolution[9][3]
kds_dyn_ratio_lower_error_all = table_dynamical_mass_evolution[9][4]
kds_dyn_ratio_upper_error_all = table_dynamical_mass_evolution[9][5]
kds_dyn_ratio_with_sigma_all = table_dynamical_mass_evolution[9][6]
kds_dyn_ratio_with_sigma_lower_error_all = table_dynamical_mass_evolution[9][7]
kds_dyn_ratio_with_sigma_upper_error_all = table_dynamical_mass_evolution[9][8]
kds_dyn_ratio_number_rot = table_dynamical_mass_evolution[9][9]
kds_dyn_ratio_rot = table_dynamical_mass_evolution[9][10]
kds_dyn_ratio_lower_error_rot = table_dynamical_mass_evolution[9][11]
kds_dyn_ratio_upper_error_rot = table_dynamical_mass_evolution[9][12]
kds_dyn_ratio_with_sigma_rot = table_dynamical_mass_evolution[9][13]
kds_dyn_ratio_with_sigma_lower_error_rot = table_dynamical_mass_evolution[9][14]
kds_dyn_ratio_with_sigma_upper_error_rot = table_dynamical_mass_evolution[9][15]
kds_dyn_ratio_number_disp = table_dynamical_mass_evolution[9][16]
kds_dyn_ratio_disp = table_dynamical_mass_evolution[9][17]
kds_dyn_ratio_lower_error_disp = table_dynamical_mass_evolution[9][18]
kds_dyn_ratio_upper_error_disp = table_dynamical_mass_evolution[9][19]
kds_dyn_ratio_with_sigma_disp = table_dynamical_mass_evolution[9][20]
kds_dyn_ratio_with_sigma_lower_error_disp = table_dynamical_mass_evolution[9][21]
kds_dyn_ratio_with_sigma_upper_error_disp = table_dynamical_mass_evolution[9][22]

ax[0].errorbar(kds_redshift,
            kds_dyn_ratio_all,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[0].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_all - kds_dyn_ratio_lower_error_all,
      kds_dyn_ratio_all - kds_dyn_ratio_lower_error_all,
      kds_dyn_ratio_all + kds_dyn_ratio_upper_error_all,
      kds_dyn_ratio_all + kds_dyn_ratio_upper_error_all],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

ax[1].errorbar(kds_redshift,
            kds_dyn_ratio_rot,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_rot - kds_dyn_ratio_lower_error_rot,
      kds_dyn_ratio_rot - kds_dyn_ratio_lower_error_rot,
      kds_dyn_ratio_rot + kds_dyn_ratio_upper_error_rot,
      kds_dyn_ratio_rot + kds_dyn_ratio_upper_error_rot],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

ax[2].errorbar(kds_redshift,
            kds_dyn_ratio_disp,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_disp - kds_dyn_ratio_lower_error_disp,
      kds_dyn_ratio_disp - kds_dyn_ratio_lower_error_disp,
      kds_dyn_ratio_disp + kds_dyn_ratio_upper_error_disp,
      kds_dyn_ratio_disp + kds_dyn_ratio_upper_error_disp],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/dyn_mass_evolution.png')
plt.close('all')

# make the dynamical mass evolution plots

lines = {'linestyle': 'None'}
plt.rc('lines', **lines)

fig, ax = plt.subplots(1, 3, sharey=True, figsize=(18,7))


ax[0].set_ylabel(r'M$\boldsymbol{_{vir}}$/M$\boldsymbol{_{\star}}$',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

ax[0].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[0].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[0].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[0].spines.itervalues()]
ax[0].set_xlim(-0.3,3.9)
ax[0].set_ylim(-4.0,8.0)
ax[0].minorticks_on()

xa = ax[0].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

ax[1].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[1].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[1].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[1].spines.itervalues()]
ax[1].set_xlim(-0.3,3.9)
ax[1].set_ylim(-4.0,8.0)
ax[1].minorticks_on()


xa = ax[1].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))


ax[2].set_xlabel(r'redshift',
              fontsize=30,
              fontweight='bold',
              labelpad=15)

# tick parameters 
ax[2].tick_params(axis='both',
                   which='major',
                   labelsize=26,
                   length=12,
                   width=4)
ax[2].tick_params(axis='both',
                   which='minor',
                   labelsize=26,
                   length=6,
                   width=4)

[i.set_linewidth(4.0) for i in ax[2].spines.itervalues()]
ax[2].set_xlim(-0.3,3.9)
ax[2].set_ylim(-4.0,8.0)
ax[2].minorticks_on()


xa = ax[2].get_xaxis()
xa.set_major_locator(MaxNLocator(integer=True))

p1_all = ax[0].errorbar(dynamo_redshift,
            dynamo_dyn_ratio_with_sigma_all,
            ecolor='cornflowerblue',
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([dynamo_redshift-0.05,
      dynamo_redshift+0.05,
      dynamo_redshift+0.05,
      dynamo_redshift-0.05],
     [dynamo_dyn_ratio_with_sigma_all - dynamo_dyn_ratio_with_sigma_lower_error_all,
      dynamo_dyn_ratio_with_sigma_all - dynamo_dyn_ratio_with_sigma_lower_error_all,
      dynamo_dyn_ratio_with_sigma_all + dynamo_dyn_ratio_with_sigma_upper_error_all,
      dynamo_dyn_ratio_with_sigma_all + dynamo_dyn_ratio_with_sigma_upper_error_all],
     'cornflowerblue',
      edgecolor='cornflowerblue',
      lw=2,
      alpha=0.25)


p1_rot = ax[1].errorbar(dynamo_redshift,
            dynamo_dyn_ratio_with_sigma_rot,
            ecolor='cornflowerblue',
            marker='p',
            markersize=10,
            markerfacecolor='cornflowerblue',
            markeredgecolor='cornflowerblue',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([dynamo_redshift-0.05,
      dynamo_redshift+0.05,
      dynamo_redshift+0.05,
      dynamo_redshift-0.05],
     [dynamo_dyn_ratio_with_sigma_rot - dynamo_dyn_ratio_with_sigma_lower_error_rot,
      dynamo_dyn_ratio_with_sigma_rot - dynamo_dyn_ratio_with_sigma_lower_error_rot,
      dynamo_dyn_ratio_with_sigma_rot + dynamo_dyn_ratio_with_sigma_upper_error_rot,
      dynamo_dyn_ratio_with_sigma_rot + dynamo_dyn_ratio_with_sigma_upper_error_rot],
     'cornflowerblue',
      edgecolor='cornflowerblue',
      lw=2,
      alpha=0.25)

p2_all = ax[0].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_with_sigma_all,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_with_sigma_all - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_all,
      swinbank_low_redshift_dyn_ratio_with_sigma_all - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_all,
      swinbank_low_redshift_dyn_ratio_with_sigma_all + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_all,
      swinbank_low_redshift_dyn_ratio_with_sigma_all + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_all],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_with_sigma_rot,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_with_sigma_rot - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_rot,
      swinbank_low_redshift_dyn_ratio_with_sigma_rot - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_rot,
      swinbank_low_redshift_dyn_ratio_with_sigma_rot + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_rot,
      swinbank_low_redshift_dyn_ratio_with_sigma_rot + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_rot],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(swinbank_low_redshift_redshift,
            swinbank_low_redshift_dyn_ratio_with_sigma_disp,
            ecolor='navy',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='navy',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([swinbank_low_redshift_redshift-0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift+0.05,
      swinbank_low_redshift_redshift-0.05],
     [swinbank_low_redshift_dyn_ratio_with_sigma_disp - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_disp,
      swinbank_low_redshift_dyn_ratio_with_sigma_disp - swinbank_low_redshift_dyn_ratio_with_sigma_lower_error_disp,
      swinbank_low_redshift_dyn_ratio_with_sigma_disp + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_disp,
      swinbank_low_redshift_dyn_ratio_with_sigma_disp + swinbank_low_redshift_dyn_ratio_with_sigma_upper_error_disp],
     'navy',
      edgecolor='navy',
      lw=2,
      alpha=0.25)

p2_all = ax[0].errorbar(kross_redshift,
            kross_dyn_ratio_with_sigma_all,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_with_sigma_all - kross_dyn_ratio_with_sigma_lower_error_all,
      kross_dyn_ratio_with_sigma_all - kross_dyn_ratio_with_sigma_lower_error_all,
      kross_dyn_ratio_with_sigma_all + kross_dyn_ratio_with_sigma_upper_error_all,
      kross_dyn_ratio_with_sigma_all + kross_dyn_ratio_with_sigma_upper_error_all],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(kross_redshift,
            kross_dyn_ratio_with_sigma_rot,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_with_sigma_rot - kross_dyn_ratio_with_sigma_lower_error_rot,
      kross_dyn_ratio_with_sigma_rot - kross_dyn_ratio_with_sigma_lower_error_rot,
      kross_dyn_ratio_with_sigma_rot + kross_dyn_ratio_with_sigma_upper_error_rot,
      kross_dyn_ratio_with_sigma_rot + kross_dyn_ratio_with_sigma_upper_error_rot],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(kross_redshift,
            kross_dyn_ratio_with_sigma_disp,
            ecolor='limegreen',
            marker='>',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='limegreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([kross_redshift-0.05,
      kross_redshift+0.05,
      kross_redshift+0.05,
      kross_redshift-0.05],
     [kross_dyn_ratio_with_sigma_disp - kross_dyn_ratio_with_sigma_lower_error_disp,
      kross_dyn_ratio_with_sigma_disp - kross_dyn_ratio_with_sigma_lower_error_disp,
      kross_dyn_ratio_with_sigma_disp + kross_dyn_ratio_with_sigma_upper_error_disp,
      kross_dyn_ratio_with_sigma_disp + kross_dyn_ratio_with_sigma_upper_error_disp],
     'limegreen',
      edgecolor='limegreen',
      lw=2,
      alpha=0.25)

p2_all = ax[0].errorbar(massiv_redshift,
            massiv_dyn_ratio_with_sigma_all,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_with_sigma_all - massiv_dyn_ratio_with_sigma_lower_error_all,
      massiv_dyn_ratio_with_sigma_all - massiv_dyn_ratio_with_sigma_lower_error_all,
      massiv_dyn_ratio_with_sigma_all + massiv_dyn_ratio_with_sigma_upper_error_all,
      massiv_dyn_ratio_with_sigma_all + massiv_dyn_ratio_with_sigma_upper_error_all],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

p2_rot = ax[1].errorbar(massiv_redshift,
            massiv_dyn_ratio_with_sigma_rot,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_with_sigma_rot - massiv_dyn_ratio_with_sigma_lower_error_rot,
      massiv_dyn_ratio_with_sigma_rot - massiv_dyn_ratio_with_sigma_lower_error_rot,
      massiv_dyn_ratio_with_sigma_rot + massiv_dyn_ratio_with_sigma_upper_error_rot,
      massiv_dyn_ratio_with_sigma_rot + massiv_dyn_ratio_with_sigma_upper_error_rot],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

p1_disp = ax[2].errorbar(massiv_redshift,
            massiv_dyn_ratio_with_sigma_disp,
            ecolor='olive',
            marker='v',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='olive',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([massiv_redshift-0.05,
      massiv_redshift+0.05,
      massiv_redshift+0.05,
      massiv_redshift-0.05],
     [massiv_dyn_ratio_with_sigma_disp - massiv_dyn_ratio_with_sigma_lower_error_disp,
      massiv_dyn_ratio_with_sigma_disp - massiv_dyn_ratio_with_sigma_lower_error_disp,
      massiv_dyn_ratio_with_sigma_disp + massiv_dyn_ratio_with_sigma_upper_error_disp,
      massiv_dyn_ratio_with_sigma_disp + massiv_dyn_ratio_with_sigma_upper_error_disp],
     'olive',
      edgecolor='olive',
      lw=2,
      alpha=0.25)

p5_all = ax[0].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_with_sigma_all,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[0].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_with_sigma_all - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_all,
      swinbank_high_redshift_dyn_ratio_with_sigma_all - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_all,
      swinbank_high_redshift_dyn_ratio_with_sigma_all + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_all,
      swinbank_high_redshift_dyn_ratio_with_sigma_all + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_all],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

p5_rot = ax[1].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_with_sigma_rot,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[1].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_with_sigma_rot - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_rot,
      swinbank_high_redshift_dyn_ratio_with_sigma_rot - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_rot,
      swinbank_high_redshift_dyn_ratio_with_sigma_rot + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_rot,
      swinbank_high_redshift_dyn_ratio_with_sigma_rot + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_rot],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

p4_disp = ax[2].errorbar(swinbank_high_redshift_redshift,
            swinbank_high_redshift_dyn_ratio_with_sigma_disp,
            ecolor='darkgreen',
            marker='^',
            markersize=10,
            markerfacecolor='none',
            markeredgecolor='darkgreen',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2)

ax[2].fill([swinbank_high_redshift_redshift-0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift+0.05,
      swinbank_high_redshift_redshift-0.05],
     [swinbank_high_redshift_dyn_ratio_with_sigma_disp - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_disp,
      swinbank_high_redshift_dyn_ratio_with_sigma_disp - swinbank_high_redshift_dyn_ratio_with_sigma_lower_error_disp,
      swinbank_high_redshift_dyn_ratio_with_sigma_disp + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_disp,
      swinbank_high_redshift_dyn_ratio_with_sigma_disp + swinbank_high_redshift_dyn_ratio_with_sigma_upper_error_disp],
     'darkgreen',
      edgecolor='darkgreen',
      lw=2,
      alpha=0.25)

ax[0].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_with_sigma_all,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[0].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_with_sigma_all - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_all,
      sigma_low_redshift_dyn_ratio_with_sigma_all - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_all,
      sigma_low_redshift_dyn_ratio_with_sigma_all + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_all,
      sigma_low_redshift_dyn_ratio_with_sigma_all + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_all],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

ax[1].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_with_sigma_rot,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[1].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_with_sigma_rot - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_rot,
      sigma_low_redshift_dyn_ratio_with_sigma_rot - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_rot,
      sigma_low_redshift_dyn_ratio_with_sigma_rot + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_rot,
      sigma_low_redshift_dyn_ratio_with_sigma_rot + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_rot],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

ax[2].errorbar(sigma_low_redshift_redshift,
            sigma_low_redshift_dyn_ratio_with_sigma_disp,
            ecolor='teal',
            marker='H',
            markersize=10,
            markerfacecolor='teal',
            markeredgecolor='teal',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.1 - SIGMA (SIMONS+16)}')

ax[2].fill([sigma_low_redshift_redshift-0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift+0.05,
      sigma_low_redshift_redshift-0.05],
     [sigma_low_redshift_dyn_ratio_with_sigma_disp - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_disp,
      sigma_low_redshift_dyn_ratio_with_sigma_disp - sigma_low_redshift_dyn_ratio_with_sigma_lower_error_disp,
      sigma_low_redshift_dyn_ratio_with_sigma_disp + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_disp,
      sigma_low_redshift_dyn_ratio_with_sigma_disp + sigma_low_redshift_dyn_ratio_with_sigma_upper_error_disp],
     'teal',
      edgecolor='teal',
      lw=2,
      alpha=0.25)

ax[0].errorbar(cresci_redshift,
            cresci_dyn_ratio_with_sigma_all,
            ecolor='darkorange',
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[0].fill([cresci_redshift-0.05,
      cresci_redshift+0.05,
      cresci_redshift+0.05,
      cresci_redshift-0.05],
     [cresci_dyn_ratio_with_sigma_all - cresci_dyn_ratio_with_sigma_lower_error_all,
      cresci_dyn_ratio_with_sigma_all - cresci_dyn_ratio_with_sigma_lower_error_all,
      cresci_dyn_ratio_with_sigma_all + cresci_dyn_ratio_with_sigma_upper_error_all,
      cresci_dyn_ratio_with_sigma_all + cresci_dyn_ratio_with_sigma_upper_error_all],
     'darkorange',
      edgecolor='darkorange',
      lw=2,
      alpha=0.25)

ax[1].errorbar(cresci_redshift,
            cresci_dyn_ratio_with_sigma_rot,
            ecolor='darkorange',
            marker='*',
            markersize=15,
            markerfacecolor='darkorange',
            markeredgecolor='darkorange',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.6 - SINS (Cresci+09)}')

ax[1].fill([cresci_redshift-0.05,
      cresci_redshift+0.05,
      cresci_redshift+0.05,
      cresci_redshift-0.05],
     [cresci_dyn_ratio_with_sigma_rot - cresci_dyn_ratio_with_sigma_lower_error_rot,
      cresci_dyn_ratio_with_sigma_rot - cresci_dyn_ratio_with_sigma_lower_error_rot,
      cresci_dyn_ratio_with_sigma_rot + cresci_dyn_ratio_with_sigma_upper_error_rot,
      cresci_dyn_ratio_with_sigma_rot + cresci_dyn_ratio_with_sigma_upper_error_rot],
     'darkorange',
      edgecolor='darkorange',
      lw=2,
      alpha=0.25)


ax[0].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_with_sigma_all,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[0].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_with_sigma_all - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_all,
      sigma_high_redshift_dyn_ratio_with_sigma_all - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_all,
      sigma_high_redshift_dyn_ratio_with_sigma_all + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_all,
      sigma_high_redshift_dyn_ratio_with_sigma_all + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_all],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

ax[1].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_with_sigma_rot,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[1].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_with_sigma_rot - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_rot,
      sigma_high_redshift_dyn_ratio_with_sigma_rot - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_rot,
      sigma_high_redshift_dyn_ratio_with_sigma_rot + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_rot,
      sigma_high_redshift_dyn_ratio_with_sigma_rot + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_rot],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

ax[2].errorbar(sigma_high_redshift_redshift,
            sigma_high_redshift_dyn_ratio_with_sigma_disp,
            ecolor='red',
            marker='H',
            markersize=10,
            markerfacecolor='red',
            markeredgecolor='red',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - SIGMA (SIMONS+16)}')

ax[2].fill([sigma_high_redshift_redshift-0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift+0.05,
      sigma_high_redshift_redshift-0.05],
     [sigma_high_redshift_dyn_ratio_with_sigma_disp - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_disp,
      sigma_high_redshift_dyn_ratio_with_sigma_disp - sigma_high_redshift_dyn_ratio_with_sigma_lower_error_disp,
      sigma_high_redshift_dyn_ratio_with_sigma_disp + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_disp,
      sigma_high_redshift_dyn_ratio_with_sigma_disp + sigma_high_redshift_dyn_ratio_with_sigma_upper_error_disp],
     'red',
      edgecolor='red',
      lw=2,
      alpha=0.25)

ax[0].errorbar(amaze_redshift,
            amaze_dyn_ratio_with_sigma_all,
            ecolor='saddlebrown',
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[0].fill([amaze_redshift-0.05,
      amaze_redshift+0.05,
      amaze_redshift+0.05,
      amaze_redshift-0.05],
     [amaze_dyn_ratio_with_sigma_all - amaze_dyn_ratio_with_sigma_lower_error_all,
      amaze_dyn_ratio_with_sigma_all - amaze_dyn_ratio_with_sigma_lower_error_all,
      amaze_dyn_ratio_with_sigma_all + amaze_dyn_ratio_with_sigma_upper_error_all,
      amaze_dyn_ratio_with_sigma_all + amaze_dyn_ratio_with_sigma_upper_error_all],
     'saddlebrown',
      edgecolor='saddlebrown',
      lw=2,
      alpha=0.25)

ax[1].errorbar(amaze_redshift,
            amaze_dyn_ratio_with_sigma_rot,
            ecolor='saddlebrown',
            marker='D',
            markersize=10,
            markerfacecolor='saddlebrown',
            markeredgecolor='saddlebrown',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{10.0 - AMAZE (Gnerucci+11)}')

ax[1].fill([amaze_redshift-0.05,
      amaze_redshift+0.05,
      amaze_redshift+0.05,
      amaze_redshift-0.05],
     [amaze_dyn_ratio_with_sigma_rot - amaze_dyn_ratio_with_sigma_lower_error_rot,
      amaze_dyn_ratio_with_sigma_rot - amaze_dyn_ratio_with_sigma_lower_error_rot,
      amaze_dyn_ratio_with_sigma_rot + amaze_dyn_ratio_with_sigma_upper_error_rot,
      amaze_dyn_ratio_with_sigma_rot + amaze_dyn_ratio_with_sigma_upper_error_rot],
     'saddlebrown',
      edgecolor='saddlebrown',
      lw=2,
      alpha=0.25)

ax[0].errorbar(kds_redshift,
            kds_dyn_ratio_with_sigma_all,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[0].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_with_sigma_all - kds_dyn_ratio_with_sigma_lower_error_all,
      kds_dyn_ratio_with_sigma_all - kds_dyn_ratio_with_sigma_lower_error_all,
      kds_dyn_ratio_with_sigma_all + kds_dyn_ratio_with_sigma_upper_error_all,
      kds_dyn_ratio_with_sigma_all + kds_dyn_ratio_with_sigma_upper_error_all],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

ax[1].errorbar(kds_redshift,
            kds_dyn_ratio_with_sigma_rot,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[1].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_with_sigma_rot - kds_dyn_ratio_with_sigma_lower_error_rot,
      kds_dyn_ratio_with_sigma_rot - kds_dyn_ratio_with_sigma_lower_error_rot,
      kds_dyn_ratio_with_sigma_rot + kds_dyn_ratio_with_sigma_upper_error_rot,
      kds_dyn_ratio_with_sigma_rot + kds_dyn_ratio_with_sigma_upper_error_rot],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

ax[2].errorbar(kds_redshift,
            kds_dyn_ratio_with_sigma_disp,
            ecolor='firebrick',
            marker='o',
            markersize=16,
            markerfacecolor='firebrick',
            markeredgecolor='firebrick',
            markeredgewidth=2,
            capsize=2,
            elinewidth=2,
            label=r'\textbf{9.8 - KDS (This study)}')

ax[2].fill([kds_redshift-0.05,
      kds_redshift+0.05,
      kds_redshift+0.05,
      kds_redshift-0.05],
     [kds_dyn_ratio_with_sigma_disp - kds_dyn_ratio_with_sigma_lower_error_disp,
      kds_dyn_ratio_with_sigma_disp - kds_dyn_ratio_with_sigma_lower_error_disp,
      kds_dyn_ratio_with_sigma_disp + kds_dyn_ratio_with_sigma_upper_error_disp,
      kds_dyn_ratio_with_sigma_disp + kds_dyn_ratio_with_sigma_upper_error_disp],
     'firebrick',
      edgecolor='firebrick',
      lw=2,
      alpha=0.25)

# split legend for the full samples
all_legend = []
all_legend.append([p1_all,p2_all,p3_all,p4_all,p5_all])
legend1 = ax[0].legend(all_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[0].add_artist(legend1)
ax[0].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for rotation dominated galaxies
rot_legend = []
rot_legend.append([p1_rot,p2_rot,p3_rot,p4_rot,p5_rot])
legend2 = ax[1].legend(rot_legend[0],
                    [r'\textbf{10.3 - DYNAMO (Green+14)}',
                     r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[1].add_artist(legend2)
ax[1].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

# split legend for dispersion dominated galaxies
disp_legend = []
disp_legend.append([p1_disp,p2_disp,p3_disp,p4_disp])
legend3 = ax[2].legend(disp_legend[0],
                    [r'\textbf{9.4 - MKS (Swinbank+17)}',
                     r'\textbf{9.9 - KROSS (Harrison+17)}',
                     r'\textbf{10.2 - MASSIV (Epinat+12)}',
                     r'\textbf{9.8 - MKS (Swinbank+17)}'],
                    loc=('lower left'),
                    prop={'size':11,'weight':'bold'},
                    frameon=False,
                    markerscale=0.75,
                    numpoints=1)
ax[2].add_artist(legend3)
ax[2].legend(loc='lower right',
          prop={'size':11,'weight':'bold'},
          frameon=False,
          markerscale=0.75,
          numpoints=1)

fig.tight_layout()
fig.subplots_adjust(wspace=0)
plt.show()
fig.savefig('/disk2/turner/disk1/turner/DATA/kmos_dynamics_paper_plots/PAPER_2_PLOTS/dyn_mass_with_sigma_evolution.png')
plt.close('all')
