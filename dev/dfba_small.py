from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable

import time

start_time = time.time()

# DfbaModel instance initialized with cobra model
small_model = read_sbml_model('/home/csiharath/Téléchargements/model_cancer_noDrugs.xml')
small_model.solver = "glpk"
# dfba_model = DfbaModel(small_model)


############################# Fonction de Xe ##############################
# Fonction biomass HMR
PYR = small_model.metabolites.get_by_id("PYR")
AMP = small_model.metabolites.get_by_id("AMP")
PALM = small_model.metabolites.get_by_id("PALM")


small_model.reactions.get_by_id('Vgrowth').add_metabolites({
    PYR: -10.0,
    AMP: -1.0,
    PALM: -1.0
})

small_model.objective = 'Vgrowth'
print(small_model.objective)

# Petit modèle humain
# Param pour G1 :
# PFK
small_model.reactions.get_by_id("VPFK").lower_bound = 0
small_model.reactions.get_by_id("VPFK").upper_bound = 1000
# G6PDH
small_model.reactions.get_by_id("VG6PDH").lower_bound = 0
small_model.reactions.get_by_id("VG6PDH").upper_bound = 1000
# TKT
small_model.reactions.get_by_id("VTK").lower_bound = 0
small_model.reactions.get_by_id("VTK").upper_bound = 1000
# PALM
small_model.reactions.get_by_id("VPALM").lower_bound = 0
small_model.reactions.get_by_id("VPALM").upper_bound = 1000

# solution = small_model.optimize()

# s = solution.fluxes
# dict_fluxes = s.to_frame().to_dict()['fluxes']
# dict_fluxes1 = {x:y for x,y in dict_fluxes.items() if y!=0}

############################# Dfba model definition #############################

dfba_model = DfbaModel(small_model)
dfba_model.add_objectives(["Vgrowth"], ["max"])
print(dfba_model.objectives)
#################################################################################


# instances of KineticVariable (default initial conditions are 0.0, but can be
# set here if wanted e.g. Oxygen)

X = KineticVariable("Biomass", initial_condition=1.04e-4)
Gluc = KineticVariable("EGLC")
Lac = KineticVariable("LAC")
Gln = KineticVariable("EGLN")
Amp = KineticVariable("AMP")
Pyr = KineticVariable("PYR")
Palm = KineticVariable("PALM")
AcoA = KineticVariable("ACCOA")
R5P = KineticVariable("R5P")

# add kinetic variables to dfba_model
dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln, Amp, Pyr, Palm, AcoA, R5P])
# dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln, Amp, Pyr, Palm])

# instances of ExchangeFlux
mu = ExchangeFlux("Vgrowth")
# v_G = ExchangeFlux("EX_EGLC") # glucose[x] <=> glucose[s]  HMR_9034
# v_N = ExchangeFlux("EX_EGLN") # EGLN[x] <=> EGLN[s] HMR_9063
v_G = ExchangeFlux("VHK")
v_N = ExchangeFlux("VGlnT")
v_L = ExchangeFlux("VLDH")
v_D = ExchangeFlux("VPPRibP")
v_P = ExchangeFlux("VPALM")


# add exchange fluxes to dfba_model
dfba_model.add_exchange_fluxes([mu, v_G, v_L, v_N, v_D, v_P])

# add rhs expressions for kinetic variables in dfba_model
dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("EGLC", -v_G * (X / 1000.0))
dfba_model.add_rhs_expression("EGLN", -v_N * (X / 1000.0))
dfba_model.add_rhs_expression("LAC", v_L * (X / 1000.0)) 
dfba_model.add_rhs_expression("AMP", v_D * (X / 1000.0))
dfba_model.add_rhs_expression("PALM", v_P * (X / 1000.0)) 
dfba_model.add_rhs_expression("PYR", v_G * (X / 1000.0))
dfba_model.add_rhs_expression("ACCOA", -v_P * (X / 1000.0))
dfba_model.add_rhs_expression("R5P", -v_D * (X / 1000.0))


# add lower/upper bound expressions for exchange fluxes in dfba_model together
# with expression that must be non-negative for correct evaluation of bounds

#dfba_model.add_exchange_flux_lb("Vgrowth", 2.25e-4 * (X/(4.29 + X)), X)
dfba_model.add_exchange_flux_lb("VHK", 2.25e-4 * (Gluc/(4.29 + Gluc)), Gluc)
dfba_model.add_exchange_flux_lb("VGlnT",8.7e-5 * (Gln/(1.27 + Gln)), Gln)
dfba_model.add_exchange_flux_lb("VLDH", 1.5e-3 * (Pyr/(4.29 + Pyr)), Pyr) #km lac
dfba_model.add_exchange_flux_lb("VPPRibP", 2e-8 * (R5P / (2 + R5P)), R5P)
dfba_model.add_exchange_flux_lb("VPALM", 2.95e-5 * (AcoA/(1 + AcoA)), AcoA) #

dfba_model.add_exchange_flux_ub("VHK", 10000, Gluc)
dfba_model.add_exchange_flux_ub("VGlnT",10000, Gln)
dfba_model.add_exchange_flux_ub("VLDH", 10000, Pyr) #km lac
dfba_model.add_exchange_flux_ub("VPPRibP", 10000, R5P)
dfba_model.add_exchange_flux_ub("VPALM", 10000, AcoA)

 
# add initial conditions for kinetic variables in dfba_model X (gDW/L),
# metabolites (g/L)
dfba_model.add_initial_conditions(
    {
        "Biomass": 1e-4,
        "EGLC": 1e-4,
        "EGLN": 1e-4,
        "LAC": 0,
        "PALM": 0
    }
)


# simulate model across interval t = [0.0,25.0](hours) with outputs for plotting
# every 0.1h and optional list of fluxes
concentrations, trajectories = dfba_model.simulate(
    0, 0.5, 0.01, ["Vgrowth", "VLDH", "VPPRibP", "VPALM"]
)

# print(concentrations)

# generate plots of results (in this case using plotlly)
import matplotlib.pyplot as plt
from dfba.plot.plotly import *

import plotly.io as pio

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrations)
# print(type(fig))
plt.ylim([0, 50])
fig.show()

fig = plot_trajectories(trajectories)
fig.show()

# write results to file

# concentrations.to_csv("concentrations.csv")
# trajectories.to_csv("trajectories.csv")


# # plt.rcParams["figure.figsize"] = 16, 9

# fig, ax = plt.subplots()
# concentrations.plot(ax=ax)

# plt.ylim([0, 50])
# plt.show() # plt in place of ax
# # concentrations.to_csv("concentrations.csv")

# from dfba.plot.matplotlib import *

# plt.rcParams["figure.figsize"] = 16, 9
# plt.rcParams["figure.autolayout"] = True

# plot_concentrations(concentrations)
# plt.ylim([0, 50])

print("--- %s seconds ---" % (time.time() - start_time))
# plt.show()

