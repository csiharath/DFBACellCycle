from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable

import matplotlib.pyplot as plt
from dfba.plot.plotly import *

import plotly.io as pio

import csv
from csv import writer
import time

import pandas as pd

############################# Variables Declaration ###########################

start_time = time.time()

######### Growth rates threshold ##########

xB = 0.1
alpha = 0.45
beta = 0.8
t1 = (1 + alpha) * xB
t2 = (1 + beta) * xB

############################## SBML Model for G1 ##############################

small_model_G1 = read_sbml_model('/home/csiharath/Stage/DFBACellCycle/small_model/model_mario2013_wocitG1.xml')
small_model_G1.solver = "glpk"

############################## Biomass Function ###############################

small_model_G1.objective = 'VBiomass'
print(small_model_G1.objective)


small_model_G1.reactions.get_by_id("VBiomass").lower_bound = 0
small_model_G1.reactions.get_by_id("VBiomass").upper_bound = 1


############################# dFBA Model Definition ###########################

dfba_model_G1 = DfbaModel(small_model_G1)
dfba_model_G1.add_objectives(["VBiomass"], ["max"])
print(dfba_model_G1.objectives)

############################## G1 Parameters ##################################

##### Instances of Kinetic Variables ######

X = KineticVariable("Biomass", initial_condition=1.04e-4)
Gluc = KineticVariable("EGLC")
Lac = KineticVariable("ELAC")
Gln = KineticVariable("EGLU")
#O2 = KineticVariable("O2")

### Add Kinetic variables to dfba_model ###

dfba_model_G1.add_kinetic_variables([X, Gluc, Lac, Gln])#, O2])

######## Instances of ExchangeFlux ########
mu = ExchangeFlux("VBiomass")
v_G = ExchangeFlux("EX_EGLC") 
v_N = ExchangeFlux("EX_EGLU")
v_L = ExchangeFlux("EX_ELAC")
#v_O = ExchangeFlux("EX_O2")


#### Add Exchange Fluxes to dfba_model ####

dfba_model_G1.add_exchange_fluxes([mu, v_G, v_L, v_N])#, v_O])

#### Expressions for Kinetic Variables ####
dfba_model_G1.add_rhs_expression("Biomass", mu * X)
dfba_model_G1.add_rhs_expression("EGLC", v_G * (X / 1000.0))
dfba_model_G1.add_rhs_expression("EGLU", v_N * (X / 1000.0))
dfba_model_G1.add_rhs_expression("ELAC", v_L * (X / 1000.0)) 
# dfba_model_G1.add_rhs_expression("O2", v_O * (X / 1000.0)) 


# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

dfba_model_G1.add_exchange_flux_lb("EX_EGLC", 10 * (Gluc/(3 + Gluc)), Gluc)
# dfba_model_G1.add_exchange_flux_lb("VBiomass", 2.25e-4 * (X/(4.29 + X)), X)
# dfba_model_G1.add_exchange_flux_lb("EX_EGLU",5 * (Gln/(0.0027 + Gln)), Gln)
# dfba_model_G1.add_exchange_flux_lb("EX_ELAC", 1 * (Lac/(0.1 + Lac)), Lac) 
# dfba_model_G1.add_exchange_flux_lb("EX_O2", 1, O2) 


# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######

dfba_model_G1.add_initial_conditions(
    {
        "Biomass": 0.1,
        "EGLC": 1.18,
        "EGLU": 0.1,
        "ELAC": 0#,
        #"O2": 0.16
    }
)

############ Simulate Model Across Interval t = [0.0,11.0](hours) #############
##### with outputs for plotting every 0.001h and optional list of fluxes ######

concentrationsG1, trajectories = dfba_model_G1.simulate(
    0, 11, 0.001, ["VBiomass","EX_EGLC", "EX_EGLU", "EX_ELAC"] #, "EX_O2"]
)


####### Cell Cycle Dhases Defined by Threshold ######

# out = True

# for index, row in concentrationsG1.iterrows():
#     if float(row[1]) >= t1:
#         if out:
#             S_init_list = [float(x) for x in row]
#             i = index
#             out = False

# concentrationsG1 = concentrationsG1[concentrationsG1.index <= i]
# name = concentrationsG1.columns.tolist()
# print(S_init_list)

# concentrationsG1.to_csv("concentrationsG1.csv")

# S_init = dict(zip(name, S_init_list))
# print(S_init)

#####################################################

concentrationsG1.to_csv("concentrationsG1.csv")

########## Gets latest concentrations values of Kinetic Variables ###############
####################### to use as initial condition next ########################

with open('concentrationsG1.csv', 'r') as csvfile:
    names = csvfile.readline()
    names = names.split('\n')[0].split(',')
    lastline = csvfile.readlines()[-1]
    lastline = lastline.split('\n')[0].split(',')
    val = int(lastline[0]) + 1
    S_init_list = [float(x) for x in lastline]

S_init = dict(zip(names, S_init_list))

############ generate plots of results ##############

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrationsG1)

plt.ylim([0, 50])
fig.show()

############################## SBML Model for S ###############################
small_model_S = read_sbml_model('/home/csiharath/Stage/DFBACellCycle/small_model/model_mario2013_wocitS.xml')
small_model_S.solver = "glpk"

############################## Biomass Function ###############################

small_model_S.objective = 'VBiomass'
print(small_model_S.objective)


small_model_S.reactions.get_by_id("VBiomass").lower_bound = 0
small_model_S.reactions.get_by_id("VBiomass").upper_bound = 1

############################# dFBA Model Definition ###########################

dfba_model_S = DfbaModel(small_model_S)
dfba_model_S.add_objectives(["VBiomass"], ["max"])
print(dfba_model_S.objectives)

############################### S Parameters ##################################

##### Instances of Kinetic Variables ######

X = KineticVariable("Biomass")
Gluc = KineticVariable("EGLC")
Lac = KineticVariable("ELAC")
Gln = KineticVariable("EGLU")
#O2 = KineticVariable("O2")

### Add Kinetic variables to dfba_model ###

dfba_model_S.add_kinetic_variables([X, Gluc, Lac, Gln])#, O2])

######## Instances of ExchangeFlux ########
mu = ExchangeFlux("VBiomass")
v_G = ExchangeFlux("EX_EGLC") 
v_N = ExchangeFlux("EX_EGLU")
v_L = ExchangeFlux("EX_ELAC")
#v_O = ExchangeFlux("EX_O2")

#### Add Exchange Fluxes to dfba_model ####

dfba_model_S.add_exchange_fluxes([mu, v_G, v_L, v_N])#, v_O])

#### Expressions for Kinetic Variables ####
dfba_model_S.add_rhs_expression("Biomass", mu * X)
dfba_model_S.add_rhs_expression("EGLC", v_G * (X / 1000.0))
dfba_model_S.add_rhs_expression("EGLU", v_N * (X / 1000.0))
dfba_model_S.add_rhs_expression("ELAC", v_L * (X / 1000.0)) 
#dfba_model_S.add_rhs_expression("O2", v_O * (X / 1000.0))

# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

dfba_model_S.add_exchange_flux_lb("EX_EGLC", 10 * (Gluc/(5 + Gluc)), Gluc)
#dfba_model_S.add_exchange_flux_lb("EX_O2", 1, O2) 

# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######
dfba_model_S.add_initial_conditions(
    {
        "Biomass": S_init["Biomass"],
        "EGLC": S_init["EGLC"],
        "EGLU": S_init["EGLU"],
        "ELAC": S_init["ELAC"]#,
        #"O2": 0.16
    }
)

############ Simulate Model Across Interval t = [11.0,19.0](hours) ############
##### with outputs for plotting every 0.001h and optional list of fluxes ######

concentrationsS, trajectories = dfba_model_S.simulate(
    11, 19, 0.001, ["VBiomass","EX_EGLC", "EX_EGLU", "EX_ELAC", "EX_O2"] # "VPPRIBP", "VPALM", "VHK", "VAKGDH"]
)

####### Cell Cycle Dhases Defined by Threshold ######

# out = True

# for index, row in concentrationsS.iterrows():
#     if float(row[1]) >= t2:
#         if out:
#             G2_init_list = [float(x) for x in row]
#             i = index
#             out = False

# concentrationsS = concentrationsS[concentrationsS.index <= i]
# name = concentrationsS.columns.tolist()
# print(G2_init_list)

# concentrationsS.to_csv("concentrationsS.csv")

# G2_init = dict(zip(name, G2_init_list))
# print(G2_init)

#####################################################

concentrationsS.to_csv("concentrationsS.csv")

########## Gets latest concentrations values of Kinetic Variables ###############
####################### to use as initial condition next ########################
with open('concentrationsS.csv', 'r') as csvfile:
    names = csvfile.readline()
    names = names.split('\n')[0].split(',')
    lastline = csvfile.readlines()[-1]
    lastline = lastline.split('\n')[0].split(',')
    val = val + int(lastline[0]) + 1
    G2_init_list = [float(x) for x in lastline]

G2_init = dict(zip(names, G2_init_list))


######## Concatenates concentrations values #########
############## with previous phase(s) ###############

with open('concentrationsS.csv') as csvfile_S:
    csvreader = csv.reader(csvfile_S, delimiter=',', quotechar='|')
    with open('concentrationsG1.csv', 'a', newline='') as csvfile_G1:
        for row in csvreader:
            if row[0]!= '':
                row[0] = int(row[0]) + val
                writer_object = writer(csvfile_G1)
                writer_object.writerow(row)

########### generate plots of results ###############

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrationsS)
# print(type(fig))
plt.ylim([0, 50])
fig.show()

############################## SBML Model for G2 ##############################

small_model_G2 = read_sbml_model('/home/csiharath/Stage/DFBACellCycle/small_model/model_mario2013_wocitG2.xml')
small_model_G2.solver = "glpk"

############################## Biomass Function ###############################

small_model_G2.objective = 'VBiomass'
print(small_model_G2.objective)


small_model_G2.reactions.get_by_id("VBiomass").lower_bound = 0
small_model_G2.reactions.get_by_id("VBiomass").upper_bound = 1

############################# dFBA Model Definition ###########################

dfba_model_G2 = DfbaModel(small_model_G2)
dfba_model_G2.add_objectives(["VBiomass"], ["max"])
print(dfba_model_G2.objectives)

############################## G2 Parameters ##################################

##### Instances of Kinetic Variables ######

X = KineticVariable("Biomass")
Gluc = KineticVariable("EGLC")
Lac = KineticVariable("ELAC")
Gln = KineticVariable("EGLU")
#O2 = KineticVariable("O2")

### Add Kinetic variables to dfba_model ###

dfba_model_G2.add_kinetic_variables([X, Gluc, Lac, Gln])#, #O2])

######## Instances of ExchangeFlux ########
mu = ExchangeFlux("VBiomass")
v_G = ExchangeFlux("EX_EGLC") 
v_N = ExchangeFlux("EX_EGLU")
v_L = ExchangeFlux("EX_ELAC")
# v_O = ExchangeFlux("EX_O2")


#### Add Exchange Fluxes to dfba_model ####

dfba_model_G2.add_exchange_fluxes([mu, v_G, v_L, v_N])# , v_O])

#### Expressions for Kinetic Variables ####
dfba_model_G2.add_rhs_expression("Biomass", mu * X)
dfba_model_G2.add_rhs_expression("EGLC", v_G * (X / 1000.0))
dfba_model_G2.add_rhs_expression("EGLU", v_N * (X / 1000.0))
dfba_model_G2.add_rhs_expression("ELAC", v_L * (X / 1000.0)) 
#dfba_model_G2.add_rhs_expression("O2", v_O * (X / 1000.0))

# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

dfba_model_G2.add_exchange_flux_lb("EX_EGLC", 10 * (Gluc/(5 + Gluc)), Gluc)
# dfba_model_G2.add_exchange_flux_lb("EX_O2", 1, O2) 

# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######

dfba_model_G2.add_initial_conditions(
    {
        "Biomass": G2_init["Biomass"],
        "EGLC": G2_init["EGLC"],
        "EGLU": G2_init["EGLU"],
        "ELAC": G2_init["ELAC"]#,
        #"O2": 0.16
    }
)

############ Simulate Model Across Interval t = [0.0,11.0](hours) #############
##### with outputs for plotting every 0.001h and optional list of fluxes ######

concentrationsG2, trajectories = dfba_model_G2.simulate(
    19, 23, 0.001, ["VBiomass","EX_EGLC", "EX_EGLU", "EX_ELAC", "EX_O2"]
)

concentrationsG2.to_csv("concentrationsG2.csv")

######## Concatenates concentrations values #########
############## with previous phase(s) ###############

with open('concentrationsG2.csv') as csvfile_G2:
    csvreader = csv.reader(csvfile_G2, delimiter=',', quotechar='|')
    with open('concentrationsG1.csv', 'a', newline='') as csvfile_G1:
        for row in csvreader:
            if row[0]!= '':
                row[0] = int(row[0]) + val
                writer_object = writer(csvfile_G1)
                writer_object.writerow(row)

########### generate plots of results ###############

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrationsG2)
# print(type(fig))
plt.ylim([0, 50])
fig.show()


###### generate plots of results for a cycle ########

concentrations = pd.read_csv("concentrationsG1.csv", index_col=[0])

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrations)
# print(type(fig))
plt.ylim([0, 50])
fig.show()

print(type(concentrationsG1))

######
# trajectories.to_csv("trajectories.csv")

# fig = plot_trajectories(trajectories)
# fig.show()

print("--- %s seconds ---" % (time.time() - start_time))
