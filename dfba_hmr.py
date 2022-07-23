from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable

import matplotlib.pyplot as plt
from dfba.plot.plotly import *

import plotly.io as pio

import time

############################# Variables Declaration ###########################

start_time = time.time()

############################## SBML Model for G1 ##############################

# model_recon = read_sbml_model('/home/csiharath/Stage/DFBACellCycle/HMRdatabase2_00_S.xml')
model_recon_G1 = read_sbml_model('/home/csiharath/Téléchargements/Recon3D_G1.xml')

model_recon_G1.solver = "glpk"
# dfba_model = DfbaModel(model_recon_G1)

############################## Biomass Function ###############################

model_recon_G1.objective = 'BIOMASS_reaction'
print(model_recon_G1.objective)


model_recon_G1.reactions.get_by_id("BIOMASS_reaction").lower_bound = 0
model_recon_G1.reactions.get_by_id("BIOMASS_reaction").upper_bound = 1

model_recon_G1.reactions.get_by_id("HEX1").lower_bound = 0.1
model_recon_G1.reactions.get_by_id("LDH_L").upper_bound = -0.1


# model_recon_G1.reactions.get_by_id("EX_atp_e").lower_bound = 1
# model_recon_G1.reactions.get_by_id("EX_atp_e").upper_bound = 0

# model_recon_G1.reactions.get_by_id("EX_lac__L_e").lower_bound = 1
model_recon_G1.reactions.get_by_id("GLCt2_2").lower_bound = 0.1
model_recon_G1.reactions.get_by_id("EX_glc__D_e").upper_bound = -0.1
# model_recon_G1.reactions.get_by_id("EX_lac__L_e").lower_bound = 0.1

# solution = model_recon_G1.optimize()
# # print(solution.objective_value)
# # print(solution.status)
# # print(solution.shadow_prices)
# # model_recon_G1.summary()
# s = solution.fluxes
# dict_fluxes = s.to_frame().to_dict()['fluxes']
# dict_fluxes1 = {x:y for x,y in dict_fluxes.items() if y!=0}

# for react in dict_fluxes1:
#     if 'EX_' in react:
#         print(f"{react} : {dict_fluxes1[react]}")

# print(dict_fluxes["PYK"])

# for f in dict_fluxes1:
#     # print(model_recon_G1.reactions.get_by_id(f).reactants[0])
#     for reactant in model_recon_G1.reactions.get_by_id(f).reactants:
#         # print(reactant.id + 'm01845c')
#         if "pyr" in reactant.id:
#             print(model_recon_G1.reactions.get_by_id(f))
#             print(dict_fluxes1[f])
#     for product in model_recon_G1.reactions.get_by_id(f).products:
#             if "pyr" in product.id:
#                 print("\nreversible:")
#                 print(model_recon_G1.reactions.get_by_id(f))
#                 print(dict_fluxes1[f])

############################# dFBA Model Definition ###########################

dfba_model = DfbaModel(model_recon_G1)
dfba_model.add_objectives(["BIOMASS_reaction"], ["max"])
print(dfba_model.objectives)

############################## G1 Parameters ##################################


##### Instances of Kinetic Variables ######

# X = KineticVariable("Biomass", initial_condition=1.04e-4)
# Gluc = KineticVariable("EGLC")
# Lac = KineticVariable("ELAC")
# Gln = KineticVariable("EGLN")
# #O2 = KineticVariable("O2")

X = KineticVariable("Biomass")
Gluc = KineticVariable("glc__D_e")
Lac = KineticVariable("lac__L_e")
Gln = KineticVariable("gln__L_e")
# O2 = KineticVariable("o2_e")

### Add Kinetic variables to dfba_model ###

dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln])#, O2])

######## Instances of ExchangeFlux ########
# mu = ExchangeFlux("biomass_components")
# v_G = ExchangeFlux("EX_EGLC") 
# v_N = ExchangeFlux("EX_EGLN")
# v_L = ExchangeFlux("EX_ELAC")
# #v_O = ExchangeFlux("EX_O2")

mu = ExchangeFlux("BIOMASS_reaction")
v_G = ExchangeFlux("EX_glc__D_e") 
v_N = ExchangeFlux("EX_gln__L_e")
v_L = ExchangeFlux("EX_lac__L_e")
# v_O = ExchangeFlux("EX_o2_e")


#### Add Exchange Fluxes to dfba_model ####

dfba_model.add_exchange_fluxes([mu, v_G, v_L, v_N])#, v_O])

#### Expressions for Kinetic Variables ####

# dfba_model.add_rhs_expression("Biomass", mu * X)
# dfba_model.add_rhs_expression("EGLC", v_G * (X / 1000.0))
# dfba_model.add_rhs_expression("EGLN", v_N * (X / 1000.0))
# dfba_model.add_rhs_expression("ELAC", v_L * (X / 1000.0)) 
# #dfba_model.add_rhs_expression("O2", v_O * (X / 1000.0)) 

dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("glc__D_e", v_G * (X / 1000.0))
dfba_model.add_rhs_expression("gln__L_e", v_N * (X / 1000.0))
dfba_model.add_rhs_expression("lac__L_e", v_L * (X / 1000.0)) 
# dfba_model.add_rhs_expression("o2_e", v_O * (X / 1000.0))


# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

# #dfba_model.add_exchange_flux_lb("VBiomass", 2.25e-4 * (X/(4.29 + X)), X)
# dfba_model.add_exchange_flux_lb("EX_EGLC", 10 * (Gluc/(3 + Gluc)), Gluc)
# dfba_model.add_exchange_flux_lb("EX_EGLN",5 * (Gln/(0.0027 + Gln)), Gln)
# # dfba_model.add_exchange_flux_lb("EX_ELAC", 1 * (Lac/(0.1 + Lac)), Lac) 
# #dfba_model.add_exchange_flux_lb("EX_O2", 1, O2) 

dfba_model.add_exchange_flux_lb("EX_glc__D_e", 10 * (Gluc/(3 + Gluc)), Gluc)
# dfba_model.add_exchange_flux_ub("EX_glc__D_e", 0.5 )
dfba_model.add_exchange_flux_lb("EX_gln__L_e", 5 * (Gln/(0.0027 + Gln)), Gln)
# dfba_model.add_exchange_flux_lb("EX_o2_e", 1, O2) 

# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######

dfba_model.add_initial_conditions(
    {
        "Biomass": 0.1,
        "glc__D_e": 1.18,
        "gln__L_e": 0.1,
        "lac__L_e": 0#,
        #"o2_e": 0.16
    }
)


############ Simulate Model Across Interval t = [0.0,11.0](hours) #############
##### with outputs for plotting every 0.001h and optional list of fluxes ######

# concentrationsG1, trajectories = dfba_model.simulate(
#     0, 11, 0.001, ["biomass_components","EX_EGLC", "EX_EGLN", "EX_ELAC"] # "VPPRIBP", "VPALM", "VHK", "VAKGDH"]
# )
concentrationsG1, trajectories = dfba_model.simulate(
    0, 11, 0.001, ["BIOMASS_reaction","EX_glc__D_e", "EX_gln__L_e", "EX_lac__L_e"]#, "EX_o2_e"] # "VPPRIBP", "VPALM", "VHK", "VAKGDH"]
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
print(S_init)


############ generate plots of results ##############

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrationsG1)
# print(type(fig))
plt.ylim([0, 50])
fig.show()

fig = plot_trajectories(trajectories)
fig.show()

# S_init = {'': 11000.0, 'time': 10.999999999999343, 'Biomass': 6390.679700629999, 'glc__D_e': 0.5409420299370008, 'gln__L_e': 0.1, 'lac__L_e': 24.16111829250855}
# print(S_init)

"""############################## SBML Model for S ###############################

model_recon = read_sbml_model('/home/csiharath/Téléchargements/Recon3D_S.xml')

model_recon.solver = "glpk"

############################## Biomass Function ###############################

model_recon.objective = 'BIOMASS_reaction'
print(model_recon.objective)


model_recon.reactions.get_by_id("BIOMASS_reaction").lower_bound = 0
model_recon.reactions.get_by_id("BIOMASS_reaction").upper_bound = 1


model_recon.reactions.get_by_id("G6PDH1er").lower_bound = 0.1
# model_recon.reactions.get_by_id("LDH_L").upper_bound = -0.1

# model_recon.reactions.get_by_id("EX_lac__L_e").lower_bound = 1
# model_recon.reactions.get_by_id("GLCt2_2").lower_bound = 0.1
# model_recon.reactions.get_by_id("EX_glc__D_e").upper_bound = -0.1
# model_recon.reactions.get_by_id("EX_lac__L_e").lower_bound = 0.1


solution = model_recon.optimize()
# print(solution.objective_value)
# print(solution.status)
# print(solution.shadow_prices)
# model_recon_G1.summary()
s = solution.fluxes
dict_fluxes = s.to_frame().to_dict()['fluxes']
dict_fluxes1 = {x:y for x,y in dict_fluxes.items() if y!=0}

for react in dict_fluxes1:
    if 'EX_' in react:
        print(f"{react} : {dict_fluxes1[react]}")

print(dict_fluxes["PRPPS"])

for f in dict_fluxes1:
    # print(model_recon_G1.reactions.get_by_id(f).reactants[0])
    for reactant in model_recon.reactions.get_by_id(f).reactants:
        # print(reactant.id + 'm01845c')
        if "pyr" in reactant.id:
            print(model_recon.reactions.get_by_id(f))
            print(dict_fluxes1[f])
    for product in model_recon.reactions.get_by_id(f).products:
            if "pyr" in product.id:
                print("\nreversible:")
                print(model_recon.reactions.get_by_id(f))
                print(dict_fluxes1[f])

############################# dFBA Model Definition ###########################

dfba_model = DfbaModel(model_recon)
dfba_model.add_objectives(["BIOMASS_reaction"], ["max"])
print(dfba_model.objectives)

############################### S Parameters ##################################

##### Instances of Kinetic Variables ######

X = KineticVariable("Biomass")
Gluc = KineticVariable("glc__D_e")
Lac = KineticVariable("lac__L_e")
Gln = KineticVariable("gln__L_e")
# O2 = KineticVariable("o2_e")

### Add Kinetic variables to dfba_model ###

dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln])#, O2])

######## Instances of ExchangeFlux ########

mu = ExchangeFlux("BIOMASS_reaction")
v_G = ExchangeFlux("EX_glc__D_e") 
v_N = ExchangeFlux("EX_gln__L_e")
v_L = ExchangeFlux("EX_lac__L_e")
# v_O = ExchangeFlux("EX_o2_e")


#### Add Exchange Fluxes to dfba_model ####

dfba_model.add_exchange_fluxes([mu, v_G, v_L, v_N])#, v_O])

#### Expressions for Kinetic Variables ####
dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("glc__D_e", v_G * (X / 1000.0))
dfba_model.add_rhs_expression("gln__L_e", v_N * (X / 1000.0))
dfba_model.add_rhs_expression("lac__L_e", v_L * (X / 1000.0)) 
# dfba_model.add_rhs_expression("o2_e", v_O * (X / 1000.0))


# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

dfba_model.add_exchange_flux_lb("EX_glc__D_e", 10 * (Gluc/(3 + Gluc)), Gluc)
# dfba_model.add_exchange_flux_ub("EX_glc__D_e", 0.5 )
# dfba_model.add_exchange_flux_lb("EX_gln__L_e", 5 * (Gln/(0.0027 + Gln)), Gln)
# dfba_model.add_exchange_flux_lb("EX_o2_e", 1, O2) 

# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######
dfba_model.add_initial_conditions(
    {
        "Biomass": S_init["Biomass"],
        "glc__D_e": S_init["glc__D_e"],
        "gln__L_e": S_init["gln__L_e"],
        "lac__L_e": S_init["lac__L_e"]#,
        #"o2_e": 0.16
    }
)


############ Simulate Model Across Interval t = [11.0,19.0](hours) ############
##### with outputs for plotting every 0.001h and optional list of fluxes ######

concentrationsS, trajectories = dfba_model.simulate(
    11, 19, 0.001, ["BIOMASS_reaction","EX_glc__D_e", "EX_gln__L_e", "EX_lac__L_e"]#, "EX_o2_e"] # "VPPRIBP", "VPALM", "VHK", "VAKGDH"]
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
    val = int(lastline[0]) + 1
    G2_init_list = [float(x) for x in lastline]

G2_init = dict(zip(names, G2_init_list))
print(G2_init)


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

model_recon = read_sbml_model('/home/csiharath/Téléchargements/Recon3D_G2.xml')

model_recon.solver = "glpk"

############################## Biomass Function ###############################

model_recon.objective = 'BIOMASS_reaction'
print(model_recon.objective)


model_recon.reactions.get_by_id("BIOMASS_reaction").lower_bound = 0
model_recon.reactions.get_by_id("BIOMASS_reaction").upper_bound = 1

# model_recon.reactions.get_by_id("EX_glcn_e").lower_bound = 0
# model_recon.reactions.get_by_id("EX_glcn_e").upper_bound = 0

# model_recon.reactions.get_by_id("r0354").lower_bound = 10

# model_recon.reactions.get_by_id("EX_atp_e").lower_bound = 10
# model_recon.reactions.get_by_id("EX_atp_e").upper_bound = 0


# model_recon.reactions.get_by_id("GLCt2_2").lower_bound = 10
# model_recon.reactions.get_by_id("EX_glc__D_e").upper_bound = -0.1


# solution = model_recon.optimize()

# s = solution.fluxes
# dict_fluxes = s.to_frame().to_dict()['fluxes']
# dict_fluxes1 = {x:y for x,y in dict_fluxes.items() if y!=0}

# for react in dict_fluxes1:
#     if 'EX_' in react or react == "r0354":
#         print(f"{react} : {dict_fluxes1[react]}")

# for f in dict_fluxes1:
#     # print(model_recon.reactions.get_by_id(f).reactants[0])
#     for reactant in model_recon.reactions.get_by_id(f).reactants:
#         # print(reactant.id + 'm01845c')
#         if "glc__D" in reactant.id:
#             print(model_recon.reactions.get_by_id(f))
#             print(dict_fluxes1[f])
#     for product in model_recon.reactions.get_by_id(f).products:
#             if "glc__D" in product.id:
#                 print("\nreversible:")
#                 print(model_recon.reactions.get_by_id(f))
#                 print(dict_fluxes1[f])
#     # if model_recon.reactions.get_by_id(f).lower_bound < 0:
#     #     for product in model_recon.reactions.get_by_id(f).products:
#     #         if "glc__D" in product.id:
#     #             print("\nreversible:")
#     #             print(model_recon.reactions.get_by_id(f))
#     #             print(dict_fluxes1[f])

############################# dFBA Model Definition ###########################

dfba_model = DfbaModel(model_recon)
dfba_model.add_objectives(["BIOMASS_reaction"], ["max"])
print(dfba_model.objectives)

############################## G2 Parameters ##################################

##### Instances of Kinetic Variables ######

X = KineticVariable("Biomass")
Gluc = KineticVariable("glc__D_e")
Lac = KineticVariable("lac__L_e")
Gln = KineticVariable("gln__L_e")
# O2 = KineticVariable("o2_e")

### Add Kinetic variables to dfba_model ###

dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln])#, O2])

######## Instances of ExchangeFlux ########

mu = ExchangeFlux("BIOMASS_reaction")
v_G = ExchangeFlux("EX_glc__D_e") 
v_N = ExchangeFlux("EX_gln__L_e")
v_L = ExchangeFlux("EX_lac__L_e")
# v_O = ExchangeFlux("EX_o2_e")


#### Add Exchange Fluxes to dfba_model ####

dfba_model.add_exchange_fluxes([mu, v_G, v_L, v_N])#, v_O])

#### Expressions for Kinetic Variables ####
dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("glc__D_e", v_G * (X / 1000.0))
dfba_model.add_rhs_expression("gln__L_e", v_N * (X / 1000.0))
dfba_model.add_rhs_expression("lac__L_e", v_L * (X / 1000.0)) 
# dfba_model.add_rhs_expression("o2_e", v_O * (X / 1000.0))


# Add lower/upper Bound Expressions for exchange fluxes 
# Expressions must be non-negative for correct evaluation of bounds

# dfba_model.add_exchange_flux_lb("EX_glc__D_e", 10 * (Gluc/(3 + Gluc)), Gluc)
# dfba_model.add_exchange_flux_ub("EX_glc__D_e", 0.5 )
dfba_model.add_exchange_flux_lb("EX_gln__L_e", 5 * (Gln/(0.0027 + Gln)), Gln)
# dfba_model.add_exchange_flux_lb("EX_o2_e", 1, O2) 

# Initial Conditions of Kinetic Variables #
####### X (gDW/L), metabolites (g/L) ######

dfba_model.add_initial_conditions(
    {
        "Biomass": G2_init["Biomass"],
        "glc__D_e": G2_init["glc__D_e"],
        "gln__L_e": G2_init["gln__L_e"],
        "lac__L_e": G2_init["lac__L_e"]#,
        #"o2_e": 0.16
    }
)


# simulate model across interval t = [0.0,25.0](hours) with outputs for plotting
# every 0.1h and optional list of fluxes
concentrationsG2, trajectories = dfba_model.simulate(
    19, 23, 0.001, ["BIOMASS_reaction","EX_glc__D_e", "EX_gln__L_e", "EX_lac__L_e"]#, "EX_o2_e"] # "VPPRIBP", "VPALM", "VHK", "VAKGDH"]
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
# fig.show()"""

print("--- %s seconds ---" % (time.time() - start_time))