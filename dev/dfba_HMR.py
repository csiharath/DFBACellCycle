from os.path import dirname, join, pardir

from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable

import time

start_time = time.time()

# DfbaModel instance initialized with cobra model
fba_model_HMR = read_sbml_model("/home/csiharath/Stage/model/HMRdatabase2_00.xml")
fba_model_HMR.solver = "glpk"
# dfba_model = DfbaModel(fba_model_HMR)

############################# Metabolites utilisees #############################

fba_model_HMR.metabolites.get_by_id('m01968c').id = "G6Pc"
fba_model_HMR.metabolites.get_by_id('m01968r').id = "G6Pr"
fba_model_HMR.metabolites.get_by_id('m01845c').id = "F6P"
fba_model_HMR.metabolites.get_by_id('m01939c').id = "GAP"
fba_model_HMR.metabolites.get_by_id('m02696c').id = "PEP"
fba_model_HMR.metabolites.get_by_id('m02819c').id = "PYR"
fba_model_HMR.metabolites.get_by_id('m02845c').id = "R5P"
fba_model_HMR.metabolites.get_by_id('m01761c').id = "X5P"
fba_model_HMR.metabolites.get_by_id('m01261c').id = "ACoA"
fba_model_HMR.metabolites.get_by_id('m01587c').id = "CIT"
fba_model_HMR.metabolites.get_by_id('m01306c').id = "AKG"
fba_model_HMR.metabolites.get_by_id("m02943c").id = "SUC"
fba_model_HMR.metabolites.get_by_id("m02439c").id = "MAL"
fba_model_HMR.metabolites.get_by_id("m02633c").id = "OXA"
fba_model_HMR.metabolites.get_by_id("m01716c").id = "LAC"
fba_model_HMR.metabolites.get_by_id("m02674c").id = "PALM"
fba_model_HMR.metabolites.get_by_id("m01974c").id = "GLU"
fba_model_HMR.metabolites.get_by_id("m01638c").id = "ALA"
fba_model_HMR.metabolites.get_by_id("m02579c").id = "NH4"
fba_model_HMR.metabolites.get_by_id("m01965x").id = "GLCex"
fba_model_HMR.metabolites.get_by_id("m01974x").id = "GLUex"
fba_model_HMR.metabolites.get_by_id("m01975x").id = "GLNex"
fba_model_HMR.metabolites.get_by_id("m01371c").id = "ATP"
fba_model_HMR.metabolites.get_by_id("m01285c").id = "ADP"
fba_model_HMR.metabolites.get_by_id("m01334c").id = "AMP"
fba_model_HMR.metabolites.get_by_id("m02552c").id = "NAD"
fba_model_HMR.metabolites.get_by_id("m02553c").id = "NADH"
fba_model_HMR.metabolites.get_by_id("m02554c").id = "NADP"
fba_model_HMR.metabolites.get_by_id("m02555c").id = "NADHP"
fba_model_HMR.metabolites.get_by_id("m01721c").id = "DNA"

fba_model_HMR.metabolites.get_by_id("m02040p").id = "H2O"
fba_model_HMR.metabolites.get_by_id("m02678p").id = "Palm_coa"
fba_model_HMR.metabolites.get_by_id("m01597p").id = "COA"
fba_model_HMR.metabolites.get_by_id("m02674p").id = "Palm_ac"

############################## Reactions utilisees ##############################

fba_model_HMR.reactions.get_by_id("HMR_4394").id = "HK"
fba_model_HMR.reactions.get_by_id("HMR_4381").id = "PGI_rev"
fba_model_HMR.reactions.get_by_id("HMR_4379").id = "PFK1"
fba_model_HMR.reactions.get_by_id("HMR_4375").id = "PFK2_rev"
fba_model_HMR.reactions.get_by_id("HMR_4391").id = "PFK3"
fba_model_HMR.reactions.get_by_id("HMR_4373").id = "PGK1_rev"
fba_model_HMR.reactions.get_by_id("HMR_4368").id = "PGK2"
fba_model_HMR.reactions.get_by_id("HMR_4365").id = "PGK3"
fba_model_HMR.reactions.get_by_id("HMR_4363").id = "PGK4"
fba_model_HMR.reactions.get_by_id("HMR_4358").id = "PK"

fba_model_HMR.reactions.get_by_id("HMR_4306").id = "G6PDH1"
fba_model_HMR.reactions.get_by_id("HMR_4625").id = "G6PDH2"
fba_model_HMR.reactions.get_by_id("HMR_4474").id = "G6PDH3"
fba_model_HMR.reactions.get_by_id("HMR_4352").id = "G6PDH4_rev_EP1"
fba_model_HMR.reactions.get_by_id("HMR_4052").id = "AMP"
fba_model_HMR.reactions.get_by_id("HMR_7160").id = "DNA"
#fba_model_HMR.reactions.get_by_id("HMR_4352").id = "EP1"
fba_model_HMR.reactions.get_by_id("HMR_4477").id = "EP2_rev"
fba_model_HMR.reactions.get_by_id("HMR_4501").id = "TKT1"
fba_model_HMR.reactions.get_by_id("HMR_4565").id = "TKT2"


fba_model_HMR.reactions.get_by_id("HMR_4145").id = "CS"
fba_model_HMR.reactions.get_by_id("HMR_4456").id = "CITS1"
fba_model_HMR.reactions.get_by_id("HMR_4586").id = "CITS2"
fba_model_HMR.reactions.get_by_id("HMR_5297").id = "AKGDH1"
fba_model_HMR.reactions.get_by_id("HMR_4152").id = "AKGDH2"

fba_model_HMR.reactions.get_by_id("HMR_4139").id = "MLD"
fba_model_HMR.reactions.get_by_id("HMR_4143").id = "PC"

fba_model_HMR.reactions.get_by_id("HMR_4388").id = "LDH" # HMR_4281(p) 

fba_model_HMR.reactions.get_by_id("HMR_4149").id = "ACL"
fba_model_HMR.reactions.get_by_id("HMR_0709").id = "Vpalm"


############################# Fonction de Biomasse ##############################
########################### Coefficient pour phase G1 ###########################

G6P = fba_model_HMR.metabolites.get_by_id("G6Pc")
PYR = fba_model_HMR.metabolites.get_by_id("PYR")
ACCOA = fba_model_HMR.metabolites.get_by_id('ACoA')
R5P = fba_model_HMR.metabolites.get_by_id("R5P")
AMP = fba_model_HMR.metabolites.get_by_id("AMP")
PALM = fba_model_HMR.metabolites.get_by_id("PALM")
ATP = fba_model_HMR.metabolites.get_by_id("ATP")
X = fba_model_HMR.metabolites.get_by_id("temp001x")
ADP = fba_model_HMR.metabolites.get_by_id("ADP")
Pi = fba_model_HMR.metabolites.get_by_id("m02751c")
DNA = fba_model_HMR.metabolites.get_by_id("DNA")

fba_model_HMR.reactions.get_by_id('biomass_components').add_metabolites({
    PYR: -10.0,
    DNA: -1.0,
    PALM: -1.0,
    ATP: -1.0,
    X: 1.0,
    ADP: 1.0,
    Pi: 1.0
})

fba_model_HMR.objective = 'biomass_components'
# fba_model_HMR.reactions.biomass_components

############################# Dfba model definition #############################

dfba_model = DfbaModel(fba_model_HMR)

#################################################################################


# instances of KineticVariable (default initial conditions are 0.0, but can be
# set here if wanted e.g. Oxygen)

X = KineticVariable("Biomass", initial_condition=1.04e-4)
Gluc = KineticVariable("Glucose")
Lac = KineticVariable("Lactate")
Gln = KineticVariable("Glutamine")
Dna = KineticVariable("DNA")

# G6p = KineticVariable("Glucose 6-phosphate")
# F6p = KineticVariable("Fructose 6-phosphate")
# Gap = KineticVariable("Glyceraldehyde 3-phosphate")
# Pep = KineticVariable("Phosphoenolpyruvate")

Pyr = KineticVariable("Pyruvate")

# R5p = KineticVariable("ribose 5-phosphate")
# X5p = KineticVariable("xylulose 5-phosphate")

AcoA = KineticVariable("Acetyl-CoA")

# Cit = KineticVariable("Citrate")
# Akg = KineticVariable("2-Oxoglutarate")
# Suc = KineticVariable("Succinate")
# Mal = KineticVariable("Malate")
# Oxa = KineticVariable("Oxaloacetate")
Palm = KineticVariable("Palmitate")
# Glu = KineticVariable("Glutamate")
# Ala = KineticVariable("Alanine")
# Nh4 = KineticVariable("Ammonium")
# Glu_ex = KineticVariable("Glutamate ext")

Atp = KineticVariable("ATP", initial_condition=2.63)
Adp = KineticVariable("ADP", initial_condition=0.28)
Amp = KineticVariable("AMP", initial_condition=0.197)
Nad = KineticVariable("NAD", initial_condition=0.44)
Nadh = KineticVariable("NADH", initial_condition=5.41)
Nadp = KineticVariable("NADP", initial_condition=0.222)
Nadph = KineticVariable("NADPH", initial_condition=0.146)


# add kinetic variables to dfba_model
# dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln, Dna, Palm, Pyr, Atp, Adp, Amp, Nad, Nadh])
dfba_model.add_kinetic_variables([X, Gluc, Lac, Gln, Dna, Pyr, AcoA, Palm, Atp, Adp, Amp, Nad, Nadh, Nadp, Nadph])
# fba_model_HMR.add_kinetic_variables([G6p, F6p, Gap, Pep, Pyr, Lac, R5p, X5p, AcoA, Cit, Akg, Suc, Mal, Oxa, Palm, Glu, Ala, Nh4, Glu_ex, Atp, Adp, Amp, Nad, Nadh, Nadp, Nadph, X])


# instances of ExchangeFlux
mu = ExchangeFlux("biomass_components")
v_G = ExchangeFlux("EX_m01965x") # glucose[x] <=> glucose[s]  HMR_9034
v_N = ExchangeFlux("EX_m01975x") # glutamine[x] <=> glutamine[s] HMR_9063
v_L = ExchangeFlux("LDH")
v_D = ExchangeFlux("DNA")
v_P = ExchangeFlux("Vpalm")


# add exchange fluxes to dfba_model
dfba_model.add_exchange_fluxes([mu, v_G, v_L, v_N, v_D, v_P])

# add rhs expressions for kinetic variables in dfba_model
dfba_model.add_rhs_expression("Biomass", mu * X)
dfba_model.add_rhs_expression("Glucose", -v_G * 16.0 * X / 1000.0)
dfba_model.add_rhs_expression("Glutamine", -v_N * 16.0 * X / 1000.0)
dfba_model.add_rhs_expression("Lactate", v_L * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("DNA", v_D * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("Palmitate", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("Pyruvate", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("Acetyl-CoA", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("ATP", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("ADP", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("AMP", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("NAD", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("NADH", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("NADP", v_P * 150.13 * X / 1000.0)
dfba_model.add_rhs_expression("NADPH", v_P * 150.13 * X / 1000.0)


# add lower/upper bound expressions for exchange fluxes in dfba_model together
# with expression that must be non-negative for correct evaluation of bounds
# dfba_model.add_exchange_flux_lb(
#     "GLCex", 10.5 * (Gluc / (0.0027 + Gluc)) * (1 / (1 + Eth / 20.0)), Gluc
# )
dfba_model.add_exchange_flux_lb(
    "EX_m01965x", 10.5 * (Gluc / (0.0027 + Gluc)) * (1 / (1 / 20.0)), Gluc
)
# dfba_model.add_exchange_flux_lb(
#     "Glutamine",
#     6.0 * (Gln / (0.0165 + Gln)) * (1 / (1 + Eth / 20.0)) * (1 / (1 + Gluc / 0.005)),
#     Gln,
# )
dfba_model.add_exchange_flux_lb("Vpalm", (2.95e-5 * (AcoA/(1e-2 + AcoA)) - Palm))#  * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp))))\
#     - 7.35e-5 * mu\
#     - mu * Palm)

# add initial conditions for kinetic variables in dfba_model biomass (gDW/L),
# metabolites (g/L)
dfba_model.add_initial_conditions(
    {
        "Biomass": 1.04e-4,
        "Glucose": 1,
        "Glutamine": 0,
        "Lactate": 0,
        "Palmitate": 0
    }
)


# simulate model across interval t = [0.0,25.0](hours) with outputs for plotting
# every 0.1h and optional list of fluxes
concentrations, trajectories = dfba_model.simulate(
    0.0, 25.0, 0.1, ["biomass_components", "LDH", "DNA", "Vpalm"] #Glc, Xyl, Eth
)


# generate plots of results (in this case using plotlly)
from dfba.plot.plotly import *

import plotly.io as pio

pio.templates.default = "plotly_white"

fig = plot_concentrations(concentrations)
fig.show()

fig = plot_trajectories(trajectories)
fig.show()

# write results to file

concentrations.to_csv("concentrations.csv")
trajectories.to_csv("trajectories.csv")


print("--- %s seconds ---" % (time.time() - start_time))