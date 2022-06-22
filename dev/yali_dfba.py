import cobra
import sys
from dynamicFBA import dynamicFBA

file_path = sys.argv[1]
#Init cobrapy

cobra_config = cobra.Configuration()
# cobra_config.solver = "cplex"

model = cobra.io.read_sbml_model(file_path)

# Setting glucose uptake
model.reactions.get_by_id('1714').lower_bound = -0.60
#print(model.reactions.get_by_id('1714').lower_bound)

# Adding a reaction
reaction  = cobra.Reaction('AOX')
reaction.upper_bound = 0.01
reaction.add_metabolites({model.metabolites.get_by_id('m468'):-1,\
    model.metabolites.get_by_id('m64'):-0.25,\
    model.metabolites.get_by_id('m471'):1,\
    model.metabolites.get_by_id('m26'):0.5})

# print(reaction.reaction)
# model.add_reaction('m468[C_mi] + 0.25 m64[C_mi] => m471[C_mi] + 0.5 m26[C_mi]')

model.add_reaction(reaction)
# print(model.reactions.get_by_id('AOX').bounds)

model.objective = {model.reactions.get_by_id('2111'): 1 }
print(model.objective)
print(model.reactions.get_by_id('2111'))

# uptake = model.lower_bound
# print(uptake)

# Configuration for function dynamicFBA
substrateRxns = ['1654', '1714'] # NH4 - glucose
initConcentration = [55, 288] # NH4 and Glucose concentration (mM)
initBiomass = 0.06 # initial biomass concentration (gDCW/L)
timeStep = 1
nSteps = 171

plotRxns = ['1654','1714','1687','EXC_OUT_isocitrate'] # plot NH4 - Glucose - Citrate
