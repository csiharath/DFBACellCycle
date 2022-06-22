from cobra.test import create_test_model
import cobra

cobra_config = cobra.Configuration()

# Define solver to use
cobra_config.solver = "glpk"

# model = create_test_model("textbook")

# Load a sbml model
# model = cobra.io.read_sbml_model("/home/csiharath/Stage/ecoli_core_model.xml", use_fbc_package=False)
# cobra.io.write_sbml_model(model, "output")

model = cobra.io.read_sbml_model("/home/csiharath/Stage/iJR904.xml")


print(f"{len(model.reactions)} Reactions in the model\n------------------")
for r in model.reactions[:10]:
    print(f"{r.id}: {r.reaction}, associated with {r.gene_name_reaction_rule}")
print()

print(f"{len(model.metabolites)} Metabolites in the model\n------------------")
for m in model.metabolites[:10]:
    print(f"{m.id}: {m.name}")
print()

print(f"{len(model.genes)} Genes in the model\n------------------")
for g in model.genes[:10]:
    print(f"{g.id}: {g.name}")

# Param pour G1 :
model.reactions.get_by_id('PFK').lower_bound = -1000
model.reactions.get_by_id('PFK').upper_bound = 1000
print(model.reactions.get_by_id('PFK').bounds)
#FB3? 
model.reactions.get_by_id('G6PDH2r').lower_bound = 0
model.reactions.get_by_id('G6PDH2r').upper_bound = 0
model.reactions.get_by_id('TKT1').lower_bound = 0
model.reactions.get_by_id('TKT1').upper_bound = 0
# model.reactions.get_by_id('TKT2').lower_bound = 0
# model.reactions.get_by_id('TKT2').upper_bound = 0
model.reactions.get_by_id('G6PDH2r').lower_bound = 0
model.reactions.get_by_id('G6PDH2r').upper_bound = 0
#enzyme lipidique

# Model objective
print(f"objective : {model.objective}")

# Stoichiometric matrix of the studied model
stoichmat = cobra.util.create_stoichiometric_matrix(model)
print(stoichmat.shape)


solution = model.optimize()
print(solution.objective_value)
print(solution.status)
# print(solution.fluxes)
# print(solution.shadow_prices)
# print(model.summary())