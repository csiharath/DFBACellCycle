"""

Organism -> _Escherichia coli str. K-12 substr. MG1655_
Model in http://bigg.ucsd.edu/models/iJR904
"""

from lib2to3.pgen2.token import AT
from cobra.io import read_sbml_model

from dfba import DfbaModel, ExchangeFlux, KineticVariable


# DfbaModel instance initialized with cobra model
fba_model = read_sbml_model("iJR904.xml")
fba_model.solver = "glpk"
dfba_model_G1 = DfbaModel(fba_model)

# instances of KineticVariable (default initial conditions are 0.0)
G6p = KineticVariable("Glucose 6-phosphate")
F6p = KineticVariable("Fructose 6-phosphate")
Gap = KineticVariable("Glyceraldehyde 3-phosphate")
Pep = KineticVariable("Phosphoenolpyruvate")
Pyr = KineticVariable("Pyruvate")
Lac = KineticVariable("Lactate")
R5p = KineticVariable("ribose 5-phosphate")
X5p = KineticVariable("xylulose 5-phosphate")
AcoA = KineticVariable("Acetyl-CoA")
Cit = KineticVariable("Citrate")
Akg = KineticVariable("2-Oxoglutarate")
Suc = KineticVariable("Succinate")
Mal = KineticVariable("Malate")
Oxa = KineticVariable("Oxaloacetate")
Palm = KineticVariable("Palmitate")
Glu = KineticVariable("Glutamate")
Ala = KineticVariable("Alanine")
Nh4 = KineticVariable("Ammonium")
Glu_ex = KineticVariable("Glutamate ext")
Atp = KineticVariable("ATP", initial_condition=2.63)
Adp = KineticVariable("ADP", initial_condition=0.28)
Amp = KineticVariable("AMP", initial_condition=0.197)
Nad = KineticVariable("NAD", initial_condition=0.44)
Nadh = KineticVariable("NADH", initial_condition=5.41)
Nadp = KineticVariable("NADP", initial_condition=0.222)
Nadph = KineticVariable("NADPH", initial_condition=0.146)
X = KineticVariable("Biomass", initial_condition=1.04e-4)


# add kinetic variables to dfba_model_G1
dfba_model_G1.add_kinetic_variables([G6p, F6p, Gap, Pep, Pyr, Lac, R5p, X5p, AcoA, Cit, Akg, Suc, Mal, Oxa, Palm, Glu, Ala, Nh4, Glu_ex, Atp, Adp, Amp, Nad, Nadh, Nadp, Nadph, X])

# instances of ExchangeFlux
mu = ExchangeFlux("BiomassEcoli")
# v_GP = ExchangeFlux("EX_g6p(e)")
# v_P = ExchangeFlux("EX_pyr(e)")
# v_L = ExchangeFlux("EX_lac(e)") # Lac_D lac_L ?
# v_C = ExchangeFlux("EX_cit(e)")
# v_AK = ExchangeFlux("EX_akg(e)")
# v_S = ExchangeFlux("EX_succ(e)")
# v_M = ExchangeFlux("EX_mal_L(e)") #ou just mal ?
# v_G = ExchangeFlux("EX_glu_L(e)") # ou juste glu ?
# v_Al = ExchangeFlux("EX_ala")# D ou L ?
# v_N = ExchangeFlux("EX_nh4(e)")
# v_Am = ExchangeFlux("EX_amp(e)")
# v_Nd = ExchangeFlux("EX_nad(e)")

v_O2 = ExchangeFlux("Ex_v_O2(e)")
v_Gln = ExchangeFlux("Ex_Gln(e)")
v_Glc = ExchangeFlux("Ex_Glc(e)")

# add exchange fluxes to dfba_model_G1
# dfba_model_G1.add_exchange_fluxes([mu, v_GP, v_P, v_L, v_C, v_AK, v_S, v_M, v_G, v_Al, v_N, v_Am, v_Nd])
dfba_model_G1.add_exchange_fluxes([mu])

# add rhs expressions for kinetic variables in dfba_model_G1

# Donner la formule de Mu comme dans les données du papier????? 
dfba_model_G1.add_rhs_expression("Biomass", mu * X)
dfba_model_G1.add_rhs_expression("Glucose 6-phosphate", (1e-2 * ((v_Glc * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (3 * (1 + (1/9.04e-2) * (Amp/Atp)) + v_Glc * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 1.09e-2/(1.09e-2 + G6p))\
    - (7.6e-4 * (G6p/(1.83e-2 + G6p)) * (2.357e-1/(2.357e-1 + Pep)) - 3.5e-4 * (F6p/(1e-2 + F6p)))\
    - (1.5e-5 * (G6p/(1.83e-2 + G6p)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - 8.79e-5 * mu - mu * G6p)
dfba_model_G1.add_rhs_expression("Fructose 6-phosphate", (7.6e-4 * (G6p/(1.83e-2 + G6p)) * (2.357e-1/(2.357e-1 + Pep)) - 3.5e-4 * (F6p/(1e-2 + F6p)))\
    - (1.5e-3 * ((F6p * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (1e-2 * (1 + (1/9.04e-2) * (Amp/Atp)) + F6p * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 471.3376/(471.3376 + Cit))\
    + 2 * (2.3e-5 * (R5p/(4.46e-2 + R5p)) * (X5p/(7e-3 + X5p)))\
    - mu * F6p)
dfba_model_G1.add_rhs_expression("Glyceraldehyde 3-phosphate", 2 * (1.5e-3 * ((F6p * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (1e-2 * (1 + (1/9.04e-2) * (Amp/Atp)) + F6p * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 471.3376/(471.3376 + Cit))\
    - (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    + (2.3e-5 * (R5p/(4.46e-2 + R5p)) * (X5p/(7e-3 + X5p)))\
    - mu * Gap)
dfba_model_G1.add_rhs_expression("Phosphoenolpyruvate", (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - (1.4e-3 * ((Pep * (1 + (1.745/(4.1113 * 3.121)) * F6p)) / (5e-2 * (1 + (1/3.121) * F6p) + Pep * (1 + (1/(4.1113 * 3.121)) * F6p))) * ((Adp/Atp)/(0.2 + (Adp/Atp))))\
    - mu * Pep)
dfba_model_G1.add_rhs_expression("Pyruvate", (1.4e-3 * ((Pep * (1 + (1.745/(4.1113 * 3.121)) * F6p)) / (5e-2 * (1 + (1/3.121) * F6p) + Pep * (1 + (1/(4.1113 * 3.121)) * F6p))) * ((Adp/Atp)/(0.2 + (Adp/Atp))))\
    - (6e-9 * ((Pyr * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (0.1 * (1 + (1/9.04e-2) * (Amp/Atp)) + Pyr * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Nadh/Nad)/(8.3 + Nadh/Nad)))\
    - (7.1e-5 * (Pyr/(0.1 + Pyr)) * ((Nad/Nadh)/4.8e-2 * (Nad/Nadh)))\
    - (1.76e-4 * (Pyr/(0.1 + Pyr)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    + (1.3e-5 * (Mal/(0.5 + Mal)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - (1.7e-4 * (Glu/(1e-2 + Glu)) * (Pyr/(0.1 + Pyr)))\
    - mu * Pyr)
dfba_model_G1.add_rhs_expression("Lactate", (6e-9 * ((Pyr * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (0.1 * (1 + (1/9.04e-2) * (Amp/Atp)) + Pyr * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Nadh/Nad)/(8.3 + Nadh/Nad)))\
    - mu * Lac)
dfba_model_G1.add_rhs_expression("ribose 5-phosphate", (1.5e-5 * (G6p/(1.83e-2 + G6p)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - (1.3e-5 * (R5p/(4.46e-2 + R5p)))\
    - (2.3e-5 * (R5p/(4.46e-2 + R5p)) * (X5p/(7e-3 + X5p)))\
    - 7.35e-5 * mu\
    - mu * R5p)
dfba_model_G1.add_rhs_expression("xylulose 5-phosphate", (1.3e-5 * (R5p/(4.46e-2 + R5p)))\
    -2 * (2.3e-5 * (R5p/(4.46e-2 + R5p)) * (X5p/(7e-3 + X5p)))\
    - mu * X5p)
dfba_model_G1.add_rhs_expression("Acetyl-CoA", (7.1e-5 * (Pyr/(0.1 + Pyr)) * ((Nad/Nadh)/4.8e-2 * (Nad/Nadh)))\
    - (8.8e-5 * (AcoA/(1e-2 + AcoA)) * (Oxa/(1.401e-1 + Oxa)))\
    + (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - 8 * (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp))))\
    - mu * AcoA)
dfba_model_G1.add_rhs_expression("Citrate", (8.8e-5 * (AcoA/(1e-2 + AcoA)) * (Oxa/(1.401e-1 + Oxa)))\
    - (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - (3.9e-5 * (Cit/(1.019e-1 + Cit)) * ((Nad/Nadh)/4.8e-2 * (Nad/Nadh)))\
    - mu * Cit)
dfba_model_G1.add_rhs_expression("2-Oxoglutarate", (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    + (1.2e-6 * (Glu/(1e-2 + Glu)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    + (1.7e-4 * (Glu/(1e-2 + Glu)) * (Pyr/(0.1 + Pyr)))\
    - mu * Akg)
dfba_model_G1.add_rhs_expression("Succinate", (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    - (4e-4 * (Suc/(1.465e-1 + Suc)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - mu * Suc)
dfba_model_G1.add_rhs_expression("Malate", (4e-4 * (Suc/(1.465e-1 + Suc)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - (9.5e-5 * (Mal/(0.5 * Mal)) * (((Nad/Nadh)/(4.8e-2 + Nad/Nadh))))\
    - (1.3e-5 * (Mal/(0.5 + Mal)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - mu * Mal)
dfba_model_G1.add_rhs_expression("Oxaloacetate", - (8.8e-5 * (AcoA/(1e-2 + AcoA)) * (Oxa/(1.401e-1 + Oxa)))\
    + (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    + (9.5e-5 * (Mal/(0.5 * Mal)) * (((Nad/Nadh)/(4.8e-2 + Nad/Nadh))))\
    + (1.76e-4 * (Pyr/(0.1 + Pyr)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - mu * Oxa)
dfba_model_G1.add_rhs_expression("Palmitate", (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp))))\
    - 7.35e-5 * mu\
    - mu * Palm)
dfba_model_G1.add_rhs_expression("Glutamate", (2.54e-4 * (v_Gln/(0.1 + v_Gln)))\
    - (1.8e-6 * (Glu/(1e-2 + Glu)))\
    - (1.2e-6 * (Glu/(1e-2 + Glu)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - (1.7e-4 * (Glu/(1e-2 + Glu)) * (Pyr/(0.1 + Pyr)))\
    - mu * Glu)
dfba_model_G1.add_rhs_expression("Alanine", (1.7e-4 * (Glu/(1e-2 + Glu)) * (Pyr/(0.1 + Pyr)))\
    - mu * Ala)
dfba_model_G1.add_rhs_expression("Ammonium", (2.54e-4 * (v_Gln/(0.1 + v_Gln)))\
    + (1.2e-6 * (Glu/(1e-2 + Glu)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - mu * Nh4)
dfba_model_G1.add_rhs_expression("Glutamate ext", (1.8e-6 * (Glu/(1e-2 + Glu)))\
    - mu * Glu_ex)
dfba_model_G1.add_rhs_expression("ATP", - (1e-2 * ((Glu_ex * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (3 * (1 + (1/9.04e-2) * (Amp/Atp)) + Glu_ex * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 1.09e-2/(1.09e-2 + G6p))\
    - (1.5e-3 * ((F6p * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (1e-2 * (1 + (1/9.04e-2) * (Amp/Atp)) + F6p * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 471.3376/(471.3376 + Cit))\
    + (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - (1.76e-4 * (Pyr/(0.1 + Pyr)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    + (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    - 7 * (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp))))\
    - (9.14e-4 *(Atp/(1 + Atp)))\
    + 4.2 * (1.3e-3 * ((Nadh/Nad)/(8.3 + Nadh/Nad)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * (v_O2/(1e-2 + v_O2)))\
    - (1.8e-8 * (Atp/(1 + Atp)) * (Amp/(1.5e-3 + Amp)) - 1.2e-2 * (Adp/0.2 + Adp))\
    - 1.19e-2 * mu)
dfba_model_G1.add_rhs_expression("ADP", (1e-2 * ((Glu_ex * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (3 * (1 + (1/9.04e-2) * (Amp/Atp)) + Glu_ex * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 1.09e-2/(1.09e-2 + G6p))\
    + (1.5e-3 * ((F6p * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (1e-2 * (1 + (1/9.04e-2) * (Amp/Atp)) + F6p * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Atp/Adp)/(5 + Atp/Adp)) * 471.3376/(471.3376 + Cit))\
    - (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    + (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    + (1.76e-4 * (Pyr/(0.1 + Pyr)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    + 7 * (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp))))\
    + (9.14e-4 *(Atp/(1 + Atp)))\
    - 4.2 * (1.3e-3 * ((Nadh/Nad)/(8.3 + Nadh/Nad)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * (v_O2/(1e-2 + v_O2)))\
    + (1.8e-8 * (Atp/(1 + Atp)) * (Amp/(1.5e-3 + Amp)) - 1.2e-2 * (Adp/0.2 + Adp))\
    + (1.8e-8 * (Atp/(1 + Atp)) * (Amp/(1.5e-3 + Amp)) - 1.2e-2 * (Adp/0.2 + Adp))\
    + 1.19e-2 * mu)
dfba_model_G1.add_rhs_expression("AMP", - (1.8e-8 * (Atp/(1 + Atp)) * (Amp/(1.5e-3 + Amp)) - 1.2e-2 * (Adp/0.2 + Adp)))
dfba_model_G1.add_rhs_expression("NAD", - (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    + (6e-9 * ((Pyr * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (0.1 * (1 + (1/9.04e-2) * (Amp/Atp)) + Pyr * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Nadh/Nad)/(8.3 + Nadh/Nad)))\
    - (7.1e-5 * (Pyr/(0.1 + Pyr)) * ((Nad/Nadh)/4.8e-2 * (Nad/Nadh)))\
    - (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    - (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    - 0.66 * (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    - (9.5e-5 * (Mal/(0.5 * Mal)) * (((Nad/Nadh)/(4.8e-2 + Nad/Nadh))))\
    - (1.2e-6 * (Glu/(1e-2 + Glu)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    + 2 * (1.3e-3 * ((Nadh/Nad)/(8.3 + Nadh/Nad)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * (v_O2/(1e-2 + v_O2)))\
    + 2 * (2.9e-5 * (Nadh/(0.1 + Nadh))))
dfba_model_G1.add_rhs_expression("NADH", (9.91e-4 * (Gap/(0.1 * Gap)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - (6e-9 * ((Pyr * (1 + (10.4557/(4.712e-1 * 9.04e-2)) * (Amp/Atp))) / (0.1 * (1 + (1/9.04e-2) * (Amp/Atp)) + Pyr * (1 + (1/(4.712e-1 * 9.04e-2)) * (Amp/Atp)))) * ((Nadh/Nad)/(8.3 + Nadh/Nad)))\
    + (7.1e-5 * (Pyr/(0.1 + Pyr)) * ((Nad/Nadh)/4.8e-2 * (Nad/Nadh)))\
    + (2.4e-5 * (Cit/(1.019e-1 + Cit)) * ((Atp/Adp)/(5 * (Atp/Adp))))\
    + (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    + 0.66 * (2.4e-4 * (Akg/(7.006e-1 + Akg)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)) * ((Adp/Atp)/(0.2 + Adp/Atp)))\
    + (9.5e-5 * (Mal/(0.5 * Mal)) * (((Nad/Nadh)/(4.8e-2 + Nad/Nadh))))\
    + (1.2e-6 * (Glu/(1e-2 + Glu)) * ((Nad/Nadh)/(4.8e-2 + Nad/Nadh)))\
    - 2 * (1.3e-3 * ((Nadh/Nad)/(8.3 + Nadh/Nad)) * ((Adp/Atp)/(0.2 + Adp/Atp)) * (v_O2/(1e-2 + v_O2)))\
    - 2 * (2.9e-5 * (Nadh/(0.1 + Nadh))))
#Probleme : Deux equations pour nadp ou deuxieme nadph?
dfba_model_G1.add_rhs_expression("NADP", - 2 * (1.5e-5 * (G6p/(1.83e-2 + G6p)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - (1.3e-5 * (Mal/(0.5 + Mal)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    + (1.4e-5 * (Nadph/(5.32e-2 + Nadph)))\
    + 14 * (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp)))))
dfba_model_G1.add_rhs_expression("NADPH", 2 * (1.5e-5 * (G6p/(1.83e-2 + G6p)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    + (1.3e-5 * (Mal/(0.5 + Mal)) * ((Nadp/Nadph)/2 * (Nadp/Nadph)))\
    - (1.4e-5 * (Nadph/(5.32e-2 + Nadph)))\
    - 14 * (2.95e-5 * (AcoA/(1e-2 + AcoA)) * ((Atp/Adp)/(5 * (Atp/Adp))) * ((Nadph/Nadp)/(0.5 * (Nadph/Nadp)))))


# add lower/upper bound expressions for exchange fluxes in dfba_model_G1 together
# with expression that must be non-negative for correct evaluation of bounds
dfba_model_G1.add_exchange_flux_lb("EX_Glc(e)", )
dfba_model_G1.add_exchange_flux_lb("EX_O2(e)", )
dfba_model_G1.add_exchange_flux_lb("EX_Gln(e)", )

# add initial conditions for kinetic variables in dfba_model_G1 biomass (gDW/L),
# metabolites (g/L)
# dfba_model_G1.add_initial_conditions(
#     {
#         "Biomass": 0.03,
#         "Glucose": 15.5,
#         "Xylose": 8.0,
#         "Oxygen": 0.0,
#         "Ethanol": 0.0,
#     }
# )

# simulate model across interval t = [0.0,25.0](hours) with outputs for plotting
# every 0.1h and optional list of fluxes
# concentrations, trajectories = dfba_model_G1.simulate(
#     0.0, 25.0, 0.1, ["EX_glc(e)", "EX_xyl_D(e)", "EX_etoh(e)"]
# )

# # generate plots of results (in this case using plotlly)
# from dfba.plot.plotly import *

# import plotly.io as pio

# pio.templates.default = "plotly_white"

# fig = plot_concentrations(concentrations)
# fig.show()

# fig = plot_trajectories(trajectories)
# fig.show()

# # write results to file

# concentrations.to_csv("concentrations.csv")
# trajectories.to_csv("trajectories.csv")
