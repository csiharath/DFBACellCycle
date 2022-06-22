import cobra
from math import exp
import numpy as np

def dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns=['EX_glc(e)','EX_ac(e)','EX_for(e)'], exclUptakeRxns = ['CO2 transport','O2 transport','H2O transport','H+ transport']):
    """
    Performs dynamic FBA simulation using the static optimization approach

    usage:

        dynamicFBA(model, substrateRxns, initConcentrations, initBiomass, timeStep, nSteps, plotRxns, exclUptakeRxns)
    """

    ############################ Find exchange reactions #############################

    excRxns = cobra.medium.boundary_types.find_boundary_types(model,'exchange') # Gets a list of all exchange reactions
    RxnsNames = [rxn.name for rxn in excRxns] # Creates a list of the reactions names
    # print(RxnsNames)
    keepRxns = [not(i) for i in np.isin(RxnsNames, exclUptakeRxns)] # Creates an array of booleans with reaction to keep
    excRxnsNames = [excRxns[rxn] for rxn in range(len(RxnsNames)) if keepRxns[rxn]] # Adapts list with the booleans array
    # print(len(excRxnsNames))

    ################# Figure out if substrate reactions are correct ##################



    ############################ Initialize concentrations ###########################

    subRxn = [model.reactions.get_by_id(rxn) for rxn in substrateRxns] # List of reactions for substrate initially present
    keepsub = np.isin(excRxnsNames, substrateRxns) # Creates an array of booleans with true if it's a reaction for wixhe substrate is initially here
    substrateMatchInd = []
    for ind in range(len(excRxnsNames)):
        if keepsub[ind]:
            # print(subRxn.index(excRxnsNames[ind]))
            substrateMatchInd.append(subRxn.index(excRxnsNames[ind]) + 1)
        else:
            substrateMatchInd.append(0)

    concentration = [0 for i in range(len(excRxnsNames))]
    lb_list = [rxn.lower_bound for rxn in excRxnsNames]

    # Initialize concentration accordingly to the index
    for ind in range(len(substrateMatchInd)):
        if substrateMatchInd[ind] != 0:
            concentration[ind] += initConcentrations[substrateMatchInd[ind]-1]
        elif abs(lb_list[ind]) > 0: # Deal with reactions for which there are no initial concentrations
            concentration[ind] += 1000

    biomass = initBiomass
    
    ############################### Initialize bounds ################################
   
    uptakebound = [c/biomass*timeStep for c in concentration]

    ###### Make sure bounds are not higher than what are specified in the model ######

    biomassVec = [biomass]
    fluxesMat = []

    for step in range(nSteps):
        sol = model.optimize()
        mu = sol.objective_value
        if (sol.status != 'optimal' or mu == 0):
            model.reactions.get_by_id('1714').lower_bound = -0.39 # Glucose uptake
            model.reactions.get_by_id('1992').lower_bound = -0.217 #O2
            model.reactions.get_by_id('1687').bounds = 0.1, 0.1 # citrate
            model.reactions.get_by_id('AOX').upper_bound = 0.014 #AOX reaction

            model.objective = {model.reactions.get_by_id('xLIPID'): 1 }
            sol = model.optimize()

        
        fluxes = sol.fluxes
        fluxesMat.append(fluxes)
        #uptakeFlux = sol.fluxes()
        biomass = biomass*exp(mu*timeStep)
        biomassVec.append(biomass)

        # Update concentrations
        # if mu != 0:
        #     concentrations = concentrations - uptakeFlux / (mu*biomass*(1-exp(mu*timeStep)))
        # else:
        #     concentrations = concentrations + uptakeFlux * biomass * timeStep
        