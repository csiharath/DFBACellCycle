"""

"""
import logging
from typing import Dict, List, Optional, Tuple
from cobra import DictList, Model, Object, Reaction
from optlang import symbolics

logger = logging.getLogger(__name__)

class DfbaModel(Object):
    """
    Class representation for a dynamic FBA model.
    """

    def __init__(self, cobra_object: Model) -> None:
        """Init method."""
        if isinstance(cobra_object, Model):
            self._id = id(self)
            self._cobra_model = cobra_object
            self._reactions = cobra_object.reactions
            self._objectives = []
            self._directions = []
            self._kinetic_variables = DictList()
            self._exchange_fluxes = DictList()
        else:
            raise Exception(
                "Error: must use instance of class cobra.Model " "to init DfbaModel!"
            )
    
    def add_objectives(self, objectives: List, directions: List) -> None:
        """Add objectives.

        Parameters
        -------
        objectives : list
            The list of reaction indetifiers to be added to the model as
            objectives for lexicographic optimization.
        directions : list
            The list of directions (max or min) of each objective to be
            added to the model for lexicographic optimization.


        """
        # TODO: currently not supported. Should it raise NotImplementedError?

        if type(objectives) is not list:
            objectives = [objectives]
        if type(directions) is not list:
            directions = [directions]
        if len(objectives) != len(directions):
            raise Exception(
                "Error: objectives list must be same length as directions list!"
            )
        self._objectives.extend(objectives)
        self._directions.extend(directions)


    def add_kinetic_variables(self, kinetic_variable_list: List) -> None:
        """Add kinetic variables.

        Parameters
        -------
        kinetic_variable_list : list
            The list of indetifiers of kinetic variables to be added to the
            model.

        """
        if type(kinetic_variable_list) is not list:
            kinetic_variable_list = [kinetic_variable_list]

        def existing_filter(kinetic_variable):
            if kinetic_variable.id in self.kinetic_variables:
                logger.warning(
                    f"Ignoring kinetic variable {kinetic_variable.id} since it"
                    f"already exists in the model."
                )
                return False
            return True

        pruned = DictList(filter(existing_filter, kinetic_variable_list))
        self._kinetic_variables += pruned
        DictList.sort(self._kinetic_variables)
        counter = 0
        for kinetic_variable in self._kinetic_variables:
            symbolics.Symbol.__init__(kinetic_variable, "yval[" + str(counter) + "]")
            counter += 1

    def add_exchange_fluxes(self, exchange_flux_list: List) -> None:
        """Add exchange fluxes.

        Parameters
        -------
        exchange_flux_list : list
            list of identifiers of exchange fluxes to be added to the model.

        """
        if type(exchange_flux_list) is not list:
            exchange_flux_list = [exchange_flux_list]

        def existing_filter(exchange_flux: ExchangeFlux) -> bool:
            if exchange_flux.id not in self.reactions:
                logger.warning(
                    f"Ignoring exchange flux {exchange_flux.id} since it does"
                    f"not correspond to a reaction in the model."
                )
                return False
            if exchange_flux.id in self.exchange_fluxes:
                logger.warning(
                    f"Ignoring exchange flux {exchange_flux.id} since it"
                    f"already exists in the model."
                )
                return False
            return True

        pruned = DictList(filter(existing_filter, exchange_flux_list))
        self._exchange_fluxes += pruned
        DictList.sort(self._exchange_fluxes)
        counter = 0
        for exchange_flux in self._exchange_fluxes:
            symbolics.Symbol.__init__(exchange_flux, "v" + str(counter))
            counter += 1