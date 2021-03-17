"""
Created on 11/04/2014

@author: Mark N. Read


Module represents a bacteria cell in the model. There is an abstract class, `Bacterium`, from which more concrete
bacteria classes inherit.
"""

from abc import ABCMeta, abstractmethod
import random
from gutsim.nutrients import Nutrients
import math
import numpy as np

# ================= BACTERIA PARAMETERS ====================
# When true, bacterium nutrient uptake rates are positively related to limiting nutrient level.
variable_uptake_rate = False
# Maximum rate at which a bacterium can internalize nutrient from the environment.
max_uptake_rate = 500.0
# The weight content of carbon to nitrogen in a bacterium cell.
# During the bacteria exponentially growing phase, this has been found to be 5.2g carbon for every 1g nitrogen.
# See Vrede K, Heldal M, Norland S, Bratbak G. Elemental composition (C, N, P) and cell volume of
# exponentially growing and nutrient-limited bacterioplankton. Appl Environ Microb. 2002;68(6):2965-71.
# Further support found in 
# Chrzanowski, T. H., Kyle, M., Elser, J. J., & Sterner, R. W. (1996). Element ratios and growth dynamics of bacteria 
# in an oligotrophic Canadian shield lake. Aquatic Microbial Ecology, 11(2), 119–125. 
CtoNRatio = 5.2
# Used in calculating internalised nutrient decay (utilization) rates.
mean_internalised_carbon_lifetime = 5.0  # In hours
mean_internalised_nitrogen_lifetime = 5.0  # In hours


def _calculate_death_rate_exp(limit):
    """
    Calculates the death rate of a bacteria cell with the given limiting nutrient level. As limiting resource level
    tends towards zero, the death rate increases, and then falls to reflect spore formation.
    :param limit: level of limiting nutrient
    :return:
    """
    d = np.exp(-limit * (1.0 / 5.0))
    s = np.exp(-limit * (1.0 / 3.0))
    die = d - s
    return die


def _calculate_divide_rate_logistic(limit):
    """
    Logistic growth curve (s-shaped), remains below 1.0, and above 0.0 at all times
    """
    steepness = 1.0 / 20.0   # dictates the gradient. Small values give shallower gradients.
    shift = 90			# shifts the curve relative to the x-axis. Positive values move curve right.
    return 1.0 / (1.0 + np.exp( -steepness * (limit - shift)))


death_rate_model = _calculate_death_rate_exp

# Lookup table for relating death rate to limiting resource level. 
# The equation is computationally expensive, hence the use of a lookup table.
limiting_resource_resolution = 0.01  # Resolution of limiting resource for which death rates are calculated.
max_limiting_resource_lookup_value = 500  # Maximum value for which lookup table values should be calculated.
death_rate_lookup_table = []  # The lookup table itself.
divide_rate_lookup_table = []
# Populate lookup tables.
for limit in np.arange(0.0, max_limiting_resource_lookup_value, limiting_resource_resolution):
    death_rate_lookup_table.append(death_rate_model(limit))
    divide_rate_lookup_table.append(_calculate_divide_rate_logistic(limit))


def death_rate_lookup(limit):
    """
    Uses look up table to estimate death rate based on limiting resource.
    """
    if limit > max_limiting_resource_lookup_value:
        print('WARNING: limiting resource value exceeds range of pre-calculated values : ' + str(limit))
        return death_rate_model(limit)
    scaled_limit = int(limit / limiting_resource_resolution)
    return death_rate_lookup_table[scaled_limit]


def divide_rate_lookup(limit):
    if limit > max_limiting_resource_lookup_value:
        print('WARNING: limiting resource value exceeds range of pre-calculated values : ' + str(limit))
        return _calculate_divide_rate_logistic(limit)
    scaled_limit = int(limit / limiting_resource_resolution)
    return divide_rate_lookup_table[scaled_limit]


def load_parameters(root):
    """
    Loads parameters specific to this module.
    :param root: an xml.etree.ElementTree Root object.
    """
    global max_uptake_rate
    global CtoNRatio
    global limitExcessiveInternalisation
    global alternativeNutrientUtilisation
    global mean_internalised_carbon_lifetime
    global mean_internalised_nitrogen_lifetime
    global variable_uptake_rate
    node = root.find('./bacteria/maxInternalisationRate')
    if node is not None and node.text is not None:
        max_uptake_rate = float(node.text)
    node = root.find('./bacteria/CtoNRatio')
    if node is not None and node.text is not None:
        CtoNRatio = float(node.text)
    node = root.find('./bacteria/mean_internalised_carbon_lifetime')
    if node is not None and node.text is not None:
        mean_internalised_carbon_lifetime = float(node.text)
    node = root.find('./bacteria/mean_internalised_nitrogen_lifetime')
    if node is not None and node.text is not None:
        mean_internalised_nitrogen_lifetime = float(node.text)
    node = root.find('./bacteria/variable_uptake_rate')
    if node is not None and node.text is not None:
        variable_uptake_rate = bool(int(node.text))


class Bacterium:
    """
    Abstract class that implements the core functionality of bacteria. Concrete sub-classes provide the specifics.
    """
    __metaclass__ = ABCMeta

    def __init__(self, timestep, carbon=130.0, nitrogen=25.0):
        """
        :param timestep: time delta for each simulation step, in hours.
        :param carbon: the initial quantity of carbon that this bacterium contains.
        :param nitrogen: initial quantity of nitrogen that this bacterium contains.
        :return:
        """
        self._timestep = timestep
        self._carbon = carbon
        self._nitrogen = nitrogen
        self.death_rate = 0.0
        self.divide_rate = 0.0
        self._limiting_resource = self._calculate_limiting_resource()
        self.dead = False

    @abstractmethod
    def _extract_nutrients(self, available):
        """
        Extracts nutrients from the environment and into the bacteria cell. This abstract method must be made
        concrete by subclasses to select specific nutrients from what is generally available in the environment.
        """

    @abstractmethod
    def _spawn_child(self, carbon, nitrogen):
        """ Returns another bacteria object of the same type as this one. """

    def _calculate_limiting_resource(self):
        """
        A certain ratio of carbon to nitrogen is required for growth. This calculates which of the two resources
        is the limiting resource, given the required ratio, and expresses it in terms of nitrogen. 
        E.g. if a cell needs twice as much carbon (C) as nitrogen (N), and has
        - C = 10, N = 4, then the limiting resource is 4, because nitrogen is limiting.
        - C = 15, N = 10, then the limiting resource is 7.5, because carbon is limiting.
        """
        ln = self._n_limit_equiv()
        lc = self._c_limit_equiv()
        return min(ln, lc)

    def _calculate_uptake_rate(self):
        """
        Stressed bacteria cells have a reduced capacity to uptake nutrients from their environments. 
        Function returns a value between 0 and max_uptake_rate.
        """
        rate = max_uptake_rate
        if variable_uptake_rate:
            steepness = 1.0 / 20.0
            rate = max_uptake_rate / (1.0 + np.exp(-steepness * (self._limiting_resource - 40.0)) )
        return rate

    def _c_limit_equiv(self):
        """ Returns quantity of carbon, as the limiting resource equivalent (which is just N at present). """
        return self._carbon / CtoNRatio

    def _n_limit_equiv(self):
        """ Returns quantity of nitrogen as a limiting resource equivalent (which is just N at present). """
        return self._nitrogen

    def _interalise_choices(self):
        """
        Returns two boolean values as a tuple, representing whether carbon and nitrogen sources can be internalized.
        This is based on whether the internalised quantity of each is less than 1.3 times the equivalent (given carbon
        to nitrogen ratio) of the other. 
        This is used to prevent, e.g., a carbon-limited bacterium consuming limitless quantities of nitrogen.
        """
        limN = self._n_limit_equiv()
        limC = self._c_limit_equiv()
        internaliseN = limN < 1.30 * limC
        internaliseC = limC < 1.30 * limN
        return internaliseN, internaliseC

    def _decide_divide(self, limit):
        """
        Returns true when a cell is to divide. Does not handle any of the mechanics of dividing, simply the decision
        on whether or not to do so.
        """
        self.divide_rate = divide_rate_lookup(limit)
        p = self.divide_rate * self._timestep
        return random.random() < p

    def _divide(self):
        """ Performs cellular division. """
        # Divide the nutrient contents of this cell evenly between it and the child bacterium.
        self._carbon /= 2.0
        self._nitrogen /= 2.0
        # Generate another bacterium of the same type as this one.
        child = self._spawn_child(carbon=self._carbon, nitrogen=self._nitrogen)
        # Half the nutritional uptake machinery goes with it, hence halve that too.
        return child

    def _decide_die(self, limit):
        """
        Returns true when a cell is to die. Does not handle any of the mechanics of integrating death into the
        simulation, simply the decision on whether or not a cell dies now.
        """
        self.death_rate = death_rate_lookup(limit)  # Use a lookup table, for efficiency.
        # Scale by timestep, since the probability represents the chance of dying in the last hour. 
        # If the timestep changes, the probability should not.
        p = self.death_rate * self._timestep
        die = random.random() < p
        return die

    def _perform_die(self):
        """
        Carries out the nuts and bolts of actually dying, ie, telling the vessel that this bacteria no longer needs
        to exist.
        """
        self.dead = True

    def _internalise_nutrients(self, available):
        """
        Internalises nutrients extracted from the environment (which is handled by another function) into the cell.
        Arguments:
          available - Nutrients object, the nutrients visible to the cell.
        """
        consumed = self._extract_nutrients(available)
        self._carbon += consumed.carbon_content()
        self._nitrogen += consumed.nitrogen_content()
        return consumed

    def _nutrient_decay(self):
        """
        Handles the ongoing decay of nutrients held within the cell
        """
        # Exponential decay. dN/dt = -delta . N
        # Mean lifetime (ML) used in this equation.
        # Hence, death rate (delta) = 1 / (ML).
        carbon_decay_rate = 1.0 / mean_internalised_carbon_lifetime
        nitrogen_decay_rate = 1.0 / mean_internalised_nitrogen_lifetime
        decayedC = self._carbon * carbon_decay_rate * self._timestep
        decayedN = self._nitrogen * nitrogen_decay_rate * self._timestep
        self._carbon -= decayedC
        self._nitrogen -= decayedN

    def get_carbon_nitrogen(self):
        """
        Returns a tuple, representing the total carbon and nitrogen resources of the cell (in that order).
        """
        return self._carbon, self._nitrogen

    def step(self, environment):
        """
        Called within the mouse simulator to step the bacteria through its states.
        """
        if self.dead:
            return Nutrients(), None  # Dead bacteria don't do anything.
        available = environment.clone()

        self._nutrient_decay()  # The cell consumes nutrients, at a decreasing rate (as it enters stress response).
        internalised = self._internalise_nutrients(available)  # Cell consumes resources from the environment.
        self._limiting_resource = self._calculate_limiting_resource()
        child = None
        if self._decide_divide(self._limiting_resource):
            child = self._divide()
        if self._decide_die(self._limiting_resource):
            self._perform_die()
        return internalised, child

    def string_report(self):
        return 'Int. nut.: carbon=' + str(self._carbon) + ', nitrogen=' + str(self._nitrogen)

    def carbon_limiting(self):
        """
        Returns True when carbon is the limiting resouce in this bacterium.
        """
        # Takes account of the difference in a bacterium's carbon to nitrogen requirements. 
        # Cell requires more carbon atoms than nitrogen atoms.
        return self._carbon < (CtoNRatio * self._nitrogen)

    def nitrogen_limiting(self):
        """ Returns True when nitrogen is the limiting resource in this bacterium. """
        return not self.carbon_limiting()

    def is_dead(self):
        return self.dead

    def is_exponential(self):
        """ Refers to the exponential growth stage. """
        return self._limiting_resource >= 20.0 and not self.is_dead()

    def is_stressed(self):
        """ Refers to the stressed response. Note this is not the long term stress response. """
        return 1.0 <= self._limiting_resource < 20.0 and not self.is_dead()

    def is_resistant(self):
        """ Refers to the long term stressed response, at which point bacteria are highly resistant to death. """
        return self._limiting_resource < 1 and not self.is_dead()


class Guild_rf_f(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_rf_f(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()  # Which elements can be internalised?
        if internaliseN:
            consumed.Nf = min(available.Nf, uptake_quantity)
        if internaliseC:
            consumed.Cw = min(available.Cw, uptake_quantity)
            consumed.Ci = min(available.Ci, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_rf_f: ' + self.string_report()


class Guild_rm_m(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_rm_m(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()  # Which elements can be internalised?
        if internaliseN:
            consumed.Nm = min(available.Nm, uptake_quantity)
        if internaliseC:
            consumed.Cw = min(available.Cw, uptake_quantity)
            consumed.Ci = min(available.Ci, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_rm_m: ' + self.string_report()


class Guild_pf_f(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_pf_f(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()  # Which elements can be internalised?
        if internaliseN:
            consumed.Nf = min(available.Nf, uptake_quantity)
        if internaliseC:
            consumed.Cd = min(available.Cd, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_pf_f: ' + self.string_report()


class Guild_pm_m(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_pm_m(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()	 # Which elements can be internalised?
        if internaliseN:
            consumed.Nm   = min(available.Nm, uptake_quantity)
        if internaliseC:
            consumed.Cd   = min(available.Cd, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_pm_m: ' + self.string_report()


class Guild_m_m(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_m_m(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()	 # Which elements can be internalised?
        if internaliseN:
            consumed.Nm = min(available.Nm, uptake_quantity)
        if internaliseC:
            consumed.Cm = min(available.Cm, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_m_m: ' + self.string_report()


class Guild_mf_f(Bacterium):
    def _spawn_child(self, carbon, nitrogen):
        return Guild_mf_f(timestep=self._timestep, carbon=carbon, nitrogen=nitrogen)

    def _extract_nutrients(self, available):
        consumed = Nutrients()
        uptake_quantity = self._calculate_uptake_rate() * self._timestep
        internaliseN, internaliseC = self._interalise_choices()	 # Which elements can be internalised?
        if internaliseN:
            consumed.Nf = min(available.Nf, uptake_quantity)
        if internaliseC:
            consumed.Cm = min(available.Cm, uptake_quantity)
        return consumed

    def __str__(self):
        return 'Guild_mf_f: ' + self.string_report()