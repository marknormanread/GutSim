"""
Created on 11/04/2014

@author: Mark N. Read
"""

import random
from random import shuffle
import gutsim.bacteria
import gutsim.vessel
from gutsim.nutrients import Nutrients
import stats.stats as stats


# Determines the scale of simulation. It defines how much nutrient a single bacteria can internalise. 
# Hence, the scale of observation (or alternatively, how many bacteria  single simulated bacterium really represents).
scale_factor = 50.0 * 1000000.0
# Upper limit on gut (cecum & colon) content. 
# The tradeoff of defecation against microbial growth yields variable dry weight contents. 
# In Alison's DSS study (Nov 2018), cecal weights were recorded as being around 0.5-1.5g. 
# The colon is about half the volume of the cecum, making the MAXIMUM total (modelled, colon + cecum) gut weight 2.25g
# Most of these measurements represent gut contents; the epithelial layers are minimal in weight. 
# These include the water fraction of digesta. Hard to estimate this. 
# Alison's study shows water content of mouse feces to be 20-65% by weight, with 50 being average. 
# The cecum will be higher than this, and accounts for around 2/3 of the total weight under consideration. 
# We are assuming the average water content, by weight, of all cecal & colon contents to be 70%. 
# Adjust accordingly for dry weight.
gut_carrying_dry_weight = 2.25 * 0.3 
# In Alison's DSS study (Nov 2018), individual faecal pellets were recorded as being 15-100mg, usually at 20mg.
# Diets we have been modelling do contain a fair bit of cellulose. Hence, assume 50mg=0.05g. 
# The higher figures are for high-fibre diets, and two pellets can merge together.  
# Hence, each faecal pellet is 3.3% of gut contents, and we think they defecate around 3 times an hour. 
# Note that defecation rates are highly variable with diet, circadian rhythm, and how scared/stressed the animal is. 
# Alison's study shows water content of mouse feces to be 20-65% by weight, with 50 being average. 
# Quantity of GIT contents excreted in an hour (assuming timestep unit is proportion of 1 hour). 
excretion_rate = 0.10  # Proportion of gut contents / hour.
# Length of gut from cecum to anus, in cm. Again, measured from Alison's DSS mouse study. 
# This varied from 6 to 10cm, with the longer values for high fibre diets. 
gut_length = 8.0  # cm
# Units are events per hour. 60 min, peristalsis wave hits every 4.8 minutes. Based on:
# Roberts RR, Murphy JF, Young HM, Bornstein JC. Development of colonic motility in the neonatal mouse-studies using
# spatiotemporal maps. Am J Physiol-Gastr L. 2007;292(3):G930-G8.
peristalsis_rate = (60 / 4.6)  # = ~13 / hour
# Standard deviation of distance bacteria moved by resulting from peristalsis. Given in cm/peristalsis event.
bacteria_peristalsis_rate = 0.15
# Standard deviation of distance nutrients are diffused by, as a result of peristalsis. cm/peristalsis event.
nutrient_peristalsis_rate = 0.15
# Of the nutrients input to the cecum & colon, what proportion remains there? Does not consider small intestine. 
# This organic matter takes the form of microbial biomass and terminal metabolites that are not (and whose precursors
# have not) been absorbed by the host. 
non_absorbable_organic_matter_generation_rate = 0.7

# Initial bacteria population sizes.
init_Grf = 0
init_Grm = 0
init_Gpf = 0
init_Gpm = 0
init_Gmm = 0
init_Gmf = 0
inoculation_rate = 0.000003


def load_parameters(root):
    """
    Loads parameters specific to this module.
    :param root: an xml.etree.ElementTree Root object.
    """
    global scale_factor
    global bacteria_peristalsis_rate
    global nutrient_peristalsis_rate
    global init_Grf
    global init_Grm
    global init_Gpf
    global init_Gpm
    global init_Gmm
    global init_Gmf
    global inoculation_rate
    global gut_carrying_dry_weight
    global excretion_rate
    # read in parameters
    node = root.find('./mouse/scale_factor')
    if node is not None and node.text is not None:
        scale_factor = float(node.text) * 1000000.0
    node = root.find('./mouse/bacteria_peristalsis_rate')
    if node is not None and node.text is not None:
        bacteria_peristalsis_rate = float(node.text)
    node = root.find('./mouse/nutrient_peristalsis_rate')
    if node is not None and node.text is not None:
        nutrient_peristalsis_rate = float(node.text)
    node = root.find('./mouse/init_Grf')
    if node is not None and node.text is not None:
        init_Grf = int(node.text)
    node = root.find('./mouse/init_Grm')
    if node is not None and node.text is not None:
        init_Grm = int(node.text)
    node = root.find('./mouse/init_Gpf')
    if node is not None and node.text is not None:
        init_Gpf = int(node.text)
    node = root.find('./mouse/init_Gpm')
    if node is not None and node.text is not None:
        init_Gpm = int(node.text)
    node = root.find('./mouse/init_Gmm')
    if node is not None and node.text is not None:
        init_Gmm = int(node.text)
    node = root.find('./mouse/init_Gmf')
    if node is not None and node.text is not None:
        init_Gmf = int(node.text)
    node = root.find('./mouse/inoculation_rate')
    if node is not None and node.text is not None:
        inoculation_rate = float(node.text)
    node = root.find('./mouse/gut_carrying_weight')
    if node is not None and node.text is not None:
        gut_carrying_dry_weight = float(node.text)
    node = root.find('./mouse/excretion_rate')
    if node is not None and node.text is not None:
        excretion_rate = float(node.text)


def excess_digesta(digesta_weight):
    if digesta_weight > gut_carrying_dry_weight:
        return digesta_weight - gut_carrying_dry_weight
    return 0.0


def set_scale_factor(val):
    """ Use this to adjust the scale factor for simulating bacteria. """
    global scale_factor
    scale_factor = val


class Mouse(gutsim.vessel.Vessel):
    """
    Represents and simulates the colon of a single mouse. Maintains the bacterial cultures, nutrient pools that are
    replenished as the mouse eats a particular diet.

    :param nutrientInputs: a list of regimes by which nutrients are provided to the bacteria.
    """
    def __init__(self, nutrientInputs,
                daily_wheatstarch_kj=None, daily_dextrinised_kj=None, daily_sucrose_kj=None,
                daily_nitrogen_kj=None, daily_fat_kj=None, daily_chow_g=None, cageID=None, dietCode=None,
                *args, **kwargs):
        super(Mouse, self).__init__(*args, **kwargs)
        self.cageID = cageID
        self.dietCode = dietCode
        
        self.nutrientInputs = nutrientInputs
        self.daily_wheatstarch_kj 	= daily_wheatstarch_kj
        self.daily_dextrinised_kj 	= daily_dextrinised_kj
        self.daily_sucrose_kj 		= daily_sucrose_kj
        self.daily_nitrogen_kj 		= daily_nitrogen_kj
        self.daily_fat_kj			= daily_fat_kj
        self.daily_chow_g           = daily_chow_g
        self.logger = None
        # Used for recording how many bugs have been excreted. This is read and reset to zero externally. 
        # Don't use for anything critical within the simulation.
        self.bacteria_excreted_logging = 0

    def initialise_nutrient_Store(self, timestep, run_time):
        """
        Sets up the nutrient store along the GIT such that it appears to have been running for ages.

        :param run_time: given in hours.
        """
        self._nutrientStores = []
        currentTime = 0.0
        while currentTime < run_time:
            currentTime += timestep
            # Nutrients added
            nutStore = Nutrients()
            for ni in self.nutrientInputs:
                nutStore.add(ni.getNutrients(currentTime, timestep))
            self._nutrientStores.insert(0, nutStore)  # Prepend to start of gut. 


    def execute(self, endtime, timestep=0.1, verbose=False, bacteria_graphing_start_time=float('inf'), graphPath=None):
        """
        Called to execute this simulation of a mouse. Execution ceases after `endtime` hours.

        """
        def inoculum():
            """ Inputs new bacteria of each guild into the GIT. """
            rateDay = inoculation_rate * scale_factor
            rateHour = rateDay / 24.0
            # This gives a probability of insertion of a bug every time step.
            rateTimestep = rateHour * timestep
            new_bacteria = []  # Compile list of new bacteria, which is shuffled and then placed at start of bacteria
            # if rateTimestep > 1.0:
            #     raise Exception('Cannot insert more than 1 bug per timestep (with current code).')
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Grf != 0:
                new_bacteria.append(gutsim.bacteria.Guild_rf_f(timestep=timestep))
                bugs_to_add -= 1
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Grm != 0:
                new_bacteria.append(gutsim.bacteria.Guild_rm_m(timestep=timestep))
                bugs_to_add -= 1
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Gpf != 0:
                new_bacteria.append(gutsim.bacteria.Guild_pf_f(timestep=timestep))
                bugs_to_add -= 1
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Gpm != 0:
                new_bacteria.append(gutsim.bacteria.Guild_pm_m(timestep=timestep))
                bugs_to_add -= 1
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Gmm != 0:
                new_bacteria.append(gutsim.bacteria.Guild_m_m(timestep=timestep))
                bugs_to_add -= 1
            bugs_to_add = rateTimestep
            while random.random() <= bugs_to_add and init_Gmf != 0:
                new_bacteria.append(gutsim.bacteria.Guild_mf_f(timestep=timestep))
                bugs_to_add -= 1
            shuffle(new_bacteria)  # Shuffle the additional bugs.
            new_bacteria.extend(self._bacteria)  # Insert new bacteria at start.
            self._bacteria = new_bacteria


        def bacteria_peristalsis():
            """
            Represents peristaliss-driven re-ordering of bacteria in gut.  The distance that the shuffling attempts to
            move each bacteria by is drawn from a gaussian distribution. The spread of the distribution can be altered
            to be rate-limited by the timestep.

            Useful literature: Roberts RR, Murphy JF, Young HM, Bornstein JC. Development of colonic motility in the
            neonatal mouse-studies using spatiotemporal maps. Am J Physiol-Gastr L. 2007;292(3):G930-G8.
            """
            # Calculate the standard deviation of locations that bacteria may be shuffled by due to peristalsis events.
            # Peristalsis presents a probability distribution representing the post-peristalsis location each bacteria
            # could end up in. This is represented by a gaussian curve, centred on each bacterium's current location.
            # start by calculating how many peristalsis events happen per time step. Gaussian standard dev represents
            # effect of shuffling, so we convert from distance in gut to number of bacteria places (they are ordered).
            # Lastly, time is accounted for.
            distance_per_bacteria = gut_length / len(self._bacteria)
            sigma = bacteria_peristalsis_rate / distance_per_bacteria  # Num of bacteria representing the desired dist.
            diffusion_rate = peristalsis_rate * timestep  # How many peristalsis events occur in each time step.
            sigma *= diffusion_rate  # Adjust for number of peristalsis events in each time step.
            # Start the shuffling.
            shuffled = []  # Will contain tuples: (post-shuffle movement, bacteria object).
            for index, bug in enumerate(self._bacteria):
                pos = float(index) + random.gauss(0, sigma)  # Current position + some random quantity. 
                shuffled.append( (pos, bug) )
            # Sorts by the first item in the tuples that the list contains.
            sort = sorted(shuffled, key=lambda tmp: tmp[0])
            self._bacteria = [tup[1] for tup in sort]  # Pull out the bacteria objects.


        def nutrient_peristalsis():
            """
            Moves nutrients between nutrient stores; a low level of mixing resulting from peristalsis. 
            """
            # Gaussian diffusion represent diffusion as a result of peristalsis, IE, how much diffusion results from
            # each peristalsis event. This variable indicates how many peristalsis events are represented in the current
            # time step (time step may represent more than one peristalsis event).
            diffusion_rate = peristalsis_rate * timestep
            sigma_time_adjusted = nutrient_peristalsis_rate * diffusion_rate
            num_nutrient_stores = len(self._nutrientStores)
            length_per_nutrient_store = gut_length / num_nutrient_stores
            # Create a new store, which will be populated with nutrients.
            new_nutrient_store = [Nutrients() for _ in self._nutrientStores]

            # Location of this store along the gut. Conceptualized as follows. Nutrient sources ordered sequentially.
            # Their locations defined as cx. For integration, need to know the lower and upper boundaries, which
            # are given in factors of l (l=length of nutrient source in the gut, and can vary with number of
            # nutrient sources).
            #
            # 0l           1l           2l           3l           4l           5l
            # -------------------------------------------------------------------
            # |    c0      |    c1      |    c2      |    c3      |    c4      | ...
            # -------------------------------------------------------------------
            #
            # Normalized gaussian curve integrates to 1, hence its integral describes the diffusion of nutrients
            # into a specified length of gut.
            # Pre-compute the diffusion quantities by index distance. This is an efficiency saving. Since the
            # distribution changes only for the number of nutrient stores, which only changes between time steps,
            # the same distribution applies for every nutrient store in a timestep.
            # Need compute only one half of the distribution, as it is symmetrical.
            diffusion_by_index_difference = [None] * len(self._nutrientStores)
            # Declared here for efficiency. The lower bound for nutrient source n+1 = the supper bound for nutrient
            # source n (so don't have to re-compute it).
            sink_upper_distance = None
            for k in range(len(self._nutrientStores)):
                loc_source = 0.5 * length_per_nutrient_store
                if sink_upper_distance is None:  # No boundaries have been calculated for first nutrient source.
                    sink_lower_distance = (k * length_per_nutrient_store) - loc_source
                else:  # Adopt upper boundary of the previous nutrient source location.
                        sink_lower_distance = sink_upper_distance
                sink_upper_distance = ((k + 1) * length_per_nutrient_store) - loc_source
                integral_lower = stats.phi_gauss(sink_lower_distance, mu=0, sigma=sigma_time_adjusted)
                integral_upper = stats.phi_gauss(sink_upper_distance, mu=0, sigma=sigma_time_adjusted)
                proportion_diffused = integral_upper - integral_lower
                diffusion_by_index_difference[k] = proportion_diffused
            # Diffusion only performed over ranges where a substanital difference in concentration is incurred. 
            # This saves computational expense, and makes little difference to simulation dynamics. 
            # Here we ascertain the distance over which (meaningful) diffusion will be performed.
            negligible_diffusion_distance = 0
            negligible_diffusion_threshold = 0.05 / len(self._nutrientStores)
            for i, diff in enumerate(diffusion_by_index_difference):
                if diff < negligible_diffusion_threshold:   # Not worth performing diffusion beyond this spatial range. 
                    negligible_diffusion_distance = i
                    break

            # Scan through each existing nutrient store, and distribute nutrients
            for i, n_source in enumerate(self._nutrientStores):
                # Proportion of nutrient store contents not moved anywhere because Gaussian stretches from -Inf to Inf.
                # This will be retained in the current nutrient store location.
                remaining = 1.0
                # Start from difference of 1, 0 = the source, handled separately below.
                # Diffusion distribution is symmetrical, so calculate both sides of i (getting progressively further
                # away) simultaneously. Is more efficient, as not duplicating calculations.
                for index_difference in range(1, negligible_diffusion_distance):
                    proportion_diffused = diffusion_by_index_difference[index_difference]   # use lookup table.
                    move = n_source.proportion_query(proportion_diffused)
                    i_low = i - index_difference
                    if i_low >= 0:
                        new_nutrient_store[i - index_difference].add(move)
                        remaining -= proportion_diffused
                    i_high = i + index_difference
                    if i_high <= len(self._nutrientStores) - 1:
                        new_nutrient_store[i + index_difference].add(move)
                        remaining -= proportion_diffused
                # Any nutrient not already dealt with (Gaussian stretches from -Inf to Inf) is retained
                move = n_source.proportion_query(remaining)
                new_nutrient_store[i].add(move)
            self._nutrientStores = new_nutrient_store


        self.logger = gutsim.vessel.Logger(vessel=self, bacteria_graphing_start_time=bacteria_graphing_start_time, 
                                           graphPath=graphPath)
        self.initialise_nutrient_Store(timestep, 24.0)
        # Initialise bacteria in the simulation
        for _ in range(init_Grf):
            self._bacteria.append(gutsim.bacteria.Guild_rf_f(timestep=timestep))
        for _ in range(init_Grm):
            self._bacteria.append(gutsim.bacteria.Guild_rm_m(timestep=timestep))
        for _ in range(init_Gpf):
            self._bacteria.append(gutsim.bacteria.Guild_pf_f(timestep=timestep))
        for _ in range(init_Gpm):
            self._bacteria.append(gutsim.bacteria.Guild_pm_m(timestep=timestep))
        for _ in range(init_Gmm):
            self._bacteria.append(gutsim.bacteria.Guild_m_m(timestep=timestep))
        for _ in range(init_Gmf):
            self._bacteria.append(gutsim.bacteria.Guild_mf_f(timestep=timestep))
        shuffle(self._bacteria)	  # Ensure bacteria are well mixed prior to simulation.

        print('------- Nutrient inputs --------')
        for regime in self.nutrientInputs:
            print('hourly nutrient inputs:')
            print(regime.regime_string())
        print('--------------------------------')

        current_time = 0.0
        # Cumulative count of how much digesta and bacteria should have been excreted. This is stored as float, and
        # excretion takes place when threshold of a nutrient store of sufficient size (or a whole bacteria) is reached.
        digesta_to_excrete = 0.0
        bacteria_to_excrete_cumulative = 0.0
        self.logger.bacterGraphs = bacteria_graphing_start_time
        self.logger.log(self._bacteria, current_time)
        
        # Step through time
        while current_time < endtime:
            current_time += timestep
            # Nutrients added. New nutrients (from upper GI tract and mucin) are added to new object, and added to front
            # of nutrient store list. 
            # Note that mucin is then redistributed along the length of the gut.
            nutStore = Nutrients()
            for ni in self.nutrientInputs:
                if ni.is_active(current_time):
                    nutStore.add(ni.getNutrients(current_time, timestep))
            # Clone the nutrient record for logging.
            nutrient_input = nutStore.clone()
            # Put new nutrient record at start of GIT.
            self._nutrientStores.insert(0, nutStore)
            # Secrete mucus evenly along the entire GIT.
            Cm = nutStore.Cm
            Nm = nutStore.Nm
            nutStore.Cm = 0.0
            nutStore.Nm = 0.0
            Cmi = float(Cm) / len(self._nutrientStores)	 # Secrete mucus uniformly along GIT.
            Nmi = float(Nm) / len(self._nutrientStores)
            for ni in self._nutrientStores:
                ni.add_quantities(Cm=Cmi, Nm=Nmi)

            if verbose:
                totNutrients = Nutrients()
                for n in self._nutrientStores:
                    totNutrients.add(n)
                print('ct = ' + str(current_time) + ' ' + self.logger.string_state() + '\n\t' + str(totNutrients))

            inoculum()
            if len(self._bacteria) > 0:  # No point shuffling nothing.
                bacteria_peristalsis()
            if len(self._nutrientStores) > 0:
                nutrient_peristalsis()

            # Step each bacteria. Each bacteria is assigned a nutrient store from which to extract nutrients. 
            # There are likely multiple bacteria feeding from each store.
            bPerNS = float(len(self._bacteria)) / len(self._nutrientStores)
            # Build up a new list of bacteria, representing the GIT, as the current list is iterated.
            nextGIT = []
            # Scale nutrients to the desired level (affects total number of bacteria). 
            # This is reversed after bacteria have consumed and grown (see below)
            nutStore.scale_by_factor(scale_factor)
            for ns in self._nutrientStores:
                ns.scale_by_factor(scale_factor)

            for i, b in enumerate(self._bacteria):
                # Calculate which index in the nutrients store this bacteria corresponds to. 
                # Int used here always rounds down. 
                # Even if bPerNS is not a whole number, bacteria will broadly still access the correct nutrient store.
                iNS = int(i / bPerNS)
                # Returns nutrients consumed, and either a new bacteria object or None.
                consumed, child = b.step(self._nutrientStores[iNS])
                self._nutrientStores[iNS].subtract(consumed)
                non_absorbable_organic_matter = consumed.total_fermentable_quantity() * non_absorbable_organic_matter_generation_rate
                self._nutrientStores[iNS].add_quantities(other_organic_matter=non_absorbable_organic_matter)
                # New bacteria go in front of the current bacterium.
                if child is not None:
                    nextGIT.append(child)
                nextGIT.append(b)  # Add the current bacterium also.
            # Update for next time step.
            self._bacteria = nextGIT

            # Scale nutrient quantities back to grams.
            for ns in self._nutrientStores:
                ns.scale_by_factor(1.0 / scale_factor)

            # Perform mouse defecation.
            # Awake for first 12 hours of day (nighttime; mice are nocturnal).
            # Having to convert to ints here, because floats and modulo operator are incompatible.
            if int(current_time * 10) % 240 < 120:
                excrete_timestep = excretion_rate * timestep
            else:
                # During the light phase (less activity), excretion rate drops to half of active phase. 
                # This gives 2/3 vs 1/3 pattern of total excretion across these phases.
                excrete_timestep = excretion_rate * 0.5 * timestep
            digesta_weight = sum([ns.total_digesta_weight() for ns in self._nutrientStores])
            # If the GIT's carrying capacity is exceeded, increase the excretion rate.
            if excess_digesta(digesta_weight):
                excrete_timestep *= 5.0  # Value of 5 is arbitrarily chosen. 
            # Keep a cumulative count of how much digesta should have been excreted. When the quantity of nutrient
            # in the last nutrient store is greater than this cumulative count, delete the nutrient store.
            digesta_to_excrete += digesta_weight * excrete_timestep
            # Can't continue to delete if there's nothing left
            while len(self._nutrientStores) > 0 and \
                    self._nutrientStores[-1].total_digesta_weight() <= digesta_to_excrete:
                digesta_to_excrete -= self._nutrientStores[-1].total_digesta_weight()
                del(self._nutrientStores[-1])

            # Quantity to be excreted
            bacteria_to_excrete_cumulative += excrete_timestep * len(self._bacteria)
            bacteria_excreted = int(bacteria_to_excrete_cumulative)
            if bacteria_excreted > 0:  # Safety. list[-0:] = [] will delete the entire list.
                self.bacteria_excreted_logging += bacteria_excreted  # Used for logging.
                bacteria_to_excrete_cumulative -= float(bacteria_excreted)
                self._bacteria[-bacteria_excreted:] = []

            # Perform logging of bacterial numbers.
            self.logger.log(self._bacteria, current_time, nutrient_input=nutrient_input)

    def stringBacteriaState(self):
        for b in self._bacteria:
            print(str(b))


