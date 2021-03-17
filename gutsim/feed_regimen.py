"""
Created on 28/05/2014

@author: Mark Read
"""
import math
from gutsim.nutrients import Nutrients


def si_absorption_protein(quantity):
    """
    Returns the quantity of nutrient extracted by the small intestine. `Quantity` is daily consumption in grams.
    """
    # SI absoprtion rate varies with the quantity of nutrient eaten. For low quantities of protein, the SI absorbs high
    # proportions. When protein quantities increase, the SI absorbs less. An increasing exponential decay function is
    # used to model proportion absorbed by the SI given the quantity of nutrient consumed.
    #
    # Function looks like this: y = c(1-e^(-kx)), where y is the quantity absorbed by the SI, k is a scaling factor
    # and x is the quantity eaten. (y/x) gives the proportion absorbed  ( p(x) ).
    #
    # k can be calculated directly, if desired values for x, y and c are known. The following equation is used to 
    # calibrate k based on 'standard' diets consumed by ileostomates, translated to the scales of a mouse. 
    # k = -ln(1-(y/c))/x  			note that ln is changed to log base m, if some value m is used intead of e.
    #
    # There are further restrictions on this curve. It cannot breach the y side of the line y=x. To do so would mean
    # the SI has absorbed more nutrient than is eaten. The curve must also go through the point (x,y), and approach the
    # asymptote 'c'. Note that these restrictions are not universally possible, in particular if the point (x,y) lies
    # close to the y=x line, then the curve will breach y=x. I believe that equations of the form outlined above can
    # only 'bend' so far, and the bend must be a smooth one according to e^x.
    #
    # The equation is calibrated such that 12% of the protein consumed
    # in a "standard" diet (currently 14:57:29-4) is absorbed. The shape of the curve around this single point is
    # arbitrarily determined, as we have no data against which to calibrate.
    #
    c = 2.0	 # The asymptote
    m = math.e
    k = -math.log(1.0 - ((0.37 * 0.88) / c), m) / 0.37
    y = c * (1.0 - math.exp(-k * quantity))			# this is the absolute quantity absorbed, NOT the proportion.
    return y


def si_absorption_wheatstarch(quantity):
    """ Returns quantity of nutrient extracted by the small instestine. `Quantity` is daily consumption in grams. """
    # SI absoption rate varies with the quantity of nutrient eaten. For low quantities of carbs, the SI absorbs high
    # proportions. When carb quantities increase, the SI absorbs less. The form of the equation is described above for
    # protein. It is calibrated such that some % of the carbs consumed
    # in a "standard" diet (currently 14:57:29-4) are absorbed. The shape of the curve around this single point is
    # arbitrarily determined, as we have no data against which to calibrate.
    #
    # Reported figures range from 2-11%, mostly in humans. Potato starch is reported to be low, around 3%, however this
    # is cooked. Mouse feed probably is not cooked, hence these polysaccharides are less degraded. We assume that 30%
    # of wheatstarch is resistant starch.
    c = 2.0				# the asymptote
    m = math.e
    k = -math.log(1.0 - ((0.34 * 0.70) / c), m) / 0.34
    y = c * (1.0 - math.exp(-k * quantity))			# this is the absolute quantity absorbed, NOT the proportion.
    return y


def si_absorption_dextrinised(quantity):
    """	Returns quantity of nutrient extracted by the small instestine. `Quantity` is daily consumption in grams. """
    # SI absoprtion rate varies with the quantity of nutrient eaten. For low quantities of carbs, the SI absorbs high
    # proportions. When carb quantities increase, the SI absorbs less. The form of the equation is described above for
    # protein. It is calibrated such that some % of the carbs consumed
    # in a "standard" diet (currently 14:57:29-4) are absorbed. The shape of the curve around this single point is
    # arbitrarily determined, as we have no data against which to calibrate.
    #
    # Reported figures range from 2-11%, mostly in humans. Potato starch is reported to be low, around 3%, however this
    # is cooked. Mouse feed probably is not cooked, hence these polysaccharides are less degraded. We assume that
    # nearly all dextrinised can be absorbed through the small intestine, with 6% reaching the colon.
    c = 2.0				# the asymptote
    m = math.e
    k = -math.log(1.0 - ((0.24 * 0.94) / c), m) / 0.24
    y = c * (1.0 - math.exp(-k * quantity))			# this is the absolute quantity absorbed, NOT the proportion.
    return y


def si_absorption_sucrose(quantity):
    """ Returns quantity of nutrient extracted by the small intestine. `Quantity` is daily consumption in grams. """
    # assume that the SI absorbs ALL the sucrose.
    return quantity


class NutrientInput(object):
    """ Superclass from which all other nutritional sources or diet regimens inherit. """
    def __init__(self, regimen_start_h, regimen_end_h):
        self.daily_si_absorbed = None  # NutrientInput object
        self.daily_colic_input = 0  # NutrientInput object
        self.si_absorption_prop = 0  # NutrientInput object
        self.daily_intake = None   # Assigned by concrete class initialization.
        self._regimen_start_h = regimen_start_h
        self._regimen_end_h = regimen_end_h

    def is_active(self, time_h):
        """ Determines whether this regmine is active at the time specified. """
        return self._regimen_start_h <= time_h < self._regimen_end_h

    def calculate_SI_behaviour(self, dailyCw, dailyCd, dailyCs, dailyNf):
        """ Daily quantity of nutrients absorbed by the small intestine. """
        self.daily_si_absorbed = Nutrients(Cw=si_absorption_wheatstarch(dailyCw),
                                           Cd=si_absorption_dextrinised(dailyCd),
                                           Cs=si_absorption_sucrose(dailyCs),
                                           Nf=si_absorption_protein(dailyNf))
        self.daily_colic_input = Nutrients(Cw=dailyCw - self.daily_si_absorbed.Cw,
                                           Cd=dailyCd - self.daily_si_absorbed.Cd,
                                           Cs=dailyCs - self.daily_si_absorbed.Cs,
                                           Nf=dailyNf - self.daily_si_absorbed.Nf)
        # Conditional assignments to avoid div by zero
        prop_Cw = self.daily_si_absorbed.Cw/self.daily_intake.Cw if self.daily_intake.Cw != 0 else 0
        prop_Cd = self.daily_si_absorbed.Cd/self.daily_intake.Cd if self.daily_intake.Cd != 0 else 0
        prop_Cs = self.daily_si_absorbed.Cs/self.daily_intake.Cs if self.daily_intake.Cs != 0 else 0
        prop_Nf = self.daily_si_absorbed.Nf/self.daily_intake.Nf if self.daily_intake.Nf != 0 else 0
        self.si_absorption_prop = Nutrients(Cw=prop_Cw, Cd=prop_Cd, Cs=prop_Cs, Nf=prop_Nf)

    @staticmethod
    def _daytime_hourly(daily_quantity):
        """
        Hourly intake of nutrient during an hour in the day. Mice eat 1/3 of food during day.
        Jensen, T. L., Kiersgaard, M. K., Sorensen, D. B., & Mikkelsen, L. F. (2013).
        Fasting of mice: a review. Laboratory Animals, 47(4), 225-40.
        http://doi.org/10.1177/0023677213501659
        """
        return daily_quantity * (1.0 / 3.0) / 12.0

    @staticmethod
    def _nighttime_hourly(daily_quantity):
        """
        Hourly intake of nutrient during an hour of the night. Mice eat 2/3rds of food during night.
        Jensen, T. L., Kiersgaard, M. K., Sorensen, D. B., & Mikkelsen, L. F. (2013).
        Fasting of mice: a review. Laboratory Animals, 47(4), 225-40.
        http://doi.org/10.1177/0023677213501659
        """
        return daily_quantity * (2.0 / 3.0) / 12.0


class Constant(NutrientInput):
    """
    Represents a feeding regimen where mouse is fed at a constant rate over 24h.
    """
    def __init__(self, dailyCw, dailyCd, dailyCs, dailyNf, dailyCellulose, regimen_start_h=0., regimen_end_h=math.inf):
        """
        The daily* args specify the quantities of wheatstarch, dextrinised starch and protein that are consumed by the
        mouse each day.
        """
        super(Constant, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.daily_intake = Nutrients(Cw=dailyCw, Cd=dailyCd, Cs=dailyCs, Nf=dailyNf)
        self.dailyCw = dailyCw
        self.dailyCd = dailyCd
        self.dailyCs = dailyCs
        self.dailyNf = dailyNf
        self.dailyCellulose = dailyCellulose
        self.colicCw = dailyCw - si_absorption_wheatstarch(dailyCw)	 # daily nutrient quantities secreted into colon
        self.colicCd = dailyCd - si_absorption_dextrinised(dailyCd)
        self.colicCs = dailyCs - si_absorption_sucrose(dailyCs)
        self.colicNf = dailyNf - si_absorption_protein(dailyNf)
        self.colicCellulose = dailyCellulose
        # daily quantity of nutrients absorbed by the small intestine.
        self.calculate_SI_behaviour(dailyCw=dailyCw, dailyCd=dailyCd, dailyCs=dailyCs, dailyNf=dailyNf)

    def getNutrients(self, time, timestep):
        """
        Returns the nutrient input from the small intestine to the colin given the specified current time. The
        timestep is also supplied as the period of time that has passed. Note that the timestep should be given as a
        proportion of an hour, not a day.
        """
        time = float(time)
        day_cycle = time % 24.0
        if day_cycle < 12.0:  # Nighttime.
            inputCw = self._nighttime_hourly(self.colicCw)
            inputCd = self._nighttime_hourly(self.colicCd)
            inputNf = self._nighttime_hourly(self.colicNf)
            inputCellulose = self._nighttime_hourly(self.colicCellulose)
        else:  # Daytime
            inputCw = self._daytime_hourly(self.colicCw)
            inputCd = self._daytime_hourly(self.colicCd)
            inputNf = self._daytime_hourly(self.colicNf)
            inputCellulose = self._daytime_hourly(self.colicCellulose)
        # Convert from hourly rates to whatever the time-step is.
        inputCw *= timestep
        inputCd *= timestep
        inputNf *= timestep
        inputCellulose *= timestep
        return Nutrients(Cw=inputCw, Cd=inputCd, Nf=inputNf, non_ferm=inputCellulose)

    def regime_string(self):
        return ('\nAd lib feeding. Daily intakes (g/day): {:f} Cw; {:f} Cd; {:f} Cs; {:f} Nf. Starting at {:f}h, '
                'ending at {:f}. 66% consumption at night, remainder in day.') \
            .format(self.dailyCw, self.dailyCd, self.dailyCs, self.dailyNf, self._regimen_start_h, self._regimen_end_h)


class FiveTwo(NutrientInput):
    def __init__(self, Cw5, Cd5, Cs5, Nf5, cellulose5, fastProp, fast1, fast2, regimen_start_h=0., regimen_end_h=math.inf):
        """
        Arguments:
          Cw5, Cd5, Nf5 - the quantities of wheatstarch, dextrinised starch and protein consumed on a non-diet day.
          fast1 - which day of the week (numbered 1-7) constitutes the first fast day
          fast2 - same but for second fast day.
          fastProp - the proportion of a normal diet eaten on a fasting day (between 0 and 1).
        """
        super(FiveTwo, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.daily_intake = Nutrients(Cw=Cw5, Cd=Cd5, Cs=Cs5, Nf=Nf5)

        self._dailyCw5 = Cw5 - si_absorption_wheatstarch(Cw5)
        self._dailyCd5 = Cd5 - si_absorption_dextrinised(Cd5)
        self._dailyCs5 = Cs5 - si_absorption_sucrose(Cs5)
        self._dailyNf5 = Nf5 - si_absorption_protein(Nf5)
        self._dailyCellulose5 = cellulose5

        # Calculate diet day nutrient secretions into the colon.
        Cw2 = fastProp * Cw5
        Cd2 = fastProp * Cd5
        Nf2 = fastProp * Nf5
        cellulose2 = fastProp * cellulose5
        self._dailyCw2 = Cw2 - si_absorption_wheatstarch(Cw2)
        self._dailyCd2 = Cd2 - si_absorption_dextrinised(Cd2)
        self._dailyNf2 = Nf2 - si_absorption_protein(Nf2)
        self._dailyCellulose2 = cellulose2
        # Start time is inclusive, end time is exclusive. Ie, if end = 24h, then on hour 24 feeding stops.
        self._fast1Start = 24.0 * (fast1-1)	 # minus 1 here so first day of week is zero.
        self._fast1End = (24.0 * fast1)      # 1 hour before the start of the next day.
        self._fast2Start = 24.0 * (fast2-1)	 # minus 1 here so first day of week is zero.
        self._fast2End = (24.0 * fast2)  	 # 1 hour before the start of the next day.
        # daily quantity of nutrients absorbed by the small intestine.
        self.calculate_SI_behaviour(dailyCw=Cw5, dailyCd=Cd5, dailyCs=Cs5, dailyNf=Nf5)

    def _in_fast1(self, cycleTime):
        """ Returns true when current time falls on diet day 1. """
        return self._fast1Start <= cycleTime < self._fast1End

    def _in_fast2(self, cycleTime):
        """ Returns true when current time falls on diet day 2. """
        return self._fast2Start <= cycleTime < self._fast2End

    @staticmethod
    def _daytime_hourly(daily_quantity):
        """ Hourly intake of nutrient during an hour in the day. Mice eat 1/3 of food during day. """
        return daily_quantity * (1.0 / 3.0) / 12.0

    @staticmethod
    def _nighttime_hourly(daily_quantity):
        """ Hourly intake of nutrient during an hour of the night. Mice eat 2/3rds of food during night. """
        return daily_quantity * (2.0 / 3.0) / 12.0

    def getNutrients(self, time, timestep):
        """
        Returns the nutrient input from the small intestine to the colin given the specified current time. The
        timestep is also supplied as the period of time that has passed. Note that the timestep should be given as a
        proportion of an hour, not a day.
        """
        cycle_week = time % 168.0  # How far through the week the current time is, in hours.
        # Defaults, adjusted as needed. 
        dailyCw = self._dailyCw5  # Assume on a non-fast day to begin with.
        dailyCd = self._dailyCd5
        dailyNf = self._dailyNf5
        dailyCellulose = self._dailyCellulose5
        if self._in_fast1(cycle_week) or self._in_fast2(cycle_week):  # Check for 2-day diet.
            dailyCw = self._dailyCw2
            dailyCd = self._dailyCd2
            dailyNf = self._dailyNf2
            dailyCellulose = self._dailyCellulose2

        cycle24 = time % 24.0
        if cycle24 < 12.0:   # Nighttime
            hourlyCw = self._nighttime_hourly(dailyCw)
            hourlyCd = self._nighttime_hourly(dailyCd)
            hourlyNf = self._nighttime_hourly(dailyNf)
            hourlyCellulose = self._nighttime_hourly(dailyCellulose)
        else:  # Daytime
            hourlyCw = self._daytime_hourly(dailyCw)
            hourlyCd = self._daytime_hourly(dailyCd)
            hourlyNf = self._daytime_hourly(dailyNf)
            hourlyCellulose = self._daytime_hourly(dailyCellulose)
        inputCw = hourlyCw * timestep
        inputCd = hourlyCd * timestep
        inputNf = hourlyNf * timestep
        inputCellulose = hourlyCellulose * timestep
        return Nutrients(Cw=inputCw, Cd=inputCd, Nf=inputNf, non_ferm=inputCellulose)

    def regime_string(self):
        return  'Diet-derived nutrients entering the colon, following small intestinal absorption.\n' + \
                'These are daily figures:\n' + \
                'Ordinary day: Cw=' + str(self._dailyCw5) + ', Cd=' + str(self._dailyCd5) + ', Nf=' + str(self._dailyNf5) + '\n' + \
                'Fasting day: Cw=' + str(self._dailyCw2) + ', Cd=' + str(self._dailyCd2) + ', Nf=' + str(self._dailyNf2)


class TRF(NutrientInput):
    """
    Represents a feeding regimen where mouse is fed every 24 hours, but the duration of time for which nutrients pass
    from the small intestine into the colon van be varied.
    """
    def __init__(self, dailyCw, dailyCd, dailyCs, dailyNf, dailyCellulose, dailyCi, dailyCi_inverse, 
                 start_time, colicInputDuration, regimen_start_h=0., regimen_end_h=math.inf):
        """
        The daily* args specify the quantities of wheatstarch, dextrinised startch and protein that are consumed by the
        mouse each day. ColicInputDuration specifies the duration of time for which nutrients are passed from the small
        intestine into the colon.
        """
        super(TRF, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.daily_intake = Nutrients(Cw=dailyCw, Cd=dailyCd, Cs=dailyCs, Nf=dailyNf, Ci=dailyCi)
        self.colicInputDuration = colicInputDuration
        self.start_time = start_time
        self.end_time = (self.start_time + colicInputDuration) % 24.0  # modulo in case feeding period spans 24h cycle
        self.dailyCw = dailyCw
        self.dailyCd = dailyCd
        self.dailyCs = dailyCs
        self.dailyNf = dailyNf
        self.colicCw = dailyCw - si_absorption_wheatstarch(dailyCw)	 # Daily quantities of nutrients secreted into colon.
        self.colicCd = dailyCd - si_absorption_dextrinised(dailyCd)
        self.colicCs = dailyCs - si_absorption_sucrose(dailyCs)
        self.colicNf = dailyNf - si_absorption_protein(dailyNf)        
        self.colicCi = dailyCi  # Inulin is non-digestible
        self.colicCi_inverse = dailyCi_inverse - 0.0    # Inulin is non-digestible.
        self.colicCellulose = dailyCellulose
        self.hourlyCw = self.colicCw / self.colicInputDuration
        self.hourlyCd = self.colicCd / self.colicInputDuration
        self.hourlyNf = self.colicNf / self.colicInputDuration
        self.hourlyCi = self.colicCi / self.colicInputDuration
        self.hourlyCi_inverse = self.colicCi_inverse / (24.0 - self.colicInputDuration)
        self.hourlyCellulose = self.colicCellulose / self.colicInputDuration

        # Daily quantity of nutrients absorbed by the small intestine.
        self.calculate_SI_behaviour(dailyCw=dailyCw, dailyCd=dailyCd, dailyCs=dailyCs, dailyNf=dailyNf)

    def getNutrients(self, time, timestep):
        """
        Returns the nutrient input from the small intestine to the colin given the specified current time. The
        timestep is also supplied as the period of time that has passed. Note that the timestep should be given as a
        proportion of an hour, not a day.
        """
        feedingPeriod = 24.0
        time = float(time)
        cycleTime = time % 24.0	 # How far through the feeding cycle we are
        # Default cases if no input from small intestine
        inputCw = 0.0
        inputCd = 0.0
        inputNf = 0.0
        inputCi = 0.0
        inputCi_inverse = 0.0
        inputCellulose = 0.
        # Two cases. 1) start time is smaller than end time (easy).
        # 2) the feeding period spans the 24h modulo point, in which case the end time will be smaller than the start
        # time; feeding is defined by being greater than the start time or smaller than the end time.
        if self.start_time <= cycleTime < self.end_time or \
            self.end_time <= self.start_time and (self.start_time <= cycleTime or cycleTime < self.end_time):

            inputCw = self.hourlyCw * timestep
            inputCd = self.hourlyCd * timestep
            inputNf = self.hourlyNf * timestep
            inputCi = self.hourlyCi * timestep
            inputCellulose = self.hourlyCellulose * timestep
        else:  # Fasting time
            inputCi_inverse = self.hourlyCi_inverse * timestep

        inputCi += inputCi_inverse   # Combine at this point, and deliver to mouse.
        return Nutrients(Cw=inputCw, Cd=inputCd, Nf=inputNf, Ci=inputCi, non_ferm=inputCellulose)

    def regime_string(self):
        return '\nTRF. Daily intakes (g/day): {:f} Cw; {:f} Cd; {:f} Cs; {:f} Nf. Starting at {:f}h, ending at {:f}h' \
            .format(self.dailyCw, self.dailyCd, self.dailyCs, self.dailyNf, self._regimen_start_h, self._regimen_end_h)


class AlternateDay(NutrientInput):
    """
    2 day cycle: ad lib for a day, and then a day of fasting. The dietary restriction of fasting days can be specified.
    Inulin fibre can also be administered through this mechanism, in two ways. 1) inulin given with food. 2) inulin
    given when fasting. Quantity of inulin is either all or nothing, no partial restriction for fasting periods.
    """
    def __init__(self, dailyCw, dailyCd, dailyCs, dailyNf, dailyCellulose, dailyCi, dailyCi_inverse, fastProp,
                 regimen_start_h=0., regimen_end_h=math.inf):
        """
        The daily* args specify the quantities of wheat starch, dextrinised starch and protein that are consumed by the
        mouse each day. ColicInputDuration specifies the duration of time for which nutrients are passed from the small
        intestine into the colon.
        """
        super(AlternateDay, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.daily_intake = Nutrients(Cw=dailyCw, Cd=dailyCd, Cs=dailyCs, Nf=dailyNf, Ci=dailyCi, 
                                      non_ferm=dailyCellulose)

        self.colicCw = dailyCw - si_absorption_wheatstarch(dailyCw)	 # Daily quantities of nutrients secreted into colon.
        self.colicCd = dailyCd - si_absorption_dextrinised(dailyCd)
        self.colicCs = dailyCs - si_absorption_sucrose(dailyCs)
        self.colicNf = dailyNf - si_absorption_protein(dailyNf)
        self.colicCi = dailyCi  # Inulin is non-digestible
        self.colicCi_inverse = dailyCi_inverse  # Inulin is non-digestible.
        self.colicCellulose = dailyCellulose  # Non-digestible.
        self.fastProp = fastProp  # Alter the proportion of nutrients available between 0 and 1
        # daily quantity of nutrients absorbed by the small intestine.
        self.calculate_SI_behaviour(dailyCw=dailyCw, dailyCd=dailyCd, dailyCs=dailyCs, dailyNf=dailyNf)

    def getNutrients(self, time, timestep):
        """
        Returns the nutrient input from the small intestine to the colin given the specified current time. The
        timestep is also supplied as the period of time that has passed. Note that the timestep should be given as a
        proportion of an hour, not a day.
        """
        time = float(time)
        day_cycle = time % 24.0
        two_day_cycle = time % 48.0
        # First, check if day or night time.
        if day_cycle < 12.0:  # Nighttime.
            inputCw = self._nighttime_hourly(self.colicCw)
            inputCd = self._nighttime_hourly(self.colicCd)
            inputNf = self._nighttime_hourly(self.colicNf)
            inputCi = self._nighttime_hourly(self.colicCi)
            inputCi_inverse = self._nighttime_hourly(self.colicCi_inverse)
            inputCellulose = self._nighttime_hourly(self.colicCellulose)
        else:  # Daytime
            inputCw = self._daytime_hourly(self.colicCw)
            inputCd = self._daytime_hourly(self.colicCd)
            inputNf = self._daytime_hourly(self.colicNf)
            inputCi = self._daytime_hourly(self.colicCi)
            inputCi_inverse = self._daytime_hourly(self.colicCi_inverse)
            inputCellulose = self._daytime_hourly(self.colicCellulose)
        # Next, check if current time is fasting day
        if 24.0 <= two_day_cycle:
            inputCw *= self.fastProp
            inputCd *= self.fastProp
            inputNf *= self.fastProp
            inputCi = 0.0   # on fasting days, give no inulin.
            inputCellulose += self.fastProp
        else:  # Feeding day.
            inputCi_inverse = 0.0  # On feeding days, give no inulin.

        # Lastly, convert from hourly rates to whatever the time-step is.
        inputCw *= timestep
        inputCd *= timestep
        inputNf *= timestep
        inputCi *= timestep
        inputCi_inverse *= timestep
        inputCellulose *= timestep

        inputCi += inputCi_inverse   # Combine at this point, and deliver to mouse.
        return Nutrients(Cw=inputCw, Cd=inputCd, Nf=inputNf, Ci=inputCi, non_ferm=inputCellulose)

    def regime_string(self):
        return ('Feed quantities, daily: ' +
                'Cw=' + str(self.colicCw) + ' Cd=' + str(self.colicCd) + ' Nf=' + str(self.colicNf) +
                ' Ci=' + str(self.colicCi) + ' Cii=' + str(self.colicCi_inverse) + '\n' +
                'diet days: ' +
                'Cw=' + str(self.colicCw * self.fastProp) + ' Cd=' + str(self.colicCd * self.fastProp) +
                ' Nf=' + str(self.colicNf * self.fastProp)+ ' Ci=' + str(self.colicCi))


class InulinAdLib(NutrientInput):
    """
    Constant inulin supply to the mouse, modulated by mouse activity.
    """
    def __init__(self, dailyCi, regimen_start_h=0., regimen_end_h=math.inf):
        super(InulinAdLib, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.colicCi = dailyCi

    def getNutrients(self, time, timestep):
        time = float(time)
        day_cycle = time % 24.0
        if day_cycle < 12.0:    # nighttime, active phase
            inputCi = self._nighttime_hourly(self.colicCi)
        else:   # daytime, subdued phase
            inputCi = self._daytime_hourly(self.colicCi)
        inputCi *= timestep
        return Nutrients(Ci=inputCi)

    def regime_string(self):
        """ Summarises this regime's state as a string """
        # return '\ndaily inulin administration = {x}'.format(x=self.colicCi)
        return '\nInulin adlib. {:f} colic inulin g/day. Starting at {:f}h, ending at {:f}h'\
            .format(self.colicCi, self._regimen_start_h, self._regimen_end_h)


class InulinTimed(NutrientInput):
    """
    Time restricted feeding type of inulin administration. Not modulated by active (light) or dark phases.
    """
    def __init__(self, dailyCi, start_time, duration, regimen_start_h=0., regimen_end_h=math.inf):
        super(InulinTimed, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.duration = duration
        self.colicCi = dailyCi   # inulin is non-digestible by mouse.
        self.start_time = start_time
        self.end_time = (self.start_time + duration) % 24.0
        self.hourlyCi = self.colicCi / self.duration

    def getNutrients(self, time, timestep):
        """
        Returns the prebiotic input given the specified current time. The timestep is also supplied as the period of time that has passed.
        Note that the timestep should be given as a proportion of an hour, not a day.
        """
        time = float(time)
        cycleTime = time % 24.0	 # how far through the feeding cycle we are
        # default cases if no input from small intestine
        inputCi = 0.0
        # two cases. 1) start time is smaller than end time (easy).
        # 2) the feeding period spans the 24h modulo point, in which case the end time will be smaller than the start
        # time; feeding is defined by being greater than the start time or smaller than the end time.
        if self.start_time <= cycleTime < self.end_time or \
            self.end_time <= self.start_time and (self.start_time <= cycleTime or cycleTime < self.end_time):
            inputCi = self.hourlyCi * timestep
        return Nutrients(Ci=inputCi)

    def regime_string(self):
        """ Summarises this regimen's state as a string """
        return '\nTFR Inulin. {:f} Colic Ci g/day. Starting at {:f}h, ending at {:f}h' \
            .format(self.colicCi,self._regimen_start_h, self._regimen_end_h)


class Mucin(NutrientInput):
    """ Represents the secretion of mucin in the gut. """
    def __init__(self, dailyCm, dailyNm, regimen_start_h=0., regimen_end_h=math.inf):
        super(Mucin, self).__init__(regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h)
        self.daily_intake = Nutrients(Cm=dailyCm,  Nm=dailyNm)
        # convert daily mucus secretions into hourly quantities.
        self.dailyCm = dailyCm
        self.dailyNm = dailyNm
        self.hourlyCm = dailyCm / 24.0
        self.hourlyNm = dailyNm / 24.0
        self.calculate_SI_behaviour(dailyCw=0.0, dailyCd=0.0, dailyCs=0.0, dailyNf=0.0)
        self.daily_colic_input.Cm = dailyCm
        self.daily_colic_input.Nm = dailyNm

    def getNutrients(self, time, timestep):
        inputCm = self.hourlyCm * timestep
        inputNm = self.hourlyNm * timestep
        return Nutrients(Cm=inputCm, Nm=inputNm)  # nutrients not supplied default to zero.

    def regime_string(self):
        """ Summarises this regime's state as a string """
        return '\nMucin. {:f} Colic Cm g/day; {:f} Nm g/day. Starting at {:f}h, ending at {:f}h' \
            .format(self.dailyCm, self.dailyNm, self._regimen_start_h, self._regimen_end_h)


