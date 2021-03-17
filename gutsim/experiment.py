"""
Created on 15/04/2014

@author: markread
"""
import pandas.io.excel
import xml.etree.ElementTree as ET
import os
import sys  
import random
import shutil
# project-originating imports
import gutsim.feed_regimen
from gutsim.mouse import Mouse
import gutsim.bacteria
from gutsim.mouse import set_scale_factor
from gutsim.nutrients import Nutrients
import gutsim.nutrients as nutrients


def profile():
    """ Used to profile the code's execution. """
    sys.argv = "launch.py -m 0 -d results/24-3wkNewGut -dc 24-24 -end 80.0 -seed 0 -f -s 50 -lei -v ".split()
    main()


def calculate_daily_mucin_secretions(accWater=True):
    """
    Method encapsulates the calculations for mucin secretion, in terms of both nigrogen and carbon. Secretion figures
    are given as daily quantities (24 hours). Ouptut is in grams, assuming volume to weight ratio is similar to water.
    """
    # This paper suggests mucin is 80% carb.
    # Bansil, R., and Turner, B.S. (2006) Mucin structure, aggregation, physiological
    # functions and biomedical applications.
    mucin_carb_prop = 0.8  # Proportion of mucin that is carb
    mucin_prot_prop = 1.0 - mucin_carb_prop
    # This is in microns / hour. Examining a cross section of the epithelium, this is the height of the mucin
    # layer an hour after it has been removed.
    # Ermund etal, studies of mucus in mouse stomach.... pt 1. 2013.
    mucin_secretion_rate = 82  # in microns/hour
    # large and small intestine areas taken from Wolczuk - Morphometric Characteristics..., 2011
    li_area = 11.787  # cm^2
    si_area = 31.975  # cm^2. Not currently used.
    # Work in cm^2 and cm^3. Convert mucin secretion rate to cm/hour; 24 hours in day, and multiply by
    # area of secreting epithelial layer. This is in grams, assuming 1cm^3=1g.
    daily_mucus = (li_area) * 0.0001 * mucin_secretion_rate * 24
    # 90-99% of mucus is water, only the rest is composed of glycoprotein. See reference:
    # Macfarlane, S et al. Colonization of Mucin by Human Intestinal Bacteria and Establishment of Biofilm
    # Communities in a Two-Stage Continuous Culture System
    daily_mucins = daily_mucus
    if accWater:
        daily_mucins = daily_mucus * 0.1

    dailyNm = (daily_mucins * mucin_prot_prop)
    dailyCm = daily_mucins * mucin_carb_prop
    print('daily mucin glycan = {}'.format(dailyCm))
    print('daily mucin polypeptide = {}'.format(dailyNm))
    return dailyCm, dailyNm


class Diet:
    """
    Stores information pertaining to a particular diet that a mouse is fed (independently of how much the mouse
    actually ate, that is accounted for elsewhere).
    """
    def __init__(self, code, nitrogen_g, sucrose_g, wheatstarch_g, dextrinised_g, fat_g, cellulose_g,
                nitrogen_kj, sucrose_kj, wheatstarch_kj, dextrinised_kj, fat_kj):
        self.code = code
        self.nitrogen_g = nitrogen_g
        self.sucrose_g = sucrose_g
        self.wheatstarch_g = wheatstarch_g
        self.dextrinised_g = dextrinised_g
        self.fat_g = fat_g
        self.cellulose_g = cellulose_g
        self.nitrogen_kj = nitrogen_kj
        self.wheatstarch_kj = wheatstarch_kj
        self.dextrinised_kj = dextrinised_kj
        self.sucrose_kj = sucrose_kj
        self.fat_kj = fat_kj


def prepare_experiment(mucinWater, dietCode='24:3', **kwargs):
    """
    Returns a list of Mouse objects (`mice`), which are ready to be simulated.
    Parameters:
      mucinWater - when set to true, mucin nutrient content adjusted to account for water.
      dietCode - string that represents the feeding regimen.
    """
    # Read in the diet compositions
    diets_raw = pandas.read_excel('diet-nutrient-breakdown-SSB.xlsx', 'diet-compositions', header=1, index_col=0)
    codes = diets_raw.index.dropna()  # Diet codes have been adopted as the index. 
    diets = {}
    for code in codes:
        r = diets_raw.loc[code]  # This looks up a row by the first column item
        diets[code] = Diet(
            code=code,
            # These three values are grams per kilogram of chow.
            nitrogen_g=r['Nigrogen g/Kg'],
            sucrose_g=r['Sucrose g/Kg'],
            wheatstarch_g=r['Wheatstarch g/Kg'],
            dextrinised_g=r['Dextrinised Starch g/Kg'],
            fat_g=r['Soyabean oil g/Kg'],
            cellulose_g=r['Cellulose g/Kg'],
            # These values are KJ per kilogram of chow
            nitrogen_kj=r['Nigrogen KJ/Kg'],
            sucrose_kj=r['Carbon, sucrose KJ/Kg'],
            wheatstarch_kj=r['Carbon, wheatstarch KJ/Kg'],
            dextrinised_kj=r['Carbon, dextrinised starch KJ/Kg'],
            fat_kj=r['Fat KJ/Kg']
        )
    daily_mucinC_g, daily_mucinN_g = calculate_daily_mucin_secretions(mucinWater)

    miceRAW = pandas.read_excel('diet-nutrient-breakdown-SSB.xlsx', 'mice', index_col='Mouse No')
    mice = {}
    print(miceRAW.index)
    # Read in the mouse food consumption (in weight), diet, and other mouse related metadata.
    # for m in range(miceRAW.shape[0]):   # number of rows
    for mouse_num, row in miceRAW.iterrows():        
        cageID = row.loc['Cage No']
        code = row.loc['Diet-code']
        consumed_g = row.loc['Dry weight food eaten g/mouse/cage/d']  # chow consumed by mouse, in grams.

        # Calculate how many grams and KJ of each macronutrient were consumed by each mouse, based on diet composition
        # and quantity of chow consumed. Note the division by 1000, this is because the diet compositions are expressed
        # as grams and KJ per Kg chow, but the mouse consumption figures are given in grams. Hence the conversion of
        # diet compositions from grams or KJ per Kg chow to grams or KJ per GRAMS chow.
        daily_nitrogen_g   	 = diets[code].nitrogen_g * consumed_g / 1000
        daily_sucrose_g      = diets[code].sucrose_g * consumed_g / 1000
        daily_wheatstarch_g  = diets[code].wheatstarch_g * consumed_g / 1000
        daily_dextrinised_g  = diets[code].dextrinised_g * consumed_g / 1000
        daily_fat_g          = diets[code].fat_g * consumed_g / 1000
        daily_cellulose_g    = diets[code].cellulose_g * consumed_g / 1000

        daily_nitrogen_kj    = diets[code].nitrogen_kj * consumed_g / 1000
        daily_sucrose_kj     = diets[code].sucrose_kj * consumed_g / 1000
        daily_wheatstarch_kj = diets[code].wheatstarch_kj * consumed_g / 1000
        daily_dextrinised_kj = diets[code].dextrinised_kj * consumed_g / 1000
        daily_fat_kj         = diets[code].fat_kj * consumed_g / 1000

        nutrientInputs = []     # list of feeding regimes to be supplied to mouse.
        # constant input over 24h
        if dietCode == 'const':
            nutrientInputs.append(gutsim.feed_regimen.Constant(
                                    dailyCw=daily_wheatstarch_g, dailyCd=daily_dextrinised_g, dailyCs=daily_sucrose_g,
                                    dailyNf=daily_nitrogen_g, dailyCellulose=daily_cellulose_g))
        # time restricted feeding.
        if dietCode == 'trf':
            nutrientInputs.append(gutsim.feed_regimen.TRF(
                                  start_time=kwargs['trf-start'], colicInputDuration=kwargs['trf-duration'],
                                  dailyCw=daily_wheatstarch_g, dailyCd=daily_dextrinised_g, dailyCs=daily_sucrose_g,
                                  dailyNf=daily_nitrogen_g,
                                  dailyCellulose=daily_cellulose_g,
                                  dailyCi=kwargs['inulin_chow'], dailyCi_inverse=kwargs['inulin_chow_inverse']))
        if dietCode == '5-2':
            nutrientInputs.append(gutsim.feed_regimen.FiveTwo(
                                    Cw5=daily_wheatstarch_g, Cd5=daily_dextrinised_g, Cs5=daily_sucrose_g,
                                    Nf5=daily_nitrogen_g,
                                    cellulose5=daily_cellulose_g,
                                    fastProp=kwargs['ftdq'], fast1=kwargs['ft1'], fast2=kwargs['ft2']))

        if dietCode == '48-24':
            nutrientInputs.append(gutsim.feed_regimen.AlternateDay(
                                    dailyCw=daily_wheatstarch_g, dailyCd=daily_dextrinised_g, dailyCs=daily_sucrose_g,
                                    dailyNf=daily_nitrogen_g,
                                    dailyCellulose=daily_cellulose_g,
                                    fastProp=kwargs['addq'],
                                    dailyCi=kwargs['inulin_chow'], dailyCi_inverse=kwargs['inulin_chow_inverse']))
        if kwargs['inulin_adlib'] != 0.0:
            nutrientInputs.append(gutsim.feed_regimen.InulinAdLib(dailyCi=kwargs['inulin_adlib']))
        if kwargs['inulin_timed_dose'] != 0.0:
            nutrientInputs.append(gutsim.feed_regimen.InulinTimed(dailyCi=kwargs['inulin_timed_dose'],
                            start_time=kwargs['inulin_timed_start'], duration=kwargs['inulin_timed_duration']))

        # Add mucin secretions to the input.
        nutrientInputs.append(gutsim.feed_regimen.Mucin(dailyCm=daily_mucinC_g, dailyNm=daily_mucinN_g))
        mouse = Mouse(nutrientInputs=nutrientInputs, dietCode=code,
                    daily_wheatstarch_kj=daily_wheatstarch_kj,
                    daily_dextrinised_kj=daily_dextrinised_kj, daily_sucrose_kj=daily_sucrose_kj,
                    daily_nitrogen_kj=daily_nitrogen_kj,
                    daily_fat_kj=daily_fat_kj,
                    daily_chow_g=consumed_g,
                    ID=mouse_num, cageID=cageID)
        mice[mouse_num] = mouse
    return mice


def diet_analysis(mice):
    print('\nPERFORMING DIET ANALYSIS\n')
    diet_code_file = open('diet-analysis/diet-composition-per-mouse.csv','w')
    relNutFile = open('diet-analysis/nutrient_rel_abundances.csv','w')
    relNutFile.write('#{:8s},{:16s},{:16s},{:16s},'.format('mouse-ID', 'carb-intake-KJ/d', 'prot-intake-KJ/d', 'fat-intake-KJ/d')
                    + '{:7s},{:7s},{:7s},{:7s},{:7s},'.format('relCw', 'relCd', 'relCm', 'relNf', 'relNm')
                    + '{:7s},{:7s},{:7s},{:7s},{:7s},'.format('inCw','inCd','inCm','inNf','inNm')
                    + '{:7s},{:7s},{:7s},{:7s},{:7s},'.format('siCw','siCd','siCm','siNf','siNm')
                    + '{:7s},{:7s},{:7s},{:7s},{:7s},'.format('siPCw','siPCd','siPCm','siPNf','siPNm')
                    + '{:7s},{:7s},{:7s},{:7s},{:7s},'.format('colicCw', 'colicCd', 'colicCm', 'colicNf', 'colicNm')
                    + '{:7s},'.format('Nf:Nm')
                    + '{:14s},{:14s},{:14s},{:14},{:14},{:14s},'.format('cecum_relCw', 'cecum_relCd', 'cecum_relCmuc', 'cecum_relCcas', 'cecum_relNmuc', 'cecum_relNcas')
                    + '\n')

    diet_code_file.write('%5s, %10s, %6s, %6s, %6s, %6s, %6s\n' %
                        ('#ID', 'diet code', 'Nf_g', 'Cw_g', 'Cd_g', 'Nm_g', 'Cm_g'))
    for i, mouse_num in enumerate(mice.keys()):
        mouse = mice[mouse_num]
        # This code is used for writing nutrient relative abundances to the file system.
        ci = mouse.daily_wheatstarch_kj + mouse.daily_dextrinised_kj + mouse.daily_sucrose_kj
        pi = mouse.daily_nitrogen_kj
        fi = mouse.daily_fat_kj  # Fat intake

        # Instantiate variables, values cumulatively added based on each nutrient input regimen.
        intake = Nutrients()
        colic = Nutrients()
        si_absorb = Nutrients()
        si_prop = Nutrients()
        for nIn in mouse.nutrientInputs:
            intake.add(nIn.dailyIntake)	 # A Nutrients object representing daily generation of nutrients.
            colic.add(nIn.daily_colic_input)
            si_absorb.add(nIn.daily_si_absorbed)
            si_prop.add(nIn.si_absorption_prop)

        # These are by weight of macronutrient's carbon/nitrogen atoms.
        totC_weight = intake.carbon_content()
        totN_weight = intake.nitrogen_content()

        # These are percentages by weight of macronutrient carbon content.
        rel_carbon_Cw = 100 * nutrients.carbon_content_carbohydrate(intake.Cw) / totC_weight
        rel_carbon_Cd = 100 * nutrients.carbon_content_carbohydrate(intake.Cd) / totC_weight
        rel_carbon_Cs = 100 * nutrients.carbon_content_carbohydrate(intake.Cs) / totC_weight
        rel_carbon_Cm = 100 * nutrients.carbon_content_carbohydrate(intake.Cm) / totC_weight
        rel_nitrogen_Nf = 100 * nutrients.nitrogen_content_protein(intake.Nf) / totN_weight
        rel_nitrogen_Nm = 100 * nutrients.nitrogen_content_protein(intake.Nm) / totN_weight
        ratio_nitrogen_NFNm = intake.Nf / intake.Nm

        # These are percentages by weight of macronutrient carbon content.
        # Considers macronutrients available to bacteria, i.e. feed macronutrients post-small intestine.
        cecum_totC_weight = colic.carbon_content()  # These are by weight of macronutrient's carbon/nitrogen atoms.
        cecum_totN_weight = colic.nitrogen_content()

        print('cecum C weight ' + str(cecum_totC_weight) + ' ; total weight = ' + str(totC_weight))
        cecum_rel_carbon_Cw = 100 * nutrients.carbon_content_carbohydrate(colic.Cw) / cecum_totC_weight
        cecum_rel_carbon_Cd = 100 * nutrients.carbon_content_carbohydrate(colic.Cd) / cecum_totC_weight
        cecum_rel_carbon_Cs = 100 * nutrients.carbon_content_carbohydrate(colic.Cs) / cecum_totC_weight
        cecum_rel_carbon_mucin = 100 * \
                                 (nutrients.carbon_content_carbohydrate(colic.Cm) + nutrients.carbon_content_protein(colic.Nm)) \
                                 / cecum_totC_weight
        # Some carbon comes from casein too.
        cecum_rel_carbon_cas = 100 * nutrients.carbon_content_protein(colic.Nf) / cecum_totC_weight
        
        cecum_rel_nitrogen_mucin = 100 * nutrients.nitrogen_content_protein(colic.Nm) / cecum_totN_weight
        cecum_rel_nitrogen_cas = 100 * nutrients.nitrogen_content_protein(colic.Nf) / cecum_totN_weight

        relNutFile.write(
            '{:9d},{:10.10f},{:10.10f},{:10.10f},'.format(i, ci, pi, fi)  # Mouse index, then carb, prot and fat intakes
            # Relative abundances of carbon and nitrogen, by carbon atom weight of INTAKE/secreted quantities.
            + '{:.17f},{:.17f},{:.17f},{:.17f},{:.17f},'.format(rel_carbon_Cw, rel_carbon_Cd, rel_carbon_Cm, rel_nitrogen_Nf, rel_nitrogen_Nm)
            # Nutrients eaten
            + '{:.17f},{:.17f},{:.17f},{:.17f},{:.17f},'.format(intake.Cw, intake.Cd, intake.Cm, intake.Nf, intake.Nm)
            # Absolute quantity absorbed by the SI
            + '{:.17f},{:.17f},{:.17f},{:.17f},{:.17f},'.format(si_absorb.Cw, si_absorb.Cd, si_absorb.Cm, si_absorb.Nf, si_absorb.Nm)
            # Proportion of eaten nutrients absorbed by the SI
            + '{:.17f},{:.17f},{:.17f},{:.17f},{:.17f},'.format(si_prop.Cw, si_prop.Cd, si_prop.Cm, si_prop.Nf, si_prop.Nm)
            # Quantity reaching the colon.
            + '{:.17f},{:.17f},{:.17f},{:.17f},{:.17f},'.format(colic.Cw, colic.Cd, colic.Cm, colic.Nf, colic.Nm)
            + '{:.17f},'.format(ratio_nitrogen_NFNm)
            # Relative abundances of carbon and nitrogen (by atom weight, not total macronutrient weight), given
            # availability to bacteria. i.e., post-small intestine.
            + '{:14.2f},{:14.2f},{:14.2f},{:14.2f},{:14.2f},{:14.2f}'.format(cecum_rel_carbon_Cw, cecum_rel_carbon_Cd, cecum_rel_carbon_mucin, cecum_rel_carbon_cas, cecum_rel_nitrogen_mucin, cecum_rel_nitrogen_cas)
            + '\n')
        diet_code_file.write('%5d, %10s, %6f, %6f, %6f, %6f, %6f\n' %
                        (i, mouse.dietCode, intake.Nf, intake.Cw, intake.Cd, intake.Nm, intake.Cm) )

    diet_code_file.close()
    relNutFile.close()


def write_params_fs(filename):
    """ Writes the parameters supplied in `**kwargs` to the filesystem, in the supplied filename. """
    params = open(filename, 'w')
    for i in sys.argv:
        params.write(i + ' ')
    params.close()


def main():
    """
    The main entry point into the program.
    The following command line arguments can be provided:
      -m xxx: mouse number xxx will be executed. If absent, all mice will be executed.
      -d yyy: data will be stored under directory 'yyy', which will be created if it does not already exist.
      -anu: when specified bacteria can internalise alternative (non-preferred) nutrient sources.
      -dc xxx: diet code, specifies feeding regimen. 24-24 indicates ab libitum, 24-3 is daily where nutrient passes
        to the colon only for 3 hours out of 24. 5-2 represents five-two diet.
    """
    print('Argument List:', str(sys.argv))
    parameters_file = None              # location of parameters file, optional.
    timestep = 0.1                      # duration of time represented in each iteration of simulation, in hours.
    mouse = None						# mice are numbered sequentially from 0 to 249
    cageID = None						# these are the cage IDs from the dataset.
    seed = None
    seedSpecified = False               # if user provided seed, set to true
    subMouse = ''						# for multiple runs of the same mouse. eg, file names "mouse_0.XXX.png"
    exp_dir = None
    tag = ''							# text appended to end of standard generated data dirs. mice-data-24-3_TAG
    verbose = False						# when true, write frequenty simulation status to standard out.
    mucinWater = True					# when true, mucin nutrients are adjusted for water content.
    diet_code = None						# string, representing the feeding regimen.
    # If end time is specified, then this is duration of final hours of simulated time for which bacteria graphs
    # are to be drawn.  
    bacteria_graphs = None  # Float, time in hours prior to termination, at which point to start logging and plotting. 
    mouseGraphs = False					# draw graphs of mouse (and guild) level status.
    forceExecute = False
    ft1 = None							# first day of diet in 5-2 diet. Days of week numbered 1 to 7
    ft2 = None							# second day of diet in 5-2 diet. Days of week numbered 1 to 7
    ftdq = None							# proportion of full day's nutrient that diet days represent. 0 to 1.
    addq = None                         # proportion of full day's nutrient that diet days represent. 0 to 1, alternate day
    end_time = 168.0						# simulation end time, run for a week by default.
    dietAnalysis = False
    # specific to fibre that only Bdm group can access.
    fdmStart = None						# how far through the 24 hour cycle fibre should be input to colon
    fdmDur = None						# the duration of the 24 hour cycle that fibre is provided for
    fdmQuant = None						# the quantity of daily feed that fibre constitutes IN ADDITION to feed (grams).
    # read in command line args.
    i = 1   # index into args. Ignore the first one, since it's the name of the python module being run.
    while i < len(sys.argv):
        if sys.argv[i] == '-p':
            i += 1
            parameters_file = sys.argv[i]
        elif sys.argv[i] == '-m':		# mouse to be run.
            i += 1
            mouse = int(sys.argv[i])
        elif sys.argv[i] == '-c':		# select a particular cage number to execute (not the same as mouseID)
            i += 1
            cageID = sys.argv[i]
        elif sys.argv[i] == '-sm':		# for multiple runs of the same mouse. eg, file names "mouse_0.XXX.png"
            i += 1						# the '.' is added here, does not need to be specified.
            subMouse = '.' + sys.argv[i]
        elif sys.argv[i] == '-d':		# the experimental directory to which data should be written.
            i += 1
            exp_dir = sys.argv[i]
        elif sys.argv[i] == '-seed':    # set the random number seed.
            i += 1
            seed = int(sys.argv[i])
            seedSpecified = True
        elif sys.argv[i] == '-tag':		# text to append to the end of standard generated data folders.
            i += 1
            tag = '_' + sys.argv[i]
        elif sys.argv[i] == '-dc':		# diet code, which specifies the feeding regimen.
            i += 1
            diet_code = sys.argv[i]
        elif sys.argv[i] == '-ft1':		# for five-two diet, the first diet day. Days numbered 1-7. -1 for no diet.
            i += 1
            ft1 = int(sys.argv[i])
        elif sys.argv[i] == '-ft2':		# for five-two diet, the second diet day. Days numbered 1-7. -1 for no diet.
            i += 1
            ft2 = int(sys.argv[i])
        elif sys.argv[i] == '-ftdq':    # proportion of full day nutrients that diet days represent. 0..1.
            i += 1
            ftdq = float(sys.argv[i])
        elif sys.argv[i] == '-addq':    # proportion of full day nutrients that diet days represent. 0..1, alternate day
            i += 1
            addq = float(sys.argv[i])
        elif sys.argv[i] == '-v':		# verbose output to standard out.
            verbose = True
        elif sys.argv[i] == '-bg':		# bacteria graphs, whether to graph bacterial population dynamics
            # If end time is specified, then this is duration of final hours of simulated time for which bacteria graphs
            # are to be drawn. 
            i += 1
            bacteria_graphs = float(sys.argv[i]) 
        elif sys.argv[i] == '-mg':
            mouseGraphs = True
        elif sys.argv[i] == '-end':		# simulation termination time.
            i += 1
            end_time = float(sys.argv[i])
        elif sys.argv[i] == '-f':		# force the simulation to execute, even if data file already exists.
            forceExecute = True
        elif sys.argv[i] == '-da':		# diet analysis, creates files relating to nutrients available in colon.
            dietAnalysis = True
        elif sys.argv[i] == '-s':		# the scale factor, higher numbers lead to larger bacterial populations.
            i += 1						# supplied value is multiplied by a million.
            set_scale_factor(float(sys.argv[i]) * 1000000.0)
        i += 1

    kwargs = dict()
    kwargs['inulin_adlib'] = 0.0  # ad lib inulin
    kwargs['inulin_chow'] = 0.0  # default quantity of inulin supplement.
    kwargs['inulin_chow_inverse'] = 0.0     # default quantity of inulin supplement.
    kwargs['inulin_water'] = 0.0            # default quantity of inulin supplement.
    if diet_code == '5-2':
        if '-ft1' in sys.argv:
            i = sys.argv.index('-ft1')
            kwargs['ft1'] = int(sys.argv[i+1])
        if '-ft2' in sys.argv:
            i = sys.argv.index('-ft2')
            kwargs['ft2'] = int(sys.argv[i+1])
        if '-ftdq' in sys.argv:
            i = sys.argv.index('-ftdq')
            kwargs['ftdq'] = float(sys.argv[i+1])
    # parameters relating to time-restricted feeding.
    if diet_code == 'trf':
        i = sys.argv.index('-trd')
        avail = float(sys.argv[i + 1])
        kwargs['trf-duration'] = avail
        tag = '-' + str(avail) + tag
        print('time restricted feeding selected, feed available for ' + str(avail) + ' hours a day.')
    if diet_code == '48-24':
        if '-addq' in sys.argv:
            i = sys.argv.index('-addq')
            kwargs['addq'] = float(sys.argv[i+1])

    # reading parameters from xml file.
    if parameters_file is not None:
        tree = ET.parse(parameters_file)
        root = tree.getroot()
        # load parameters pertaining to experiment here.
        timestep = float(root.find('./timestep').text)
        node = root.find('./exp_dir')
        if node is not None and node.text is not None:
            exp_dir = node.text
        node =root.find('./end')
        if node is not None and node.text is not None:
            end_time = float(node.text)
        node = root.find('./tag')
        if node is not None and node.text is not None:
            tag = node.text
        node = root.find('./verbose')
        if node is not None and node.text is not None:
            verbose = bool(int(node.text))   # '0' or '1'
        node = root.find('./force')
        if node is not None and node.text is not None:
            forceExecute = bool(int(node.text))
        node = root.find('./mouse')
        if node is not None and node.text is not None:
            # if nothing specified, defer to command line argument.
            if node.text.lower() == 'all':
                mouse = None   # ensures all mice are executed.
            else:
                mouse = int(node.text)
        node = root.find('./output/bacteria_graphs')
        if node is not None and node.text is not None:
            bacteria_graphs = float(node.text)
        node = root.find('./output/mouse_graphs')
        if node is not None and node.text is not None:
            mouseGraphs = bool(int(node.text))
            print("mouse graphs = " + str(mouse))
        # load parameters pertaining to diet code.
        node = root.find('./diet_code')
        if node is not None and node.text is not None:
            diet_code = node.text
        if diet_code == 'trf':
            node = root.find('./diet_regimes/time_restricted/start')
            if node is not None and node.text is not None:
                kwargs['trf-start'] = float(node.text)
                tag = '-s' + node.text + tag
            node = root.find('./diet_regimes/time_restricted/duration')
            if node is not None and node.text is not None:
                kwargs['trf-duration'] = float(node.text)
                tag = '-d' + node.text + tag
        if diet_code == '5-2':
            # days are numbered 1-7. -1 means no diet day.
            node = root.find('./diet_regimes/five_two/day1')
            if node is not None and node.text is not None:
                kwargs['ft1'] = int(node.text)
            node = root.find('./diet_regimes/five_two/day2')
            if node is not None and node.text is not None:
                kwargs['ft2'] = int(node.text)
            node = root.find('./diet_regimes/five_two/diet_quantity')
            if node is not None and node.text is not None:
                kwargs['ftdq'] = float(node.text)
        if diet_code == '48-24':
            node = root.find('./diet_regimes/alternate_day/diet_quantity')
            if node is not None and node.text is not None:
                kwargs['addq'] = float(node.text)
        # inulin related parameters
        node = root.find('./diet_regimes/inulin/adlib')
        if node is not None and node.text is not None:
            kwargs['inulin_adlib'] = float(node.text)
        node = root.find('./diet_regimes/inulin/timed')
        if node is not None:
            node = root.find('./diet_regimes/inulin/timed/dose')
            kwargs['inulin_timed_dose'] = float(node.text)
            node = root.find('./diet_regimes/inulin/timed/start')
            kwargs['inulin_timed_start'] = float(node.text)
            node = root.find('./diet_regimes/inulin/timed/duration')
            kwargs['inulin_timed_duration'] = float(node.text)

        # Load parameters for other parts of the simulation. Each module handles it's own parameters.
        gutsim.mouse.load_parameters(root)
        gutsim.bacteria.load_parameters(root)
        print('parameters file specified: ' + parameters_file)
        print('experimental directory: ' + exp_dir)
    # Create experimental directory, if it does not already exist.
    if not os.path.exists(exp_dir):
        os.makedirs(exp_dir)
    # Ensure that a copy of the parameters file is placed in the experimental directory.
    if parameters_file is not None and exp_dir not in parameters_file:
        destination_file_name = parameters_file
        if '/' in parameters_file:
            destination_file_name = parameters_file[parameters_file.rindex('/')+1 :]
        if destination_file_name not in os.listdir(exp_dir):
            # Check that the parameter file is not already there.
            shutil.copy(parameters_file, exp_dir + '/' + destination_file_name)

    # Ensure parameters are correctly specified.
    if mouse is not None and cageID is not None:
        raise Exception('Cannot specify both a mouse ID and a cage ID. Specify one or the other (or neither).')
    if not exp_dir:
        raise Exception("An experimental directory must be provided, using '-d' <dir name> command line argument.")
    if not mucinWater: raise Exception('Indicate if mucin water content is accounted for; -mw 1/0.')
    if not diet_code: raise Exception('Please specify a diet code. ')

    expMice = prepare_experiment(mucinWater=mucinWater, dietCode=diet_code, **kwargs)
    if dietAnalysis:
        diet_analysis(expMice)

    # If a specific mouse was requested, execute that in place of all of the mice.
    if cageID is not None:
        print('searching for mouse by cage ID. Requested ID = ' + str(cageID))
        # Search for the mouse with the specified cage id.
        for mouse_id, m in expMice.items(): 
            if str(m.cageID) == str(cageID):
                print('  which corresponds to mouse number ' + str(i))
                mouse = mouse_id
                break

    if mouse is not None:
        expMice = {mouse: expMice[mouse]}

    # Create directories to store mouse data and graphs, if they don't already exist.
    dataPath = exp_dir + '/mice-data-' + diet_code + tag
    graphPath = exp_dir + '/mice-graphs-' + diet_code + tag
    # Execution on cluster means two executions of this code can try to create directory simultaneously. 
    # No other way of checking this.
    try:
        if not os.path.exists(dataPath):
            os.makedirs(dataPath)
        if not os.path.exists(graphPath):
            os.makedirs(graphPath)
    except Exception:
        print('Directory already exists.')

    # Simulate mice.
    for index, (mouse_id, mouse) in enumerate(expMice.items()):
        if not seedSpecified:
            seed = index
        random.seed(seed)
        if index == 0:
            # Write the parameters to the filesystem.
            paramPath = dataPath + '/parameters-' + str(mouse._ID) + subMouse  # File to write parameter values to.
            if not os.path.exists(paramPath):
                write_params_fs(filename=paramPath)
        dPath = dataPath + '/mouse-' + str(mouse._ID) + subMouse
        timeseries_path = dataPath + '/timeseries/'        
        gPath = graphPath + '/mouse-' + str(mouse._ID) + subMouse

        if forceExecute or (not os.path.exists(dPath)):
            print('Mouse vital statistics')
            print(mouse.dietCode)
            print('wheatstarch = ' + str(mouse.daily_wheatstarch_kj))
            print('dextrinised = ' + str(mouse.daily_dextrinised_kj))
            print('sucrose = ' + str(mouse.daily_sucrose_kj))
            print('fat = ' + str(mouse.daily_fat_kj))
            print('nitrogen = ' + str(mouse.daily_nitrogen_kj))
            print('carbs = ' + str(mouse.daily_wheatstarch_kj + mouse.daily_dextrinised_kj + mouse.daily_sucrose_kj))
            print('daily_chow_g = ' + str(mouse.daily_chow_g))

            # 7 days = 168 hours
            bacteria_graphs_start_record_time = end_time - bacteria_graphs
            mouse.execute(endtime=end_time, timestep=timestep, verbose=verbose,
                          bacteria_graphing_start_time=bacteria_graphs_start_record_time, graphPath=gPath)

            xml_file = dataPath + '/output-' + str(mouse._ID) + subMouse + '.xml'

            mouse.logger.writeMouseToFile(dPath, xml_file)
            os.makedirs(timeseries_path, exist_ok=True)
            mouse_timeseries_path = timeseries_path + '/mouse-' + str(mouse._ID) + subMouse + '_timeseries.csv'
            mouse.logger.writeMouseTimeseriesToFile(mouse_timeseries_path)

            file_name = graphPath + '/mouse' + str(mouse._ID) + subMouse
            title = 'mouse ' + str(mouse._ID) + ', diet code ' + mouse.dietCode
            mouse.logger.plot_microbial_loads(filename=file_name, title=title)
            mouse.logger.graph_nutrient_inputs(filename=file_name+'nutrients')
            mouse.logger.plot_spatial_distribution(file_name + '-spatial')
            mouse.logger.plot_spatial_distribution_bacterial_states(file_name + '-spatialState')
            mouse.logger.plot_spatial_nutrient_distribution()
            # Heatmaps, showing states over time and colon length.
            if bacteria_graphs_start_record_time is not None:
                mouse.logger.plot_spatiotemporal_heatmaps()
            if mouseGraphs:
                mouse.logger.graph_guild_states(prefix=file_name)
    print("experiment completed.")


if __name__ == '__main__':
    main()
