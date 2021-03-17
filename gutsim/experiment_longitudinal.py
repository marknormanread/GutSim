__author__ = 'Mark N. Read'

import gutsim.experiment  # contains utility methods used here.
import gutsim.mouse
import gutsim.bacteria
import gutsim.feed_regimen
import os
import sys
import xml.etree.ElementTree as ET
import shutil
import random
import pandas



def prepare_experiment(intake_table, mouse_id):
    """

    :param intake_table:
    :param mouse_id: string, entry in the table's 'MouseID' column. Can appear multiple times, signifying several
     distinct feeding regimens at different points in time.
    :return:
    """
    intakes = pandas.read_excel(io=intake_table, header=0)  # returns pandas dataframe
    print('mouseID = ' + str(mouse_id))
    selected_rows = intakes.loc[:, 'MouseID'] == mouse_id  # returns array of boolean
    intake_records = intakes.loc[selected_rows, :]
    nutrient_inputs = []
    for record_row in intake_records.index:
        regimen = intake_records.loc[record_row, 'Intake_pattern']
        regimen_start_h = intake_records.loc[record_row, 'Start_time_days'] * 24.
        regimen_end_h = intake_records.loc[record_row, 'End_time_exclusive_days'] * 24.
        daily_nitrogen_g = intake_records.loc[record_row, 'Protein_g/d']
        daily_wheatstarch_g = intake_records.loc[record_row, 'Wheatstarch_g/d']
        daily_dextrinised_g = intake_records.loc[record_row, 'Dextrinised_wheatstarch_g/d']
        daily_sucrose_g = intake_records.loc[record_row, 'Sucrose_g/d']
        daily_inulin_g = intake_records.loc[record_row, 'Inulin_g/d']

        # constant input over 24h
        if regimen == 'const':
            nutrient_inputs.append(gutsim.feed_regimen.Constant(
                dailyCw=daily_wheatstarch_g, dailyCd=daily_dextrinised_g, dailyCs=daily_sucrose_g,
                dailyNf=daily_nitrogen_g,
                regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h))
        # time restricted feeding.
        if regimen == 'TRF':
            trf_start_time = intake_records.loc[record_row, 'TRF_start_h']
            trf_duration = intake_records.loc[record_row, 'TRF_duration_h']
            nutrient_inputs.append(gutsim.feed_regimen.TRF(
                start_time=trf_start_time, colicInputDuration=trf_duration,
                dailyCw=daily_wheatstarch_g, dailyCd=daily_dextrinised_g, dailyCs=daily_sucrose_g,
                dailyCi=daily_inulin_g, dailyCi_inverse=0.,
                dailyNf=daily_nitrogen_g,
                regimen_start_h=regimen_start_h, regimen_end_h=regimen_end_h))
    # add mucin secretions to the input.
    daily_mucinC_g, daily_mucinN_g = gutsim.experiment.calculate_daily_mucin_secretions(accWater=True)
    nutrient_inputs.append(gutsim.feed_regimen.Mucin(dailyCm=daily_mucinC_g, dailyNm=daily_mucinN_g))
    mouse = gutsim.mouse.Mouse(nutrientInputs=nutrient_inputs, dietCode='multiple',
                               daily_wheatstarch_kj=0,  # KJ values can't be given, multiple regimens.
                               daily_dextrinised_kj=0, daily_sucrose_kj=0,
                               daily_nitrogen_kj=0,
                               daily_fat_kj=0,
                               daily_chow_g=0,
                               ID=0, cageID=mouse_id)
    return mouse


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
    parameters_file = None  # location of parameters file, optional.
    intake_table = None  # location of excel file containing breakdown of mouse's longitudinal dietary intake. Required
    timestep = 0.1  # duration of time represented in each iteration of simulation, in hours.
    mouse_id = None  # mouse unique identifier in excel spreadsheet
    seed = None
    seedSpecified = False  # if user provided seed, set to true
    subMouse = ''  # for multiple runs of the same mouse. eg, file names "mouse_0.XXX.png"
    exp_dir = None  # where to write output data
    tag = ''  # text appended to end of standard generated data dirs. mice-data-24-3_TAG
    verbose = False	 # when true, write simulation status to standard out frequently.
    mucinWater = True  # when true, mucin nutrients are adjusted for water content.
    bacteria_graphs = False
    mouseGraphs = False  # draw graphs of mouse (and guild) level status.
    force_execute = False
    end_time = 168.0  # simulation end time, run for a week by default.


    # read in command line args.
    i = 1  # index into args. Ignore the first one, since it's the name of the python module being run.
    while i < len(sys.argv):
        if sys.argv[i] == '-p':
            i += 1
            parameters_file = sys.argv[i]
        elif sys.argv[i] == '-intake':
            i += 1
            intake_table = sys.argv[i]
        elif sys.argv[i] == '-m':  # mouse ID to be run.
            i += 1
            mouse_id = sys.argv[i]
        elif sys.argv[i] == '-sm':  # for multiple runs of the same mouse. eg, file names "mouse_0.XXX.png"
            i += 1                  # the '.' is added here, does not need to be specified.
            subMouse = '_' + sys.argv[i]
        elif sys.argv[i] == '-d':  # the experimental directory to which data should be written.
            i += 1
            exp_dir = sys.argv[i]
        elif sys.argv[i] == '-seed':  # set the random number seed.
            i += 1
            seed = int(sys.argv[i])
            seedSpecified = True
        elif sys.argv[i] == '-tag':	 # text to append to the end of standard generated data folders.
            i += 1
            tag = '_' + sys.argv[i]
        elif sys.argv[i] == '-v':  # verbose output to standard out.
            verbose = True
        elif sys.argv[i] == '-bg':  # bacteria graphs, whether to graph bacterial population dyanmics.
            bacteria_graphs = True
        elif sys.argv[i] == '-mg':
            mouseGraphs = True
        elif sys.argv[i] == '-end':	 # simulation termination time.
            i += 1
            end_time = float(sys.argv[i])
        elif sys.argv[i] == '-f':  # force the simulation to execute, even if data file already exists.
            force_execute = True
        elif sys.argv[i] == '-s':  # the scale factor, higher numbers lead to larger bacterial populations.
            i += 1                 # supplied value is multiplied by a million.
            gutsim.mouse.set_scale_factor(float(sys.argv[i]) * 1000000.0)
        i += 1

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
            force_execute = bool(int(node.text))
        node = root.find('./mouse')
        if node is not None and node.text is not None:
            # if nothing specified, defer to command line argument.
            if node.text.lower() == 'all':
                mouse = None   # ensures all mice are executed.
            else:
                mouse = int(node.text)
        node = root.find('./output/bacteria_graphs')
        if node is not None and node.text is not None:
            bacteria_graphs = bool(int(node.text))
        node = root.find('./output/mouse_graphs')
        if node is not None and node.text is not None:
            mouseGraphs = bool(int(node.text))

        # load parameters for other parts of the simulation. Each module handles it's own parameters.
        gutsim.mouse.load_parameters(root)
        gutsim.bacteria.load_parameters(root)
        print('parameters file specified: ' + parameters_file)
        print('experimental directory: ' + exp_dir)
    # create experimental directory, if it does not already exist.
    if not os.path.exists(exp_dir):
        os.makedirs(exp_dir)
    # ensure that a copy of the parameters file is placed in the experimental directory.
    if parameters_file is not None and exp_dir not in parameters_file:
        destination_file_name = parameters_file
        if '/' in parameters_file:
            destination_file_name = parameters_file[parameters_file.rindex('/')+1 :]
        if destination_file_name not in os.listdir(exp_dir):
            # check that the parameter file is not already there.
            shutil.copy(parameters_file, exp_dir + '/' + destination_file_name)

    if not exp_dir:
        raise Exception("An experimental directory must be provided, using '-d' <dir name> command line argument.")

    mouse = prepare_experiment(intake_table=intake_table, mouse_id=mouse_id)

    # create directories to store mouse data and graphs, if they don't already exist.
    print(exp_dir)
    print(mouse_id)
    print(tag)
    data_path = exp_dir + '/mice-data-' + mouse_id + tag
    graph_path = exp_dir + '/mice-graphs-' + mouse_id + tag
    # execution on cluster means two executions of this code can try to create directory simultaneously. No other way of
    # checking this.
    try:
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        if not os.path.exists(graph_path):
            os.makedirs(graph_path)
    except Exception:
        print('Directory already exists.')

    if not seedSpecified:
        seed = 999
    random.seed(seed)

    # write the parameters to the filesystem.
    paramPath = data_path + '/parameters-' + str(mouse._ID) + subMouse  # file to write parameter values to.
    if not os.path.exists(paramPath):
        gutsim.experiment.write_params_fs(filename=paramPath)
    dPath = data_path + '/mouse-' + str(mouse._ID) + subMouse
    gPath = graph_path + '/mouse-' + str(mouse._ID) + subMouse

    print("here")
    print(force_execute)
    if force_execute or (not os.path.exists(dPath)):
        print('Mouse vital statistics')
        # 7 days = 168 hours
        mouse.execute(endtime=end_time, timestep=timestep, verbose=verbose,
                      bacteria_graphing_start_time=bacteria_graphs, graphPath=gPath)

        xml_file = data_path + '/output-' + str(mouse._ID) + subMouse + '.xml'

        mouse.logger.writeMouseToFile(dPath, xml_file)
        mouse.logger.writeMouseTimeseriesToFile(dPath + '_timeseries')

        file_name = graph_path + '/mouse' + str(mouse._ID) + subMouse
        title = 'mouse ' + str(mouse._ID) + ', diet code ' + mouse.dietCode
        mouse.logger.plot_microbial_loads(filename=file_name, title=title)
        mouse.logger.graph_nutrient_inputs(filename=file_name+'nutrients')
        mouse.logger.plot_spatial_distribution(file_name + '-spatial')
        mouse.logger.plot_spatial_distribution_bacterial_states(file_name + '-spatialState')
        if mouseGraphs:
            mouse.logger.graph_guild_states(prefix=file_name)
    print("experiment completed.")


if __name__ == '__main__':
    main()
