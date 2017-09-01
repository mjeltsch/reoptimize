#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Example: reoptimize.py -e 'EcorI 2' 'HindIII 1' -l 3000 -t 2 -m 3
#
# All parameters are given via the command line:
#
# -e restriction enzymes to be used, followed by the number of target
#    sites in the target DNA
# -l length of target DNA: 3000 bp
# -t incubation time: 2 hours
# -m amount of target DNA: 3 micrograms
#
# Parameters l, t and m are optional. If they are missing,
# default values will be used as follows:
#
# length: 5000 bp, incubation time: 1 hour, amount of target DNA: 1 µg
#
# Example: reoptimize.py -e 'EcoRI 2' 'HindIII 1'
# (equivalent to reoptimize.py -e 'EcoRI 2' 'HindIII 1' -l 5000 -t 1 -m 1)
#
# argparse to parse the command line arguments and options
# sqlite3 module to store the enzyme data in  local sqlite database file
import sys, os, argparse, sqlite3
from argparse import RawTextHelpFormatter

# If this is set to True, much more info will be printed out during the run
DEBUG = False

# Get the directory where this script is located
path = os.path.dirname(os.path.realpath(__file__))
#print("Path: " + path)

# Function to check whether debugging info should be printed
def debug_print(*args):
    if DEBUG == True:
        for string in args:
            print(string, sep='\n')

assay_DNA_length = {  'λ': 48502,
                    'est 1 µg λ': 48502,
                    'est 1 µg λ DNA (HindIII digest)': 48502,
                    'pXba': 22563,
                    'pXba-XbaI': 22563, #hypothetical
                    'pNEB193': 2713,
                    'pNEB193-SrfI': 2715, #hypothetical
                    'T7': 39937,
                    'est 1 µg T7': 39937,
                    'Adeno-2': 35937,
                    'pBC4': 10673,
                    'supercoiled pUC19': 2686,
                    'ΦX174 RF I': 5386,
                    'pBR322': 4361,
                    'T4 wild-type phage': 168922 }

# The main digest is done here, receives the list of enzymes from command line
def digest(enzyme, microgram, length, time):

    # This is the name of the sqlite database file, that contains all the enzyme information
    sqlite_file= path + '/REsqlite3.db'
    #  Create sqlite connection
    sqlcon = sqlite3.connect(sqlite_file)
    # In the list_of_enzyme_activities everything is stored for later evaluation
    list_of_enzyme_activities = {}
    #
    units = {}
    # Each key (= buffer) has as value a list, which consists of
    # 1. whether digestion is allowed in that buffer (1= yes, 0 = no, -1 = unknown),
    #    set initially to 1 (except FastDigest buffer, which needs to be dealt with separately)
    # 2. The % activity in that buffer (set initially to 0)
    # It will be updated for every enzyme that is included in the digest
    # and is evaluated in the end
    buffer_list = {}

    # Test whether the sqlite database file is present
    # and get a list of all tables, that have the string "uffer"
    # in their name
    with sqlcon:
        cursor = sqlcon.cursor()
        query = "SELECT COUNT(*) FROM restriction_enzyme"
        try:
            cursor.execute(query)
            result = cursor.fetchone()
            debug_print("Number of enzymes in restriction_enzyme table: %s " % result)
        except sqlcon.Error as err:
            print("Error opening database file " + sqlite_file + ". Error: " + str(err))
        query = "SELECT name FROM sqlite_master WHERE type='table'"
        try:
            cursor.execute(query)
            result = cursor.fetchall()
            debug_print(result)
            for sql_table in result:
                debug_print(sql_table[0])
                if 'buffer' in sql_table or 'Buffer' in sql_table[0]:
                    buffer_list[sql_table[0]] = [1,0]
            debug_print("buffer_list (from sqlite file): " + str(buffer_list))
        except sqlcon.Error as err:
            print("Error opening database file " + sqlite_file + ". Error: " + str(err))

    how_many_enzymes = len(enzyme)
    # Generate the first empty entry into this dictionary, which will be filled
    # during the parsing of the enzyme entry from the sqlite database
    list_of_enzyme_activities = {}
    #
    # Loop through buffers to get all activity data for the enzyme
    #
    # THIS IS THE MAIN LOOP. EVERYTHING THAT NEEDS TO BE DONE FOR EVERY
    # ENZYME INCLUDED IN THE DIGEST NEEDS TO GO INTO THIS LOOP
    #
    for enzyme_item in enzyme:
        # Separate enzyme from number of cutting sites
        # enzyme[1] is name and enzyme[2] is number of cutting sites in assay DNA
        enzyme_item = enzyme_item.split(' ')
        # Check whether the user has submitted two values for each enzyme (name + number of cuts in taget site)
        # and whether the second value is an integer.
        # If no values were submitted, assume 1 restriction site
        if len(enzyme_item) < 2:
            number_of_restriction_sites = 1
        # Try to convert the second value (= number of restriction sites in taget DNA) into an integer
        else:
            try:
                number_of_restriction_sites = int(enzyme_item[1])
                #if not isinstance(number_of_restriction_sites, int):
            except ValueError:
                sys.exit("Please indicate the number of restriction sites after each enzyme name! Example: reoptimize -e 'AflIII 2' 'HindIII 1'")

        # Truncate too long enzyme names (the database holds only 32 character long enzyme names)
        enzyme_item[0] = enzyme_item[0][:32]
        # Make request to table "restriction_enzyme".
        # Convert all enzyme names to upper case first.
        query = "SELECT enzyme_id, default_buffer, assay_DNA, assay_DNA_cuts, survival, reaction_temperature, enzyme_name, reaction_supplement, enzyme_concentration, timesaver FROM restriction_enzyme WHERE UPPER(enzyme_name) = '" + enzyme_item[0].upper() + "'"
        try:
            cursor.execute(query)
            result = cursor.fetchone()
        except sqlcon.Error as err:
            error_message = "Error getting activity data for enzyme " + enzyme_item[0] + ". Error: " + str(err) + "\nQuery was: " + query
            sys.exit(error_message)
        # Check whether the enzyme was found from the database
        if result is None:
            sys.exit("There is no data for enzyme " + enzyme_item[0] + " in the database!")
        # Store all data in specific variables to free the result list variable
        enzyme_name = result[6]
        enzyme_id = result[0]
        default_buffer = result[1]
        assay_DNA = result[2]
        assay_DNA_cuts = int(result[3])
        survival = result[4]
        reaction_temperature = result[5]
        reaction_supplement = result[7]
        debug_print("Result 5: " + str(result[5]))
        # Take the lowest enzyme concentration
        enzyme_concentration = int(result[8].split(',')[1])
        try:
            timesaver = int(result[9])
        except:
            timesaver = ''
        debug_print("timesaver: " + str(timesaver))
        # Get all activity data for each standard buffer
        # Three-dimensional dictionary!
        list_of_enzyme_activities[enzyme_name] = {}
        debug_print("list_of_enzyme_activities: " + str(list_of_enzyme_activities))
        list_of_enzyme_activities[enzyme_name]['reaction_buffers'] = {}
        for buffer in buffer_list.keys():
            # The following query gets the following data (list of 2 items):
            # %-activity in current buffer, star activity in current buffer
            query = "SELECT `" + buffer + "`.activity, `" + buffer +  "`.star_activity FROM restriction_enzyme INNER JOIN `" + buffer + "` ON restriction_enzyme.enzyme_id = `" + buffer + "`.enzyme_id WHERE restriction_enzyme.enzyme_id = '" + str(enzyme_id) + "'"
            debug_print(query)
            try:
                cursor.execute(query)
                result = cursor.fetchone()
                # Add % activity
                list_of_enzyme_activities[enzyme_name]['reaction_buffers'][buffer] = result[0]
                debug_print(str(list_of_enzyme_activities))
                # Add star activity
                #list_of_enzyme_activities[enzyme[0]][buffer].append(result[1])
                # Only allow digest, if activity equal or greater than 50%
                if result[0] < 50:
                    buffer_list[buffer][0] = 0
                debug_print(enzyme_item[0] + " activity in " + buffer + ": " + str(result[0]) + ", star activity: " + str(result[1]))
                # Add cumulatively all % activities to be able to select the best buffer
                # if several are possible
                buffer_list[buffer][1] += result[0]
                # If there is star activity (or an unknown situation), disallow digest
                if result[1] != 0:
                     buffer_list[buffer][0] = 0
            except sqlcon.Error as err:
                error_message = "Error getting activity data for enzyme " + enzyme_item[0] + ". Error: " + str(err) + "\nQuery was: " + query
                sys.exit(error_message)

        #
        # Start calculating enzyme amounts here
        #
        list_of_enzyme_activities[enzyme_name]['units'] = microgram * int(number_of_restriction_sites) * assay_DNA_length[assay_DNA] / (length * assay_DNA_cuts)

        # Calculate the reduced enzyme amounts for digests > 1 hour.
        # If the enzyme survival is unknown (= 0) or if the enzyme does
        # not support longer survival times than 1 hour (= 1), don't do
        # anything for digests longer or equal to 1 hour.
        #
        if time >= 1:
            if survival == 8:
                # Formulas obtained empirically with NEB data using Matlab (rational function) regression
                fx = (0.05461*time+1.343)/(time+0.3991)
                list_of_enzyme_activities[enzyme_name]['units'] = list_of_enzyme_activities[enzyme_name]['units'] * fx
            elif survival == 4:
                fx = (0.1601*time+1.819)/(time+0.9845)
                list_of_enzyme_activities[enzyme_name]['units'] = list_of_enzyme_activities[enzyme_name]['units'] * fx
            elif survival == 2:
                fx = (0.4081*time+2.61)/(time+2.031)
                list_of_enzyme_activities[enzyme_name]['units'] = list_of_enzyme_activities[enzyme_name]['units'] * fx
        # Linear regression for interval 0-1 hour for all timesaver enzymes using
        # NEB data (assuming, that no enzyme is consumed during this short period).
        # For all other enzymes, just assume inverse proportionality
        # For some enzymes, this leads to paradoxical results, e.g.
        # for AvrII, where a 1 hour digests needs more enzyme than a 5 minute digest
        if time < 1:
            if timesaver != '':
                # convert timesaver into hours
                timesaver = timesaver/60
                debug_print("enzyme_concentration: " + str(enzyme_concentration))
                debug_print("timesaver: " + str(timesaver))
                debug_print("list_of_enzyme_activities[enzyme_name]['units']: " + str(list_of_enzyme_activities[enzyme_name]['units']))
                mm = (enzyme_concentration-list_of_enzyme_activities[enzyme_name]['units'])/(timesaver-1)
                debug_print("mm: " + str(mm))
                bb = list_of_enzyme_activities[enzyme_name]['units'] - mm
                debug_print("bb: " + str(bb))
                list_of_enzyme_activities[enzyme_name]['units'] = mm*time+bb
            else:
                list_of_enzyme_activities[enzyme_name]['units'] = list_of_enzyme_activities[enzyme_name]['units']/time
        # Add reaction temperatures to the list for a later comparison
        list_of_enzyme_activities[enzyme_name]['reaction_temperature'] = reaction_temperature

        # Add potential supplements to the list for a later comparison
        list_of_enzyme_activities[enzyme_name]['reaction_supplement'] = reaction_supplement

    # END OF "ENZYME IN ENZYMES" LOOP
    #
    # Make a list of possible buffers where the digest is allowed
    # Criteria (already cheked above):
    # - No star activity
    # - %-activity at least 50%
    #
    # Maybe remove this list later and sort the buffer_list dictionary instead
    # and drop from buffer_list all buffers that are not allowed
    #
    possible_buffers = []
    debug_print("buffer_list: " + str(buffer_list))
    for buffer, digest_allowed in buffer_list.items():
        if digest_allowed[0] == 1:
            possible_buffers.append([buffer, digest_allowed[1]])

    # Sort the list of possible buffers according to highest cumulative activity
    # Secondary sort key: name of buffer NOT YET IMPLEMENTED
    possible_buffers = sorted(possible_buffers, key = lambda number: number[1], reverse = True)
    debug_print("possible_buffers: " + str(possible_buffers))

    # Separate the input from the result by a blank line
    print("")

    # Put all reaction temperatures into a list and check whether
    # they are all the same, print warning if not
    # Print warning if a supplement is needed!
    temperature_list = []
    for restriction_enzyme, value in list_of_enzyme_activities.items():
        temperature_list.append(value['reaction_temperature'])
        if value['reaction_supplement'] != '':
            print("Note: " + value['reaction_supplement'] + "!")
    debug_print("Temperature list: " + str(temperature_list))
    if len(set(temperature_list)) != 1:
        print("The reaction temperatures of the enzymes are dfferent and you should perform a sequential digest!")

    #
    # PRINT THE RESULTS
    #
    # Print multiple digests
    #
    if len(possible_buffers) > 1:
        print("Digest is possible in the following buffers (avaraged % activity in brackets):")
        for buffer in possible_buffers:
            print("- " + str(buffer[0]) + " (" + str(round(buffer[1]/how_many_enzymes)) + ")", end = ' ')
            #
            # Print the amount of enzymes needed
            #
            debug_print(str(list_of_enzyme_activities))
            for restriction_enzyme, value in sorted(list_of_enzyme_activities.items()):
                percentage = value['reaction_buffers'][buffer[0]]
                if value['units']*100/percentage < 1:
                    rounded_units = str(round(value['units']*100/percentage, 2))
                else:
                    rounded_units = str(round(value['units']*100/percentage, 1))
                print(restriction_enzyme + ": " + rounded_units + " units", end = ' ')
            print("")
    #
    # Print single enzyme digests
    #
    elif len(possible_buffers) > 0:
        print("Digest is possible in the following buffer (", end = '')
        if len(units) > 1:
            print("averaged ", end='')
        print("% activity in brackets): ", end = "")
        for buffer in possible_buffers:
            print(str(buffer[0]) + " (" + str(round(buffer[1]/how_many_enzymes)) + ")")
            #
            # Print the amount of enzymes needed
            #
            debug_print(str(list_of_enzyme_activities))
            for restriction_enzyme, value in sorted(list_of_enzyme_activities.items()):
                percentage = value['reaction_buffers'][buffer[0]]
                if value['units']*100/percentage < 1:
                    rounded_units = str(round(value['units']*100/percentage, 2))
                else:
                    rounded_units = str(round(value['units']*100/percentage, 1))
                print(restriction_enzyme + ": " + rounded_units + " units", end = ' ')
            print("")
    #
    # Print not recommended digests
    #
    else:
        print("This ", end='')
        if how_many_enzymes == 2:
            print("double", end='')
        elif how_many_enzymes == 3:
            print("triple", end='')
        elif how_many_enzymes > 3:
            print("multiple", end='')
        print(" digest is not recommended.")

    debug_print(possible_buffers)
    debug_print(list_of_enzyme_activities)


def run():
    # Set up command line
    parser = argparse.ArgumentParser(description='reoptimize calculates enzyme amounts and possible buffers for restriction digests of DNA.\n\nUSAGE EXAMPLES:\n\nDigest a plasmid that has two EcoRI sites and one HindIII site with EcoRI and HindIII:\nreoptimize -e \'EcoRI 2\' \'HindIII 1\'\n\nDigest in 4 hours 4 µg of a 3000-bp plasmid that has 3 EcoRI sites and 5 HindIII sites with EcoRI and HindIII:\nreoptimize -e \'EcoRI 3\' \'HindIII 5\' -t 2 -l 3000 -m 4\n\nIf you don\'t specify time, target DNA length, DNA amount and number of restriction sites\ndefault values are assumed as follows:\n1 hour, 5000 bp, 1 µg, 1 restriction site/plasmid for all enzymes used', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-e','--enzyme', help='Restriction Enzyme', required=True, nargs='+')
    parser.add_argument('-m','--microgram', help='DNA amount (in µg)', default=1, type=float, nargs='?')
    parser.add_argument('-l','--length', help='Length of target DNA (in bp)', default=5000, type=int, nargs='?')
    parser.add_argument('-t','--time', help='Digestion time (in hours)', default=1, type=float, nargs='?')
    # Parse command line arguments
    args = vars(parser.parse_args())
    debug_print(args['enzyme'], args['microgram'], args['length'], args['time'])
    # Call the main function
    digest(args['enzyme'], args['microgram'], args['length'], args['time'])

if __name__ == '__main__':
    run()
