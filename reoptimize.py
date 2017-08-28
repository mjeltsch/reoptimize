#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Example: restriction.py digest -e EcoRI 2 -e HindIII 1
#
import sys
# sqlite3 module to store the enzyme data in  local sqlite database file
import sqlite3
# Click module to implement the command line functionality
import click

# If this is set to True, much more info will be printed out during the run
DEBUG = False

# Function to check whether debugging info should be printed
def debug_print(string):
    if DEBUG == True:
        print(string)

assay_DNA_length = {  'λ': 48502,
                    'est 1 µg λ': 48502,
                    'est 1 µg λ DNA (HindIII digest)': 48502,
                    'pXba': 22563,
                    #'pXba-XbaI': 22563, #?
                    'pNEB193': 2713,
                    #'pNEB193-SrfI': 2713, #?
                    'T7': 39937,
                    'est 1 µg T7': 39937,
                    'Adeno-2': 35937,
                    'pBC4': 10673,
                    'supercoiled pUC19': 2686,
                    'ΦX174 RF I': 5386,
                    'pBR322': 4361,
                    'T4 wild-type phage': 168922 }

@click.group()
def greet():
    pass

@greet.command()

# These are the enzymes that are following the "digest" command: filename.py digest EcoRI HindIII
@click.option('--enzymes', '-e', multiple=True, prompt='Enter enzymes like this: EcoRI 1 HindIII 2')

@click.option('--microgram', '-m', nargs=1, default=1, type=float, prompt='Enter DNA amount (in µg)')

@click.option('--length', '-l', nargs=1, default=5000, type=int, prompt='Length of target DNA (in bp)')

@click.option('--time', '-t', nargs=1, default=1, type=float, prompt='Enter digestion time (in hours)')

# The main digest is done here, receives the list of enzymes from command line
def digest(enzymes, microgram, length, time):

    # This is the name of the sqlite database file, that contains all the enzyme information
    sqlite_file='REsqlite3.db'

    #  Create sqlite connection
    sqlcon = sqlite3.connect(sqlite_file)

    # This list is not yet needed (but needed for features that will be later implemented)
    list_of_enzyme_activities = {}
    units = {}

    # Each key (= buffer) has as value a list, which consists of
    # 1. whether digestion is allowed in that buffer (1= yes, 0 = no, -1 = unknown),
    #    set initially to 1 (except FastDigest buffer, which needs to be dealt with separately)
    # 2. The % activity in that buffer (set initially to 0)
    # It will be updated for every enzyme that is included in the digest
    # and is evaluated in the end
    buffer_list = {'NEBuffer 1.1': [1,0], 'NEBuffer 2.1': [1,0], 'NEBuffer 3.1': [1,0], 'CutSmart® Buffer': [1,0], 'FastDigest buffer': [-1,0]}

    # Test where the sqlite database file is present
    with sqlcon:
        cursor = sqlcon.cursor()
        query = "SELECT COUNT(*) FROM restriction_enzyme"
        try:
            cursor.execute(query)
            result = cursor.fetchone()
            debug_print("Number of enzymes in restriction_enzyme table: %s " % result)
        except sqlcon.Error as err:
            print("Error opening database file " + sqlite_file + ". Error: " + str(err))

    # The "Click" module gives a list with every character as a single list item when prompting is used!
    # When the enzymes are giving via the command line, each enzyme is a list item
    # Hence, we have to determine whether the enzymes were given as a command line
    # argument or whether they were prompted for. We do this by checking the length
    # of each individual item in the enzyme list, which should never be longer than 1
    # if the promting was used.
    #
    # enzymes = list of enzymes, that is passed by the command line module Click
    #
    for enzyme in enzymes:
        if len(enzyme) > 1:
            prompt = False
        else:
            prompt = True
    if prompt == True:
        # Make a single string from the character list,
        # that is passed by the prompt
        enzymes = ''.join(enzymes)
        # Now split at every space and rejoin in pairs of two
        enzymes = [name+" "+number for name,number in zip(*[iter(enzymes.split(' '))]*2)]
        debug_print(enzymes)
    else:
        debug_print("Enzymes: " + str(enzymes))
    how_many_enzymes = len(enzymes)
    # Generate the first empty entry into this dictionary, which will be filled
    # during the parsing of the enzyme entry from the sqlite database
    list_of_enzyme_activities = {}
    #
    # Loop through buffers to get all activity data for the enzyme
    #
    # THIS IS THE MAIN LOOP. EVERYTHING THAT NEEDS TO BE DONE FOR EVERY
    # ENZYME INCLUDED IN THE DIGEST NEEDS TO GO INTO THIS LOOP
    #
    for enzyme in enzymes:
        # Separate enzyme from number of cutting sites
        # enzyme[1] is name and enzyme[2] is number of cutting sites in assay DNA
        enzyme = enzyme.split(' ')
        debug_print("enzyme[0]: " + enzyme[0])
        debug_print("enzyme[1]: " + enzyme[1])
        list_of_enzyme_activities[enzyme[0]] = {}
        debug_print("list_of_enzyme_activities: " + str(list_of_enzyme_activities))
        # Truncate too long enzyme names (the database holds only 32 character long enzyme names)
        enzyme[0] = enzyme[0][:32]
        # Make request to table "restriction_enzyme"
        query = "SELECT enzyme_id, default_buffer, assay_DNA, assay_DNA_cuts, survival, reaction_temperature FROM restriction_enzyme WHERE enzyme_name = '" + enzyme[0] + "'"
        try:
            cursor.execute(query)
            result = cursor.fetchone()
        except sqlcon.Error as err:
            error_message = "Error getting activity data for enzyme " + enzyme[0] + ". Error: " + str(err) + "\nQuery was: " + query
            sys.exit(error_message)
        # Check whether the enzyme was found from the database
        if result is None:
            sys.exit("There is no data for enzyme " + enzyme[0] + " in the database!")
        # Store all data in specific variables to free the result list variable
        enzyme_id = result[0]
        default_buffer = result[1]
        assay_DNA = result[2]
        assay_DNA_cuts = int(result[3])
        survival = result[4]
        reaction_temperature = result[5]
        # Get all activity data for each buffer
        # Three-dimensional dictionary!
        list_of_enzyme_activities[enzyme[0]]['reaction_buffers'] = {}
        for buffer in buffer_list.keys():
            # The following query gets the following data (list of 2 items):
            # %-activity in current buffer, star activity in current buffer
            query = "SELECT `" + buffer + "`.activity, `" + buffer +  "`.star_activity FROM restriction_enzyme INNER JOIN `" + buffer + "` ON restriction_enzyme.enzyme_id = `" + buffer + "`.enzyme_id WHERE restriction_enzyme.enzyme_id = '" + str(enzyme_id) + "'"
            debug_print(query)
            try:
                cursor.execute(query)
                result = cursor.fetchone()
                # Add % activity
                list_of_enzyme_activities[enzyme[0]]['reaction_buffers'][buffer] = result[0]
                debug_print(str(list_of_enzyme_activities))
                # Add star activity
                #list_of_enzyme_activities[enzyme[0]][buffer].append(result[1])
                # Only allow digest, if activity equal or greater than 50%
                if result[0] < 50:
                    buffer_list[buffer][0] = 0
                debug_print(enzyme[0] + " activity in " + buffer + ": " + str(result[0]) + ", star activity: " + str(result[1]))
                # Add cumulatively all % activities to be able to select the best buffer
                # if several are possible
                buffer_list[buffer][1] += result[0]
                # If there is star activity (or an unknown situation), disallow digest
                if result[1] != 0:
                     buffer_list[buffer][0] = 0
            except sqlcon.Error as err:
                error_message = "Error getting activity data for enzyme " + enzyme[0] + ". Error: " + str(err) + "\nQuery was: " + query
                sys.exit(error_message)
        #
        # Start calculating enzyme amounts here
        #
        list_of_enzyme_activities[enzyme[0]]['units'] = microgram * int(enzyme[1]) * assay_DNA_length[assay_DNA] / (length * assay_DNA_cuts)

        # Calculate the reduced enzyme amounts for digests != 1 hour.
        # If the enzyme survival is unknown (= 0) or if the enzyme does
        # not support longer survival times than 1 hour (= 1), don't do
        # anything!
        #
        if survival == 8:
            # Formulas obtained empirically with NEB data using Matlab (rational function) regression
            fx = (0.05461*time+1.343)/(time+0.3991)
            list_of_enzyme_activities[enzyme[0]]['units'] = list_of_enzyme_activities[enzyme[0]]['units'] * fx
        elif survival == 4:
            fx = (0.1601*time+1.819)/(time+0.9845)
            list_of_enzyme_activities[enzyme[0]]['units'] = list_of_enzyme_activities[enzyme[0]]['units'] * fx
        elif survival == 2:
            fx = (0.4081*time+2.61)/(time+2.031)
            list_of_enzyme_activities[enzyme[0]]['units'] = list_of_enzyme_activities[enzyme[0]]['units'] * fx

    # END OF "ENZYME IN ENZYMES" LOOP
    #
    # Make a list of possible buffers where the digest is allowed
    # Criteria (already cheked above):
    # - No star activity
    # - %-activity at least 50%
    #
    possible_buffers = []
    debug_print("buffer_list: " + str(buffer_list))
    for buffer, digest_allowed in buffer_list.items():
        if digest_allowed[0] == 1:
            possible_buffers.append([buffer, digest_allowed[1]])

    # Append the default buffer of an enzyme to the list of possible buffers
    # if a single enzyme digest was requested and the buffer is not yet
    # present in the list (this is needed to display unique/non-standard buffers)
    debug_print(str(possible_buffers))
    if how_many_enzymes == 1:
        if [default_buffer, 100] not in possible_buffers:
            possible_buffers.append([default_buffer, 100])

    # Sort the list of possible buffers according to highest cumulative activity
    possible_buffers = sorted(possible_buffers, key = lambda number: number[1], reverse = True)

    #
    # PRINT THE RESULTS
    #
    # Print multiple digests
    #
    if len(possible_buffers) > 1:
        print("Digest is possible in the following buffers (avaraged % activity in brackets):")
        for buffer in possible_buffers:
            print(str(buffer[0]) + " (" + str(round(buffer[1]/how_many_enzymes)) + ")", end = ' ')
            #
            # Print the amount of enzymes needed
            #
            for restriction_enzyme, value in list_of_enzyme_activities.items():
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

@greet.command()
def test(**kwargs):
    print("This is only a test.")

if __name__ == '__main__':
    greet()
