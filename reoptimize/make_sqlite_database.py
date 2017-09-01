#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Example: restriction.py digest -e EcoRI -e HindIII
#
# glob => needed to specify filepattern *.gb when reading assay DNA sequences
import os, urllib3, shutil, re, glob
# sqlite3 module to store the enzyme data in  local sqlite database file
import sqlite3
# Click module to implement the command line functionality

# Use Biopython to calculate restriction enzyme frequencies
# for the assay DNAs

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import Restriction
from Bio import SeqIO
from Bio.Restriction import *
amb = IUPACAmbiguousDNA()

# To switch off warning due to unverified https request
urllib3.disable_warnings()

DEBUG = True

# Uncomment if you don't want to download the whole data set
#limit = '10'
#offset = '270'

def debug_print(string):
    if DEBUG == False:
        print(string)

# Test whether an enzyme is already in the database
def is_duplicate(vendor, enzyme_name, enzyme_url):
    query = "SELECT COUNT(*) FROM restriction_enzyme WHERE vendor = '" + vendor + "' AND enzyme_name = '" + enzyme_name + "' AND enzyme_url = '" + enzyme_url + "'"
    try:
        c.execute(query)
        result = c.fetchone()
        #print(query)
        #print("Result: " + str(result))
        if result[0] == 0:
            return (False, 0)
        else:
            return (True, 0)
    except sqlcon.Error as err:
        return(None, err)

def fix_enzyme_name(enzyme_name):
    #if enzyme_name[-5:] == '&reg;':
    #    enzyme_name = enzyme_name[:-5]
    #
    # The alpha designation seems to be similar to the HF
    # designation (what is the difference????)
    enzyme_name = enzyme_name.replace('<sup>&alpha;</sup>', '-alpha')
    enzyme_name = enzyme_name.replace('<sup>&reg;</sup>', ' (reg)')
    enzyme_name = enzyme_name.replace('&reg;', ' (reg)')
    # For CviKI-1
    enzyme_name = enzyme_name.replace('I-1', 'I_1')

    #enzyme_name = re.sub('\s*&reg;', ' (reg)', enzyme_name)
    return enzyme_name

# Make numbers from string activity definitions in html file
def activities(string):
    if string[-4:] == ' 0%':
        return 0
    elif string[-4:] == ' 10%':
        return 10
    elif string[-4:] == ' 25%':
        return 25
    elif string[-4:] == ' 50%':
        return 50
    elif string[-4:] == ' 75%':
        return 75
    elif string[-4:] == '100%':
        return 100
    else:
        return -1

# Remove html code and return plain text
def strip_html(string):
    return re.sub('<[^<]+?>', '', string)

# Biopython doesn't know about NEB's High Fidelity
# enzymes and the designation needs to be stripped
# before calling any Restriction method.
# The same for "-alpha"!
# Some enzymes (e.g. EcoP15I) are so new, that
# the set of restriction enzymes that biopython
# uses is too old. It can be updated manually,
# but e.g. the enzyme EcoP15I is only in the
# emboss_e.708 file, but only in the embossa_e.708
# file. What is the difference between these files?
# Compiling a new restriction dictionary with the
# embossa_e.708 file seems to work fine...
# How to do this in the final distribution?
# For Ubuntu 16.04 (biopython via ubuntu),
# the updated list goes to:
# /usr/lib/python3/dist-packages/Bio/Restriction
# Updating is done via a script, which is
# not included in the Ubuntu package,
# but in the git repo (git clone https://github.com/biopython/biopython)
# in the directory
# /home/mjeltsch/Documents/Informatics/Python/biopython/Scripts/Restriction/
#
def strip_HF_designation(string):
    tmp = re.sub('-HF(\s\(reg\))*\s*', '', string)
    return re.sub('-alpha', '', tmp)

# This is the name of the sqlite database file, that contains all the enzyme information
sqlite_file='REsqlite3.db'

# Define survival classes (0 = unknown) for enzymes
# NEB uses the +++/++/+/- designation which is not good
# for mathematical operations
survival_list = {'+++': 8, '++': 4, '+': 2, '-': 1, '': 0}

# Commands to make the database
#
# enzyme concentration is a list of sold concentrations (given in units/µl)
# timesaver: 15 or 5 (minutes needed for complete digestion of 1 µg DNA with 1 µl enzyme)
#
# For buffer tables:
#
# activity: 0 to 100%
# star activity: 1 = yes, 0 = no
#
textstring = '''DROP TABLE IF EXISTS `restriction_enzyme`\n
CREATE TABLE `restriction_enzyme` (\
`enzyme_id` mediumint(9) NOT NULL,\
`vendor` varchar(32) NOT NULL,\
`enzyme_name` varchar(32) NOT NULL,\
`enzyme_url` varchar(255) NOT NULL,\
`reaction_temperature` int(255),\
`storage_temperature` int(255),\
`default_buffer` varchar(32),\
`assay_DNA` varchar(32),\
`survival` int(2),\
`assay_DNA_cuts` mediumint(9),\
`reaction_supplement` varchar(255),\
`enzyme_concentration` varchar(12),\
`timesaver` int[2],\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `NEBuffer 1.1`\n
CREATE TABLE `NEBuffer 1.1` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `NEBuffer 2.1`\n
CREATE TABLE `NEBuffer 2.1` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `NEBuffer 3.1`\n
CREATE TABLE `NEBuffer 3.1` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `CutSmart® Buffer`\n
CREATE TABLE `CutSmart® Buffer` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `NEBuffer EcoRI`\n
CREATE TABLE `NEBuffer EcoRI` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n
DROP TABLE IF EXISTS `FastDigest buffer`\n
CREATE TABLE `FastDigest buffer` (\
`enzyme_id` mediumint(3) NOT NULL,\
`activity` int(3),\
`star_activity` BOOLEAN,\
 PRIMARY KEY (`enzyme_id`))\n'''

sqlcon = sqlite3.connect(sqlite_file)
c = sqlcon.cursor()
for line in textstring.splitlines():
    #print(text)
    c.execute(line)
    sqlcon.commit()

# Get html page of complete list of NEB restriction enzymes
http = urllib3.PoolManager()
url = "https://www.neb.com/products/restriction-endonucleases"
vendor = 'NEB'
textstring = http.request('GET', url).data.decode('utf-8')

# Count how many enzymes are inserted into the database
count = 0
# Flag to start/stop parsing the html file
start = 0
for line in textstring.splitlines():
    # Start here to search
    if line == "\t\t\t\tRestriction Endonucleases: A":
        start = 1
    # End the search here
    elif line[-47:] == "Restriction Endonuclease Buffers &amp; Diluents":
        start = 0
    #If the first 70 characters of a line match
    if (start == 1) and (line[:70] == "\t\t\t\t\t<span class=\"decorate order open\">Order</span><a href=\"/products/"):
        enzyme_url = line[70:].split('">')
        print('\n')
        enzyme_name = enzyme_url[1][:-4]
        enzyme_name = fix_enzyme_name(enzyme_name)
        enzyme_url = "https://www.neb.com/products/" + enzyme_url[0]
        # Test whether the enzyme is already in the database
        if is_duplicate(vendor, enzyme_name, enzyme_url)[0] == True:
            print(enzyme_url + ", " + enzyme_name + " is already in the database.")
        # If the enzyme is not a duplicate, enter it into the database
        else:
            query = "INSERT INTO restriction_enzyme (enzyme_id, vendor, enzyme_name, enzyme_url) VALUES (" + str(count) + ", '" + vendor + "', '" + enzyme_name + "', '" + enzyme_url + "')"
            try:
                c.execute(query)
                sqlcon.commit()
                count += 1
                print(str(count) + ". " + enzyme_url + ", " + enzyme_name)
            except sqlcon.Error as err:
                print("Error inserting enzyme " + enzyme_name + " to database. Error: " + str(err) + "\nQuery was: " + query)

print(str(count) + " enzymes inserted into the database.")

#
# PART 3: Getting enzyme data for NEB restriction enzymes
#
# DATA SOURCES
# Survival data: https://www.neb.com/tools-and-resources/usage-guidelines/restriction-endonucleases-survival-in-a-reaction
# Frequency of restriction sites in assay DNA: https://www.neb.com/tools-and-resources/selection-charts/frequencies-of-restriction-sites
# All other data: The specific enzyme page by NEB as listed in the sqlite "restriction_enzyme" table
#
# Uncomment if you don't want to download the whole data set
#limit = '10'
#offset = '0'
#
# Connect to database and retrieve all enzymes, for which we need to get the data
try:
    limit, offset
except NameError:
    query = "SELECT enzyme_id, enzyme_name, enzyme_url FROM restriction_enzyme"
else:
    query = "SELECT enzyme_id, enzyme_name, enzyme_url FROM restriction_enzyme LIMIT " + limit + " OFFSET " + offset
try:
    c.execute(query)
    result = c.fetchall()
except sqlcon.Error as err:
    print("Error connecting to database. Error: " + str(err))

# Count for how many enzymes data is retrieved and inserted into the database
count_enzymes = 0
count_buffer_entries = 0

# Set all enzyme activities and star activity to 'unknown' (= -1)
enzyme_activity = {'NEBuffer 1.1': [-1,-1], 'NEBuffer 2.1': [-1, -1], 'NEBuffer 3.1': [-1, -1], 'CutSmart® Buffer': [-1, -1], 'FastDigest buffer': [-1, -1], 'NEBuffer EcoRI': [-1, -1]}

# Get the survival table from NEB
url1 = 'https://www.neb.com/tools-and-resources/usage-guidelines/restriction-endonucleases-survival-in-a-reaction'
survivaltext = http.request('GET', url1).data.decode('utf-8')

# Get the timesaver table from NEB
url1 = 'https://www.neb.com/tools-and-resources/selection-charts/time-saver-qualified-restriction-enzymes'
timesavertext = http.request('GET', url1).data.decode('utf-8')

# For these enzymes, do not attempt to retrieve assay DNA
enzyme_blacklist_assay = ['McrBC']

# For these enzymes, do not attempt to retrieve survival data
enzyme_blacklist_survival = ['McrBC']

# For these enzymes, do not attempt to retrieve frequency data
enzyme_blacklist_frequency = []

# For these enzymes, do not attempt to timesaver data
enzyme_blacklist_timesaver = []

for enzyme in result:
    # Get html page for the specific enzyme
    url = enzyme[2]
    print("\n" + str(enzyme[0]) + ". " + enzyme[1] + ":")
    textstring = http.request('GET', url).data.decode('utf-8')
    #debug_print(textstring)
    previousline = ''
    reaction_supplement = ''
    # If no notes are found ("\t\t<li id=\"note-", used below in code), assume that there is no star activity
    # except for special buffers
    enzyme_activity['NEBuffer 1.1'][1] = 0
    enzyme_activity['NEBuffer 2.1'][1] = 0
    enzyme_activity['NEBuffer 3.1'][1] = 0
    enzyme_activity['CutSmart® Buffer'][1] = 0
    enzyme_activity['FastDigest buffer'][1] = -1
    enzyme_activity['NEBuffer EcoRI'][1] = -1
    # Make an empty enzyme concentration list (needed, since some producers
    # sell the same enzyme at different concentrations
    enzyme_concentration = []
    # Iterate through the lines of the html code
    for line in textstring.splitlines():
        # If the previous line matches
        #
        # TEMPERATURE AND SUPPLEMENTS
        #
        if previousline == "\tReaction Conditions":
            debug_print("Reaction conditions line: " + line)
            # Split the line of reaction conditions (= the line below the heading "Reaction Conditions")
            tmp = re.split('</h4><p>|<br />|</p>', line)
            enzyme_buffer = tmp[1][3:]
            # Check how many items:
            # 4 items: only temperature is given as reaction condition
            if len(tmp) == 4:
                reaction_temperature = tmp[2][12:]
            # 5 items: some additional special reaction condion is mentioned (mostly SAM)
            elif len(tmp) == 5:
                reaction_temperature = tmp[3][12:]
                reaction_supplement = tmp[2]
            # Fix the bold reaction temperatures
            if "<b>" in reaction_temperature:
                reaction_temperature = re.split('<b>|</b>', reaction_temperature)[1]
        #
        # BUFFERS
        #
        # Extract the line, where the activities in the standard buffers are listed
        elif previousline == "\tActivity in NEBuffers":
            line = line.replace('<strong>', '')
            line = line.replace('</strong>', '')
            activity = re.split('</h4>|<br />', line)
            activity.pop(0)
            debug_print("Activity extraction: " + str(activity))
            enzyme_activity['NEBuffer 1.1'][0] = activities(activity[0])
            enzyme_activity['NEBuffer 2.1'][0] = activities(activity[1])
            enzyme_activity['NEBuffer 3.1'][0] = activities(activity[2])
            enzyme_activity['CutSmart® Buffer'][0] = activities(activity[3])
            debug_print(str(enzyme_activity))
        #
        # STAR ACTIVITY
        #
        # This looks for star activity, which is given in the html file after
        # the following string in one line: "\t\t<li id=\"note-".
        # Put star activity into the enzyme_activity dictionary.
        # enzyme_activity is a dictionary where each key represents
        # one buffer and the corresponding value is a list with two elements:
        # 1. Activity of the enzyme in that buffer in % (0-100). If unknown,
        #    the value is -1
        # 2. Info whether the enzyme exhibits star activity in that buffer
        #    (yes = 1, no = 0, unknown = -1)
        if "\t\t<li id=\"note-" in line:
            if " 1.1" in strip_html(line):
                enzyme_activity['NEBuffer 1.1'][1] = 1
            else:
                enzyme_activity['NEBuffer 1.1'][1] = 0
            if " 2.1" in strip_html(line):
                enzyme_activity['NEBuffer 2.1'][1] = 1
            else:
                enzyme_activity['NEBuffer 2.1'][1] = 0
            if " 3.1" in strip_html(line):
                enzyme_activity['NEBuffer 3.1'][1] = 1
            else:
                enzyme_activity['NEBuffer 3.1'][1] = 0
            if " CutSmart" in strip_html(line):
                enzyme_activity['CutSmart® Buffer'][1] = 1
            else:
                enzyme_activity['CutSmart® Buffer'][1] = 0
            debug_print("Star activity extraction: " + str(enzyme_activity))
        # Keep the current line as variable "previousline" in memory
        # when the next line is analyzed
        #
        # ASSAY DNA
        #
        if enzyme[1] in enzyme_blacklist_assay:
            assay_DNA = 'unknown'
        else:
            if previousline == "\tUnit Definition":
                debug_print(line)
                # Maybe this search is too general and will give some false positives...
                re_result = re.search("(g of)(.*)(DNA)", line)
                # Use this less specific pattern only if the pattern above fails
                re_result2 = re.search("(g)(.*)(DNA)", line)
                if re_result is not None:
                    debug_print(str(re_result))
                    # Remove all whitespaces from the string
                    assay_DNA = re_result.group(2).strip()
                elif re_result2 is not None:
                    debug_print(str(re_result2))
                    # Remove all whitespaces from the string
                    assay_DNA = re_result2.group(2).strip()
                else:
                    assay_DNA = 'unknown'
        #
        # ENZYME  CONCENTRATION
        #
        m = re.search("(units</td><td>)(.*)( units/ml</td><td class=)", line)
        if str(m) != 'None':
            # Only write to the enzyme concentration list if the value is not
            # already in th elist
            if int(m.group(2).split(',')[0]) not in enzyme_concentration:
                enzyme_concentration.append(int(m.group(2).split(',')[0]))

        previousline = line
    #
    # RETRIEVE SURVIVAL AFTER THE MAIN ENZYME PAGE HAS BEEN SCRAPED
    #
    if enzyme[1] in enzyme_blacklist_survival:
        survival = ''
    else:
        matchline = "\t\t\t\t<td><a href=\"/products/" + enzyme[2].split("/")[-1] + "\">" + enzyme[1] + "</a></td><td>"
        debug_print(matchline)
        survival = ''
        for line in survivaltext.splitlines():
            # If the first 87 characters of the line match
            if line[:len(matchline)] == matchline:
                debug_print("survival matchline found:\n" + line)
                re_result = re.search("(<td>)(\+*|-)(</td>)", line)
                survival = re_result.group(2)
    #
    # RETRIVE TIME-SAVER STATUS AFTER THE MAIN ENZYME PAGE HAS BEEN SCRAPED
    # https://www.neb.com/tools-and-resources/selection-charts/time-saver-qualified-restriction-enzymes
    #
    if enzyme[1] in enzyme_blacklist_timesaver:
        timesaver = ''
    else:
        matchline = "\t\t\t\t<td><a href=\"/products/" + enzyme[2].split("/")[-1] + "\">" + enzyme[1] + "</a></td><td>"
        timesaver = ''
        for line in timesavertext.splitlines():
            # If the first 87 characters of the line match
            if line[:len(matchline)] == matchline:
                debug_print("timesaver matchline found:\n" + line)
                re_result = re.search("(.gif\" alt=\"Digest in )(1*5)( minutes\" Title=\")", line)
                timesaver = re_result.group(2)

    #
    # CALCULATE FREQUENCY DATA AFTER THE MAIN ENZYME PAGE HAS BEEN SCRAPED
    # This function depends on biopython
    #

    # Read assay DNA sequences
    assayDNA_list = enumerate(SeqIO.parse("assay_DNAs.fasta", "fasta"))
    for index, seq_record in assayDNA_list:
        #print(seq_record.id)
        if seq_record.description.split()[1] == "circ.":
            not_circular = True
        else:
            not_circular = False
        if enzyme[1] in enzyme_blacklist_frequency:
            frequency = ''
        else:
            if seq_record.id == assay_DNA:
                my_seq = seq_record.seq
                stripped_enzyme = strip_HF_designation(enzyme[1])
                print(stripped_enzyme)
                rb = RestrictionBatch([stripped_enzyme])
                reldict = rb.search(my_seq, linear=False)
                frequency = len(reldict[next(iter(reldict))])

    # Print all enzyme data
    print("number of " + enzyme[1]+ "-sites in " + assay_DNA + ": " + str(frequency))
    print("buffer: " + enzyme_buffer + "\nreaction temperature: " + reaction_temperature, end="\n")
    print("url: " + url)
    print("assay DNA: " + assay_DNA)
    print("other: " + reaction_supplement + "\n" if reaction_supplement != '' else '', end="")
    print("enzyme activity: " + str(enzyme_activity))
    print("survival: " + survival)
    print("reaction temperature: " + reaction_temperature)
    print("reaction supplement: " + reaction_supplement)
    conc = ''
    enzyme_concentration.sort()
    for value in enzyme_concentration:
        conc += ',' + str(value)
    print("enzyme concentration: " + conc)
    print("timesaver: " + timesaver)

    # Update enzyme database with default buffer, assay DNA, survival after the whole text has been analyzed
    query = "UPDATE restriction_enzyme SET default_buffer = '" + enzyme_buffer + "', assay_DNA = '" + assay_DNA + "', survival = " + str(survival_list[survival]) + ", assay_DNA_cuts = " + str(frequency) + ", reaction_temperature = " + reaction_temperature[:-2] + ", reaction_supplement = '" + reaction_supplement + "', enzyme_concentration = '" + conc + "', timesaver = '" + timesaver + "' WHERE enzyme_id = " + str(enzyme[0])
    try:
        c.execute(query)
        sqlcon.commit()
        count_enzymes += 1
    except sqlcon.Error as err:
        print("Error inserting enzyme " + enzyme[1] + " to database. Error: " + str(err) + "\nQuery was: " + query)

    # Loop through buffers to add all activity data to the individual buffer tables
    for key, value in enzyme_activity.items():
            debug_print(key + ": " + str(value))
            query = "INSERT INTO `" + key + "` (enzyme_id, activity, star_activity) VALUES (" + str(enzyme[0]) + ", " + str(value[0]) + ", " + str(value[1]) + ")"
            try:
                c.execute(query)
                sqlcon.commit()
                count_buffer_entries += 1
            except sqlcon.Error as err:
                print("Error inserting enzyme " + enzyme[1] + " into buffer list. Error: " + str(err) + "\nQuery was: " + query)

    # Reset enzyme activities to "unknown" for new enzyme in the next iteration
    enzyme_activity = {'NEBuffer 1.1': [-1,-1], 'NEBuffer 2.1': [-1, -1], 'NEBuffer 3.1': [-1, -1], 'CutSmart® Buffer': [-1, -1], 'FastDigest buffer': [-1, -1], 'NEBuffer EcoRI': [-1, -1]}

print("Data for " + str(count_enzymes) + "/" + str(count_buffer_entries) + " enzymes inserted into db_ddcut database.")





sqlcon.close()
