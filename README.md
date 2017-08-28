# reoptimize
Python tool to optimize reaction conditions for restriction enzyme digests of DNA

This tools lets you search for a suitable buffer for a simultaneous digest of DNA
using  an arbitrary number of enzymes. It calculates the amount of enzymes
needed (in units) and takes into account the duration of the restriction digest
(based on the NEB data on enzyme survival in a reaction). It's a very crude first
attempt (and crashes immediately when you submit some wrongly formated data, inclu-
ding capitalization), but it's going to improve slowly according to how much time
we have. Planned are also GUI interfaces for Linux, macOS and Windows.

## Input needed:

1. Enzyme name(s)
2. How many cuts does the target DNA have for each selected enzyme?
3. How long is the target DNA (base pairs)?
4. How long are you going to incubate the reaction (hours)?
5. How much DNA do you want to cut (microgram)?

## Output:

1. Possible buffers in the order of suitability (or the result "simultaneous digest not recommended")
2. Amount of each enzyme needed.


## Requirements:
Biopython, sqlite3, click

Since even the latest Biopython distribution doesn't contain all enzymes sold by NEB,
you need to manually update the Restriction_Dictionary.py manually with the file
Restriction_Dictionary.py by copying it into the Bio/Restriction folder of the
folder, where your python3 stores the python packages. On Ubuntu 16.04, this
would be /usr/lib/python3/dist-packages/Bio/Restriction.

## Files:

*reoptimize.py*
The script that does the calculations. Usage examples:

This is a double digest with AflIII and HindIII, where the target DNA
has two AflIII sites and one HindIII site:  

>./reoptimize.py digest -e 'AflIII 2' -e 'HindIII 1'


This gives all necessary parameters via the command line:

>./reoptimize.py digest -e 'EcoRI 2' -e 'HindIII 3' -l 3000 -t 4 -m 2

-l (length of target dna, in base pairs)
-t (incubation time, in hours)
-m (amount of DNA, in micrograms)

If enzymes or other parameters are omitted, the program will prompt
for them!


*make_sqlite_database.py*
This script fetches all the data for NEB enzymes from the NEB web pages and
assembles the database that is needed for the script to run. Running it
results in the dadabase file "REsqlite3.db"

*assay_DNAs.fasta*
This files contains the full DNA sequences of all assay DNAs used by NEB. We
include it here to avoid querying the "Frequency of restriction sites" table
(which is anyway incomplete).
