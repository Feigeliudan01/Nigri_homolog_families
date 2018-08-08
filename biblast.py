#!/usr/bin/python

"""----------------------------------------------------------------------------
# Code uses local functions: 
-------------------------------------------------------------------------------

CustomArgumentParser
DBconnect
executeQuery

----------------------------------------------------------------------------"""

"""----------------------------------------------------------------------------
# SUBFUNCTIONS
----------------------------------------------------------------------------"""
#------------------------------------------------------------------
# Imports
#------------------------------------------------------------------
import os, datetime, getpass
import sys
from sys import argv
sys.path.append(os.path.abspath(__file__).split("bitbucket")[0]+"bitbucket/aspmine/utils/") 

from aspmine_imports import *
import itertools


#------------------------------------------------------------------
# ARGUMENTS and setup
#------------------------------------------------------------------
startTime = datetime.now() # record runtime
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
parser.add_argument("-dbname", "-d", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-host", required=True, help="Host name")
parser.add_argument("-user", required=True, help="User name")
parser.add_argument("-passwd", required=True, help="Password")

""" TABLES """
parser.add_argument("-biblast", "-bi", required=False, type = str, default = "NULL", help="Biblast table name")
parser.add_argument("-newbiblast", "-nbi", required=False, type = str, default = "NULL", help="Biblast table name for subset")
parser.add_argument("-pident", "-id", required=False, type = int, default = 50, help="Minimum percentage alignement identity")
parser.add_argument("-sumcov", "-scov", required=False, type = int, default = 130, help="Minimum alignment coverage SUM(q_cov, h_cov)")

""" ANALYSIS """
parser.add_argument("-clean", "-c", required=False, action='store_true', help="Compile new biblast table")
parser.add_argument("-subset", "-s", required=False, action='store_true', help="Subset current biblast table")
parser.add_argument("-update", "-u", required=False, action='store_true', help="Update current biblast table")


""" SPECIES SELECTION """
parser.add_argument("-species", "-sp", nargs = '*',  required=False, default=[], action='store', help="List of species included in the table")
parser.add_argument("-section", "-sec", nargs = '*',  required=False, default=[], help="Aspergillus section selection (list)")


args = parser.parse_args()
# PRINT RUNNING DESCRIPTION 
now = datetime.now()
print '# ' + ' '.join(argv)
print '# ' + now.strftime("%a %b %d %Y %H:%M")
print '# USER: ' + getpass.getuser()
print '# CWD: ' + os.getcwd()
if os.name != "nt":
	print '# HOST: ' + ' '.join([ os.uname()[0] , os.uname()[4] , 
       	                         os.uname()[2] , os.uname()[1] ])	


""" PARSE ARGUMENTS """
dbname = args.dbname
biblast = args.biblast
newbiblast = args.newbiblast
pident = args.pident
sumcov = args.sumcov

clean = args.clean
subset = args.subset
update = args.update

species = args.species
section = args.section

loadstamp = time.strftime("%c")
#orgcheck = args.orgcheck

#------------------------------------------------------------------
# Print argument values to screen
#------------------------------------------------------------------
print "\n#--------------------------------------------------------------\n"
print "ARGUMENTS:"
print "# Database\t: ", dbname
print "# Current biblast table\t: ", biblast

if not clean and not subset and not update:
	sys.exit("\n# ERROR: Please select one of the analysis: clean (-c), subset(-s) or update (-u)")
if clean and subset:
	sys.exit("\n# ERROR: Please ONLY select one of the analysis: clean (-c), subset(-s) or update (-u)")
if clean and update:
	sys.exit("\n# ERROR: Please ONLY select one of the analysis: clean (-c), subset(-s) or update (-u)")
if update and subset:
	sys.exit("\n# ERROR: Please ONLY select one of the analysis: clean (-c), subset(-s) or update (-u)")

if clean:
	if newbiblast != "NULL":
		print "# Compile new biblast table\t: ", newbiblast
		print "\n# SUMCOV data information:"
		print "# Minimum percentage alignement identity\t: ", args.pident
		print "# Minimum alignment coverage SUM(q_cov, h_cov)\t: ", args.sumcov
		if species == None and section == None:
			print "# List of species\t: All"
			print "# List of sections\t: ", section
			species = "All"
		else:
			print "# List of species\t: ", species
			print "# List of sections\t: ", section
	else:
		sys.exit("\n# ERROR: Please add new blast table name (-nbi) when compiling new blast table")

if subset:
	if newbiblast != "NULL":
		if species == None and section == None:
			sys.exit("\n# ERROR: Please add species or section to be subsetted (-sp or -sec)")
		else:
			print "# Subset to new biblast table\t: ", newbiblast
			print "\n# BLAST selection criteria:"
			print "# Minimum percentage alignement identity\t: ", args.pident
			print "# Minimum alignment coverage SUM(q_cov, h_cov)\t: ", args.sumcov
			print "# List of species\t: ", species
			print "# List of sections\t: ", section
	else:
		sys.exit("\n# ERROR: Please add new blast table name (-nbi) when subsetting to new blast table")

if update:
	if biblast != "NULL":
		print "# Update biblast table\t: ", biblast
		print "\n# SUMCOV data information:"
		print "# Minimum percentage alignement identity\t: ", args.pident
		print "# Minimum alignment coverage SUM(q_cov, h_cov)\t: ", args.sumcov
		if species == None and section == None:
			print "# List of species\t: All"
			print "# List of sections\t: ", section
			species = "All"
		else:
			print "# List of species\t: ", species
			print "# List of sections\t: ", section
	else:
		sys.exit("\n# ERROR: Please add blast table name (-bi) when updating the table")

print "\n#--------------------------------------------------------------\n" 


'''============================================================================
# SUB PROGRAMS
============================================================================'''

#------------------------------------------------------------------
# Connection to database
#------------------------------------------------------------------
def connect_db(args):
	try:
		db = mdb.connect(host=args.host, user=args.user, passwd=args.passwd, db=args.dbname)
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return db, cursor

#------------------------------------------------------------------
# Combine execute and fetch all into one function
#------------------------------------------------------------------
def executeQuery(cursor, query):
	(columns, result) = ([],[])
	try:
		cursor.execute(query)
		result = cursor.fetchall()
		if result:
			if len(result) > 0:
				columns = map(lambda x:x[0], cursor.description) 	

	except mdb.Error, e:
		sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
	return columns, result	


def Find_all_possible_pairs(species, input_orgPair, compare_orgPair):
	# Fetch q_orgs from table_orgs
	q_orgs = [x[0] for x in input_orgPair]
	q_orgs = list(set(q_orgs))

	# Combine table q_org with input species
	species = list(set(species))
	if len(species) > 0:
		collected_species = list(set(species.extend(q_orgs)))
	else:
		collected_species = q_orgs

	# Create all possible q_org, h_org pairs from input species and table_orgs
	possible_collected_pairs = tuple(itertools.product(collected_species, collected_species))

	# Check if possible_collected_pairs are in the table_orgs:
	pair_found_in_table = list()
	pair_NOT_found_in_table = list()
	for pair in possible_collected_pairs:
		if pair in compare_orgPair:
			pair_found_in_table.append(pair)
		else:
			pair_NOT_found_in_table.append(pair)
	
	return(pair_found_in_table, pair_NOT_found_in_table)


'''============================================================================
# MAIN PROGRAM
============================================================================'''
#------------------------------------------------------------------
# Extract and compare biblast q/h_orgs
#------------------------------------------------------------------
""" FETCH BLAST ORGS (q_org and h_org) """
print "# INFO: Retrieving organism pairs from blast"
db, cursor = connect_db(args)
(column_names, blast_orgPair) = executeQuery(cursor, "SELECT DISTINCT q_org, h_org FROM blast")
db.close()

biDir_blast_orgPair = list()
oneDir_blast_orgPair = list()

# Collect bidirectional orgpairs from blast
for orgPair in blast_orgPair:
	opposite_orgPair = (orgPair[1], orgPair[0])
	if opposite_orgPair in blast_orgPair:
		biDir_blast_orgPair.append(orgPair)
	else:
		oneDir_blast_orgPair.append(orgPair)

# Find all possible pairs from blast and species input 
# and devide into those present and not present in blast
pair_found_in_blast = list()
pair_NOT_found_in_blast = list()
(pair_found_in_blast, pair_NOT_found_in_blast) = Find_all_possible_pairs(species, blast_orgPair, blast_orgPair)


""" CHECK INPUT SPECIES ARE PRESENT IN blast """
if len(species) > 0:
	q_orgs = [x[0] for x in blast_orgPair]
	q_orgs = list(set(q_orgs))

	species_not_in_blast = list()
	for q_species in species:
		if q_species not in q_orgs:
			species_not_in_blast.append(species)

	if len(species_not_in_blast):	
		print "# ERROR: These species are not in the blast table"
		for nb_species in species_not_in_blast:
			print "ERROR: ", nb_species
		sys.exit()

#------------------------------------------------------------------
# RUN SELECTED ANALYSIS
#------------------------------------------------------------------
input_orgPair = list()
compare_orgPair = list()
already_created_orgPair = list()

if update:
	""" FETCH BIBLAST ORGS (q_org and h_org) """
	db, cursor = connect_db(args)
	if cursor.execute("Show tables LIKE '%s'" %biblast):
		print "# INFO: Retrieving organism pairs from %s" %biblast
		(column_names, biblast_orgPair) = executeQuery(cursor, "SELECT DISTINCT q_org, h_org FROM %s" %biblast)
	else:
		db.close()
		sys.exit("# ERROR: The BLAST table does not exist. Please rerun with a new blast table.")
	db.close()

	input_orgPair = biDir_blast_orgPair
	compare_orgPair = biDir_blast_orgPair
	already_created_orgPair = biblast_orgPair


if clean:
	""" DROP AND CREATE NEW BIBLAST TABLE """
	db, cursor = connect_db(args)
	cursor.execute("DROP TABLE IF EXISTS %s " % biblast)
	db.close()

	input_orgPair = blast_orgPair
	compare_orgPair = biDir_blast_orgPair

if len(input_orgPair) == 0 or len(compare_orgPair) == 0:
	sys.exit("# DEBUG ERROR: input_orgPair: %s or compare_orgPair: %s empty" %(input_orgPair, compare_orgPair))


if clean and len(species) > 0:
	""" CREATE NEW TABLE IF CLEAN CHOSEN """
	
	print "# INFO: Creating - ", biblast
	query = """CREATE TABLE %s (
	`q_org` varchar(100) NOT NULL,
	`q_seqkey` varchar(100) NOT NULL,
	`h_org` varchar(100) NOT NULL,
	`h_seqkey` varchar(100) NOT NULL,
	`pident` decimal(10,2) DEFAULT NULL,
	`q_cov` decimal(10,2) DEFAULT NULL,
	`h_cov` decimal(10,2) DEFAULT NULL,
	KEY `i_qh_org_seqkey` (`q_org`,`q_seqkey`,`h_org`,`h_seqkey`),
	KEY `i_qh_org` (`q_org`,`h_org`),
	KEY `i_horg_seqkey` (`h_org`,`h_seqkey`)
	) ENGINE=InnoDB DEFAULT CHARSET=latin1""" %biblast

	db, cursor = connect_db(args)
	executeQuery(cursor, query)

	if not cursor.execute("SHOW TABLES LIKE '%s';" % biblast):
		sys.exit("# ERROR: Table not created: %s" % biblast)

	db.close()



if update or (clean and len(species) > 0):
	
	""" SELECT ORG PAIRS FROM BIBLAST/BLAST AND INPUT SPECIES """
	pair_found_in_compareTable = list()
	pair_NOT_found_in_compareTable = list()
	(pair_found_in_compareTable, pair_NOT_found_in_compareTable) = Find_all_possible_pairs(species, input_orgPair, compare_orgPair)

	""" COMPARE BIBLAST TO SELECTED POSSIBLE ORG PAIRS FOUND IN BLAST TABLE"""
	missing_orgPair = list()
	missing_orgPair = list(set(set(pair_found_in_compareTable) - set(already_created_orgPair)))


	""" INSERT MISSING ORG PAIRS """
	if len(missing_orgPair) == 0:
		print "\n# INFO: All bidirectional blast organism pair are already updated"
	else:
		print "\n# INFO: Inserting %s organism pairs - %s" %(len(missing_orgPair), biblast)
		org_count = 0
		db, cursor = connect_db(args)
	 	for orgPair in missing_orgPair:
	 		org_count += 1
	 		print "# INFO: %s of %s: %s" %(org_count, len(missing_orgPair), orgPair)
			q_org = orgPair[0]
			h_org = orgPair[1]

			startTime_0 = datetime.now() # record runtime

			bi_query = """INSERT INTO %s 
			SELECT tb.q_org, tb.q_seqkey, tb.h_org, tb.h_seqkey, tb.pident, tb.q_cov, tb.h_cov 
				FROM (
					SELECT q_org, q_seqkey, h_org, h_seqkey, pident, q_cov, h_cov 
					FROM blast 
					WHERE (q_org = "%s" AND h_org = "%s" AND pident >= %s AND (h_cov + q_cov) >= %s)) tb
				JOIN (
					SELECT q_org, q_seqkey, h_org, h_seqkey, pident, q_cov, h_cov 
					FROM blast 
					WHERE (q_org = "%s" AND h_org = "%s" AND pident >= %s AND (h_cov + q_cov) >= %s)) td
				ON (tb.q_org=td.h_org AND tb.q_seqkey=td.h_seqkey AND tb.h_org=td.q_org AND tb.h_seqkey=td.q_seqkey)
			;""" % (biblast, q_org, h_org, pident, sumcov, h_org, q_org, pident, sumcov)


			cursor.execute(bi_query)
			db.commit()	# OBS: NOT SURE I SHOULD DO THIS
			print "# INFO: Insert runtime %s: %s" %(str(datetime.now()-startTime_0), orgPair)
			
		db.close()

		if len(pair_NOT_found_in_compareTable) > 0 and update:
			pair_not_in_blast = list(set(set(pair_NOT_found_in_compareTable)-set(oneDir_blast_orgPair)))
			print "\n# INFO: These org_pairs are not found in the blast table to complete the biblast:"
			print "# INFO: (q_org, h_org)"
			for pair_noBlast in pair_not_in_blast:
				print "# INFO: %s" %str(pair_noBlast)


if clean and len(species) == 0:
	print "# INFO: Creating - ", biblast
	print "# INFO: Including all reciprocal organism pairs found in blast"
	# Create 
	startTime_1 = datetime.now() # record runtime
	create_query = "CREATE TABLE %s AS\
		SELECT ta.q_org, ta.q_seqkey, ta.h_org, ta.h_seqkey, ta.pident, ta.q_cov, ta.h_cov FROM blast AS ta join blast AS tb ON \
			(ta.q_seqkey=tb.h_seqkey AND ta.h_seqkey=tb.q_seqkey AND ta.q_org=tb.h_org AND ta.h_org=tb.q_org ) \
		WHERE (ta.pident >= %s AND ((ta.h_cov + ta.q_cov) >= %s) AND tb.pident >= %s AND ((tb.h_cov + tb.q_cov) >= %s));" % (biblast, pident, sumcov, pident, sumcov)

	db, cursor = connect_db(args)
	executeQuery(cursor, create_query)

	if not cursor.execute("SHOW TABLES LIKE '%s';" % biblast):
		sys.exit("# ERROR: Table not created: %s" % biblast)

	print "# INFO: Create table runtime %s: %s" %(str(datetime.now()-startTime_1), biblast)

	print "# INFO: Createing index"
	startTime_2 = datetime.now() # record runtime
	executeQuery(cursor, "CREATE INDEX i_qh_org_seqkey  ON %s (q_org, q_seqkey, h_org, h_seqkey);" %biblast)
	executeQuery(cursor, "CREATE INDEX i_qh_org ON %s (q_org,h_org);" %biblast)
	print "# INFO: Create index runtime %s: %s" %(str(datetime.now()-startTime_2), biblast)
	db.close()


# if len(pair_NOT_found_in_blast) > 0:
# 	print "\n# INFO: These organism pairs are not found in the blast table:"
# 	print "# INFO: (q_org, h_org)"
# 	for pair_oneDirBlast in oneDir_blast_orgPair:
# 		print "# INFO: %s" %str(pair_oneDirBlast) # Should be opposite - and is redudant

if update:
	print "\n# INFO: Checking organism pairs from the updated biblast"
	db, cursor = connect_db(args)
	(column_names, biblast_orgPair_new) = executeQuery(cursor, "SELECT DISTINCT q_org, h_org FROM %s" %biblast)
	db.close()

	not_created_pairs = list()

	for created_pair in missing_orgPair:
		if created_pair not in biblast_orgPair_new:
			not_created_pairs.append(created_pair)

	if len(not_created_pairs) > 0:
		print "# WARNING: %s organism pairs did not create entries in the updated %s table:" %(len(not_created_pairs), biblast)
		print "# WARNING: (q_org, h_org)"
		for no_pair in not_created_pairs:
			print "# WARNING: %s" %str(no_pair)
		print "# WARNING: Please check the q_seqkey/h_seqkey of the organisms, they might be uploaded incorrectly\n"
	else:
		print "# INFO: All organism pairs were inserted correctly"










