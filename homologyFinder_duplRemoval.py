#!/usr/bin/python

#--------------------------------------------------------
# IMPORTS
#--------------------------------------------------------
import os, datetime, getpass
import sys
from sys import argv
sys.path.append(os.path.abspath(__file__).split("bitbucket")[0]+"bitbucket/aspmine/utils/") 

from aspmine_imports import *

from collections import defaultdict
import datetime
from inspect import currentframe, getframeinfo
import inspect
import itertools
import time


#--------------------------------------------------------
# ARGUMENTS, SETUPS AND PRINT
startTimet_1 = datetime.datetime.now()
startTime = datetime.datetime.now().time() # record runtime
#--------------------------------------------------------
today = datetime.date.today()
parser = CustomArgumentParser(formatter_class=SmartFormatter, usage='%(prog)s -dbname [database name]')
""" DATABASE """
parser.add_argument("-dbname", required=False, default = "aspminedb", help="Database name")
parser.add_argument("-host", required=True, help="Host name")
parser.add_argument("-user", required=True, help="User name")
parser.add_argument("-passwd", required=True, help="Password")

""" TABLES TO BE USED """
parser.add_argument("--blasttable", "-btable", required=True, type = str, default = '', help="Input BLAST table used for the paralog and homolog analysis. Ex. biblast_ID[custom]_SC[custom]. ID = min. alignment identity, SC = min. collected alignment coverage [SUMCOV = min.(q_cov + h_cov)].")
parser.add_argument("--homotable", "-htable", required=True, type = str, default = '', help="Input homolog table containg cluster/family names, species and protein names. If present, new species homologs will be generated. Ex. 'homoPF_SC[custom]_proteins'")

""" SPECIES SELECTION """
parser.add_argument("--species", "-sp", nargs = '*',  required=False, default=[None], action='store', help="List of species to analyse. The JGI species names. Do NOT have any comma or quotes.")
parser.add_argument("-all", required=False, action='store_true', help="Quits the program if the homolog table contains species that are not part of the input")

""" PARSE ARGUMENTS """
args = parser.parse_args()

dbname = args.dbname
biblastTable = args.blasttable
homoTable = args.homotable

species = args.species 
all_species = args.all

if species == [None] and not all_species:
	sys.exit("# ERROR: Please select either all species (-all) or add species to input list (-sp X Y Z)")
elif species != [None] and all_species:
	sys.exit("# ERROR: Please select only one of the flags (-all or -sp X Y Z)")

#------------------------------------------------------------------
# Check and print argument values to screen
#------------------------------------------------------------------
"""" PRINT INPUT ARGUMENTS """
print "\n#--------------------------------------------------------------"
print '# ARGUMENTS:'
print "#--------------------------------------------------------------"
print "#    Database:\t\t\t\t", dbname
print "#    Biblast table:\t\t\t", biblastTable
print "#    Homolog table:\t\t\t", homoTable



#--------------------------------------------------------
# SUBFUNCTIONS
#--------------------------------------------------------
""" CONNECT TO DATABASE """
def connect_db(args, linenumber):
	try:
		db = mdb.connect(host=args.host, user=args.user, passwd=args.passwd, db=args.dbname)
	except mdb.Error, e:
		sys.exit("# ERROR line %s - %d: %s" % (linenumber, e.args[0],e.args[1]))
	try:
		cursor = db.cursor()
	except mdb.Error, e:
		sys.exit("# ERROR line %s - %d: %s" % (linenumber, e.args[0],e.args[1]))
	return db, cursor


""" COMBINE EXECUTE AND FETCH ALL """
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


""" FETCH ALL MEMBERS FROM HFAM (homoTable) AND DELETE EXISTING ONES """
def hfamMembers_homoTable(hfam_homoTable_collected, homoTable, args):
	hfamMembers_homoTable = list()

	format_hfam = ', '.join(str(hfam_entry) for hfam_entry in hfam_homoTable_collected)	# Convert integers into string format

	# Retrieve all members in hfam
	query_hfam_data = ("SELECT * from %s WHERE hfam IN (%s);" % (homoTable, format_hfam))
	(column_names, hfam_data) = executeQuery(cursor, query_hfam_data)
	hfamMembers_homoTable = list(hfam_data)

	# Delete members with the specific hfams
	query = "DELETE from %s WHERE hfam IN (%s);" % (homoTable, format_hfam)
	(column_names, delete_data) = executeQuery(cursor, query)
	db.commit()

	return (hfamMembers_homoTable)


""" CREATE NEW COLLECTED FAMILY (hfam) """
def createNewFamily(hfamMembers, upload_counter, new_hfam, new_family_to_values2insert):

	for member in hfamMembers:
		# Append entries to existing cluster
		values = (int(new_hfam), member[1], member[2], member[3])
		new_family_to_values2insert.append(values)
		upload_counter += 1

	return(upload_counter, new_family_to_values2insert)



#--------------------------------------------------------
# FINDING MISSING ORGANISMS IN BIBLAST AND HOMOLOG TABLES
#--------------------------------------------------------
create_table = False
homoTable_orgs = list()
missing_orgs = list()
biblastTable_orgs = list()
all_possible_orgPairs = list()
diff_biblast_allPossible = list()

db, cursor = connect_db(args, inspect.stack()[0][2])

""" EXTRACT BIBLAST ORGS """
if cursor.execute("Show tables LIKE '%s'" %biblastTable):
	# Extract q_org and h_org names from biblast table
	biblast_qhorg_query = "SELECT DISTINCT q_org, h_org FROM %s" %biblastTable
	(column_names, biblastTable_qhorgs) = executeQuery(cursor, biblast_qhorg_query)

	for org_pair in biblastTable_qhorgs:
		if org_pair[0] not in biblastTable_orgs:
			biblastTable_orgs.append(org_pair[0])
		if org_pair[1] not in biblastTable_orgs:
			biblastTable_orgs.append(org_pair[1])
else:
	sys.exit("# ERROR: The biblast table does not exist. Please rerun with a new blast table.")


""" EXTRACT HOMOTABLE ORGS AND MAX(HFAM) """
if cursor.execute("SHOW TABLES LIKE '%s';" % homoTable):
	# Extract org_names from homolog table
	homo_org_query = "SELECT DISTINCT org_name FROM %s" %homoTable
	(column_names, homoTable_orgs) = executeQuery(cursor, homo_org_query)
	homoTable_orgs = map(' '.join, homoTable_orgs) 	# Convert list of tuples to list of strings

	# Find the max homolog cluster number (hfam)
	(column_names, max_hfam) = executeQuery(cursor, "SELECT MAX(hfam) FROM %s;" %homoTable)
	if max_hfam[0][0] == None:
		new_hfam = int(1)
	else:
		new_hfam = int(max_hfam[0][0]+1)
else:
	create_table = True
	new_hfam = int(1)


""" CHECK COMMANDLINE INPUT SPECIES """
if species != [None]:
	# Check if all input species are in biBLAST table
	org_not_in_biblast_orgs = list(set(species) - set(biblastTable_orgs))
	if len(org_not_in_biblast_orgs) > 0:
		sys.exit("\n# ERROR: Please rerun the scrip for generating biblast updating the species:\n%s" %org_not_in_biblast_orgs)
else:
	# IF selected all orgs from biblast table
	species = biblastTable_orgs

# Combine all input species with the orgs from the homolog table
# and create all possible orgPairs
input_homoTable_species = set(homoTable_orgs + species)
all_possible_orgPairs = tuple(itertools.product(set(input_homoTable_species), set(input_homoTable_species)))


# Check that all possible orgPairs are in biblast table
diff_biblast_allPossible = list(set(all_possible_orgPairs) - set(biblastTable_qhorgs))
if len(diff_biblast_allPossible) > 0:
	print "# ERROR: These organism pairs are missing to make a complete 'all vs. all' single linkage:"
	for dif in diff_biblast_allPossible:
		print dif
	sys.exit()

""" RETRIEVE ORG ID FROM organism TABLE """
# Creating lookup lists to be able to change commandline 
# input organisms to one common organism ID for further process
try:
	# Extract organism ID, name and real name from database
	organism_query = "SELECT org_id, name, real_name, section FROM organism"
	cursor.execute(organism_query)
	organism_name = cursor.fetchall()
except mdb.Error, e:
	sys.exit("# ERROR %d: %s" % (e.args[0],e.args[1]))
db.close

# Restructure organism data from database
orgname_to_id = dict()
for element in organism_name:
	orgname_to_id[element[1]] = element[0]


""" FIND ORGS TO BE APPENDED TO HOMOTABLE """
# Find missing q_orgs from homolog table
if homoTable_orgs:
	missing_orgs = list(set(species) - set(homoTable_orgs))
else:
	missing_orgs = species

print "#    Species in Homolog table:\t\t", homoTable_orgs
if all_species:
	print "#    List of input species :\t\tall"
if species != [None]:
	print "#    List of input species -sp:\t\t", species
print "#    Species to be appended:\t\t", missing_orgs
print "#--------------------------------------------------------------\n"

#--------------------------------------------------------
# CREATE NEW MYSQL HOMOLOG TABLE FOR SINGLE LINKAGE
#--------------------------------------------------------

if create_table:
	db, cursor = connect_db(args, inspect.stack()[0][2])
	print "# INFO: Creating homolog table %s" % homoTable
	query = ("""CREATE TABLE %s (
		`hfam` int(100) NOT NULL,
		`org_id` int(100) NOT NULL,
		`org_name` varchar(100) NOT NULL,
		`protein_id` int(100) NOT NULL,
		unique `i_hfam_org_prot` (`hfam`, `org_name`, `protein_id`),
		KEY `i_org_prot` (`org_name`, `protein_id`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1""") % homoTable
	
	executeQuery(cursor, query)
	db.commit()

	if not cursor.execute("SHOW TABLES LIKE '%s';" % homoTable):
		sys.exit("# ERROR: Table not created: %s" % homoTable)

	db.close()

#--------------------------------------------------------
# APPEND, MERGE OR CREATE NEW CLUSTERS/FAMILIES
#--------------------------------------------------------

# Limit the biBLAST search to homoTable q_orgs + runned q_orgs
searchOrgs = homoTable_orgs
startTimet_2 = ""
startTimet_3 = datetime.datetime.now()

print "# INFO: Creating homolog clusters"
db, cursor = connect_db(args, inspect.stack()[0][2])
org_count = 0
for q_org in missing_orgs:
	if org_count != 0:
		print "# INFO: Iteration time: %s" %str(datetime.datetime.now()-startTimet_2)
		print "----------------------------------------------------"

	org_count += 1

	startTimet_2 = datetime.datetime.now()
	print "\n----------------------------------------------------"
	print "# INFO: Running species: %s - %s of %s" %(q_org, org_count, len(missing_orgs))
	print "----------------------------------------------------"
	
	# Set variables 
	total_seqkey_counter = 0
	upload_counter = 0
	total_upload_counter = 0
	values2insert = []
	
	# Limit the biBLAST search to homoTable q_orgs + runned q_orgs
	searchOrgs.append(q_org)

	""" RETRIEVE ALL q_seqkeys """
	# Retrieve q_seqkey from missing q_org
	q_org_seqkey_query = "SELECT DISTINCT q_seqkey FROM %s WHERE q_org = '%s';" %(biblastTable, q_org)
	(column_names, q_org_seqkey_list) = executeQuery(cursor, q_org_seqkey_query)
	q_org_seqkey_list = list(sum(q_org_seqkey_list, ())) # list of tuples to flat list


	# Find all q_seqkey biBLAST hits (h_seqkeys) in both the input and homoTable q_orgs 
	for q_seqkey in q_org_seqkey_list:
		
		# Set and reset variables 
		total_seqkey_counter += 1
		blast_homoTable_output = list()
		hfam_homoTable_collected = list()
		hfam_values2insert_collected = list()
		collect_NULLhfam = list()

		""" RETRIEVE ALL HOMOLOGS TO q_seqkey """ 
		# Combine homoTable and biblast - Output: hfam;h_org;h_seqkey
		query_blast_homoTable = """SELECT hfam, tb.*
			FROM (
			# 2.1. Select only one column pair (here h_org, h_seqkey)
			SELECT h_org, h_seqkey
			FROM (
				# 1. Get all biblast hits to q_org/q_seqkey from db 
				SELECT q_org, q_seqkey, h_org, h_seqkey
				FROM %s
				WHERE ((q_org = '%s' AND q_seqkey = %s) OR (h_org = '%s' AND h_seqkey = %s))) ta
			# 2.2. where both q_org and h_org is in the selected orgs
			WHERE (ta.q_org IN ('%s') 
			AND ta.h_org IN ('%s'))
			) tb
			# 3.1 Join hfam from homotable where h_org/h_seqkey exists
			LEFT JOIN %s
			ON h_org = org_name AND h_seqkey = protein_id
			# 3.2 group to reduce duplicates and retrieve all potential hfams per h_org/h_seqkey
			GROUP BY hfam, h_org, h_seqkey
			;""" %(biblastTable, q_org, q_seqkey, q_org, q_seqkey, "', '".join(searchOrgs), 
				"', '".join(searchOrgs), homoTable)

		(column_names, blast_homoTable_data) = executeQuery(cursor, query_blast_homoTable) 


		# Converto output from tuples of tuples to a list of lists
		# or create empty list
		if len(blast_homoTable_data) == 0:
			blast_homoTable_output = list()
		else:
			blast_homoTable_output = blast_homoTable_data

		""" FETCH ALL HFAMS FROM homoTable """
		for blasthit in blast_homoTable_output:
			hfam_homoTable = blasthit[0]
			h_org = blasthit[1]
			h_seqkey = blasthit[2]

			# Retrieve all homoTable HFAMS
			if hfam_homoTable != None: # can hfam be collected from homologs in the query?
				if hfam_homoTable not in hfam_homoTable_collected:
					hfam_homoTable_collected.append(hfam_homoTable)
			else:
				# Collect all entries with no HFAM - Exit if not a paralog
				collect_NULLhfam.append(blasthit)


		""" FETCH ALL MEMBERS FROM HFAM AND DELETE EXISTING ONES AND CREATE NEW HFAMS"""
		members_homoTable = list()
		new_family_members = list()
		new_family_to_values2insert = list()
		collect_NULLhfam_orgID = list()
		upload_counter = len(values2insert)

		# Fetch hfam members, delete existing ones and create new uploads
		if len(hfam_homoTable_collected) > 0:
			members_homoTable = hfamMembers_homoTable(hfam_homoTable_collected, homoTable, args)

		# Include org_id to collect_NULLhfam protein members
		if len(collect_NULLhfam) > 0:
			for row_entry in collect_NULLhfam:
				collect_NULLhfam_orgID.append((row_entry[0], int(orgname_to_id[row_entry[1]]), row_entry[1], int(row_entry[2])))
		
		# Combine homoTable members with biblast hfams and hits without hfams
		new_family_members = list(set(members_homoTable + collect_NULLhfam_orgID))

		# Create new family with new hfam
		if len(new_family_members) > 0:
			(upload_counter, new_family_to_values2insert) = createNewFamily(new_family_members, upload_counter, new_hfam, new_family_to_values2insert)	

		# Remove duplicates in new_family_to_values2insert and extend values2insert
		# OBS: This might be redundant - IT IS NOT (tested)
		new_family_to_values2insert_reduced = list(set(new_family_to_values2insert))
		
		# Update values2insert and hfam count
		values2insert.extend(new_family_to_values2insert_reduced)
		new_hfam += 1
		
		""" UPLOAD TO SERVER """
		# Uploading to server
		if upload_counter >= 5000 or total_seqkey_counter == len(q_org_seqkey_list):

			total_upload_counter = total_upload_counter + upload_counter
			# Prints record number that will be inserted
			if upload_counter >= 5000 : 
				print "# INFO: Inserting record number %s" % total_upload_counter
			elif total_seqkey_counter == len(q_org_seqkey_list): 
				print "# INFO: Inserting record number %s" % total_upload_counter
				total_upload_counter = 0

			# Upload into table
			try:
				query =  "INSERT IGNORE INTO %s (hfam, org_id, org_name, protein_id) values(%s);" % (homoTable, ("%s," * len(values2insert[0])).rstrip(","))
				
				cursor.executemany(query, values2insert)
				# Add changes to database
			 	db.commit()	

				upload_counter = 0 	# restart counter
				values2insert = []	# Empty list of values
				new_family_to_values2insert = []
			except mdb.Error, e:
				print values2insert
				sys.exit( "# ERROR %s load %s %d: %s" % (homoTable, q_org, e.args[0],e.args[1] ) )

""" CLOSE DATABASE """
db.close()
new_hfam += 1

if len(missing_orgs)>0:
	print "# INFO: Iteration time: %s" %str(datetime.datetime.now()-startTimet_2)
	print "----------------------------------------------------"

print "\n----------------------------------------------------"
print "# INFO: Linking time: %s" %str(datetime.datetime.now()-startTimet_3)
print "----------------------------------------------------\n"

startTimet_4 = datetime.datetime.now()

#--------------------------------------------------------
# CREATE DUPLICATE TABLE
#--------------------------------------------------------
dupl_table = homoTable+"_dupl_delete"

""" DROP TABLE IF EXISTS """
db, cursor = connect_db(args, inspect.stack()[0][2])
executeQuery(cursor, "DROP TABLE IF EXISTS %s;" % dupl_table)
db.commit()

print "# INFO: Creating protein duplication table %s" % dupl_table
query = """CREATE TABLE %s (
	`hfam` int(100) NOT NULL,
	`name` varchar(100) NOT NULL,
	KEY `i_hfam_name` (`hfam`,`name`),
	KEY `name` (`name`)
	) ENGINE=MyISAM DEFAULT CHARSET=latin1;""" %dupl_table

executeQuery(cursor, query) 

if not cursor.execute("SHOW TABLES LIKE '%s';" % dupl_table):
	sys.exit("# ERROR: Table not created: %s" % dupl_table)

PROTdupl_table_query = """INSERT %s
	SELECT tb.hfam, CONCAT(ta.org_name, ":",ta.protein_id) name
	FROM (
	SELECT org_name, protein_id
	FROM %s
	GROUP BY org_name, protein_id
	HAVING COUNT(DISTINCT hfam) > 1) ta
	JOIN %s tb
	ON (ta.org_name = tb.org_name AND ta.protein_id = tb.protein_id)
	;""" %(dupl_table, homoTable, homoTable)

executeQuery(cursor, PROTdupl_table_query) 
executeQuery(cursor, "CREATE TABLE %s AS SELECT * FROM %s;" %(homoTable+"_dupl", dupl_table)) 

db.close()


#--------------------------------------------------------
# SINGLE-LINK ALL HFAMS
#--------------------------------------------------------
""" CHECK ROW COUNT """
db, cursor = connect_db(args, inspect.stack()[0][2])

(column_names, row_count) = executeQuery(cursor, "SELECT COUNT(*) FROM %s;" % dupl_table)
print "# INFO: There are %s duplications present in %s" % (row_count[0][0], homoTable)
print "# INFO: Mergin hfams with common org/protein pairs\n"

while row_count[0][0] != 0:
	# Initiate variables
	protein_list = list()
	hfam_list = list()
	protein_list_updated = list()
	hfam_list_updated = list()
	protein_missing = list()
	hfam_missing = list()

	# Get first duplicated protein
	(column_names, dupl_protein) = executeQuery(cursor, "SELECT name FROM %s limit 1;" % dupl_table)

	# Retrieve all hfams and org/prot pairs associated with the dupl proteins hfams
	prot_hfam_query= """SELECT * FROM %s WHERE hfam IN (SELECT hfam FROM %s WHERE name = '%s')
		GROUP BY hfam, name;""" %(dupl_table, dupl_table, dupl_protein[0][0])
	(column_names, prot_hfam) = executeQuery(cursor, prot_hfam_query)

	# Save newly found org/protein pairs and their hfams to lists
	for entry in prot_hfam:
		if entry[1] not in protein_list_updated:
			protein_list_updated.append(entry[1])
		if entry[0] not in hfam_list_updated:
			hfam_list_updated.append(entry[0])
	
	# While a new protein is added find its hfams and the duplicated proteins that are in the hfams
	while set(protein_list) != set(protein_list_updated):
		missing_prots = list(set(protein_list_updated)-set(protein_list))
		protein_list = protein_list_updated
		hfam_list = hfam_list_updated

		for element in missing_prots:
			prot_hfam_query= """SELECT * FROM %s WHERE hfam IN (SELECT hfam FROM %s WHERE name = '%s')
			GROUP BY hfam, name;""" %(dupl_table, dupl_table, element)
			(column_names, prot_hfam) = executeQuery(cursor, prot_hfam_query)

			for entry in prot_hfam:
				if entry[1] not in protein_missing:
					protein_missing.append(entry[1])
				if entry[0] not in hfam_missing:
					hfam_missing.append(entry[0])

			protein_list_updated = set.union(set(protein_list_updated), set(protein_missing))
			hfam_list_updated = set.union(set(hfam_list_updated), set(hfam_missing))
	
	protein_list = set.union(set(protein_list), set(protein_list_updated))
	hfam_list =  set.union(set(hfam_list), set(hfam_list_updated))

	#--------------------------------------------------------
	# DUPLICATE DELETION
	#--------------------------------------------------------
	""" Deletion of proteins in dupl_table """
	delete_prot_query = "DELETE FROM %s WHERE name IN ('%s');" %(dupl_table, "', '".join(protein_list))
	executeQuery(cursor, delete_prot_query)
	db.commit()

	# Retrieve all members hfams in homoTable 
	ALLprotsInHfams_query = "SELECT * FROM %s WHERE hfam in (%s) GROUP BY org_name, protein_id;" %(homoTable, ", ".join([str(i) for i in hfam_list]))
	(column_names, ALLprotsInHfams) = executeQuery(cursor, ALLprotsInHfams_query)	


	""" Deletion of hfams in homoTable """	
	delete_prot_query = "DELETE FROM %s WHERE hfam IN (%s);" %(homoTable, ", ".join([str(i) for i in hfam_list]))
	executeQuery(cursor, delete_prot_query)
	db.commit()
	
	#--------------------------------------------------------
	# CREATION OF NEW HFAMS
	#--------------------------------------------------------
	new_family_to_values2insert = list()
	
	for prot_member in list(ALLprotsInHfams):
		# Append entries to existing cluster
		values = (int(new_hfam), int(prot_member[1]), prot_member[2], int(prot_member[3]))
		new_family_to_values2insert.append(values)


	# Upload into table
	try:
		query =  "INSERT IGNORE INTO %s (hfam, org_id, org_name, protein_id) values(%s);" % (homoTable, ("%s," * len(new_family_to_values2insert[0])).rstrip(","))
		
		cursor.executemany(query, new_family_to_values2insert)
		# Add changes to database
	 	db.commit()	

		new_family_to_values2insert = []	# Empty list of values
	except mdb.Error, e:
		print new_family_to_values2insert
		sys.exit( "# ERROR %s load %s %d: %s" % (homoTable, q_org, e.args[0],e.args[1] ) )

	#--------------------------------------------------------
	# CHECKING ROW COUNTS
	#--------------------------------------------------------
	(column_names, row_count) = executeQuery(cursor, "SELECT COUNT(*) FROM %s;" % dupl_table) 
	new_hfam += 1

#--------------------------------------------------------
# DELETE EMPTY DUPLICATION TABLE
#--------------------------------------------------------
executeQuery(cursor, "DROP TABLE IF EXISTS %s;" % dupl_table)
db.commit()
db.close()

print "----------------------------------------------------"
print "# INFO: Deletion duplication time: %s" %str(datetime.datetime.now()-startTimet_4)
print "----------------------------------------------------\n"

print "\n# INFO: The program is finished"
print "# INFO: Program time: %s" %str(datetime.datetime.now()-startTimet_1)

# sys.exit("Exit at line %s\nRuntime %s\n\n\n" %(currentframe().f_lineno, str(datetime.datetime.now()-startTimet_1)))
