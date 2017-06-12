# BIMM-185-Week-7

import re, numpy as np, operator
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from operator import itemgetter


TU_table_open = open("TU_table")
TU_table= TU_table_open.readlines()

# Columns:
# (1) Transcription Unit identifier assigned by RegulonDB
# (2) Transcription unit name
# (3) Operon name containing the transcription unit
# (4) Name of the gene(s) contained in the transcription unit
# (5) Promoter Name
# (6) Evidence that supports the existence of the transcription unit
# (7) Evidence confidence level (Confirmed, Strong, Weak)

#Tiffany Gassmann
#BIMM 185 Week 6

#-------------------------------OPEN-FILES------------------------------------------------------------------------------

operonset_table_open =open("operonset_table")
operonset_table = operonset_table_open.readlines()

# Columns:
# (1) Operon name
# (2) First gene-position left
# (3) Last gene-position right
# (4) DNA strand where the operon is coded
# (5) Number of genes contained in the operon
# (6) Name or Blattner number of the gene(s) contained in the operon
# (7) Evidence that support the existence of the operon's TUs
# (8) Evidence confidence level (Confirmed, Strong, Weak)

geneproduct_table_open = open("geneproduct_table")
geneproduct_table = geneproduct_table_open.readlines()
# Columns:
# (1) Gene identifier assigned by RegulonDB
# (2) Gene name
# (3) Blattner number (bnumber) of the gene
# (4) Gene left end position in the genome
# (5) Gene right end position in the genome
# (6) DNA strand where the gene is coded
# (7) Product name of the gene
# (8) Evidence that supports the existence of the gene
# (9) PMIDs list


#Reads File in and seperates each according to a tabular
def readfile_table(table):
    for i in xrange(len(table)):
        #making table into list of lists without tab
        table[i] = table[i].strip().split("\t")
    return table


def get_locus_tag(table):
    for i in xrange(len(table)):
        # making table into list of lists without tab
        table[i] = table[i].strip().split("\t")

    #B-NUMS
    B_NUM = 2
    b_nums = ""

    for row in table:
        try:
            if row[1]:
                index = len(row) -2
                if row[1]:
                    b_nums += str(row[1]).strip(" ") + ", "
        except IndexError:
            "None"

    return b_nums



open_pos = open("pos_con")
pos_con_table = open_pos.readlines()




#Formal File
def format (list1):
    return str(list1).replace('[','').replace(']','').replace("'",'').replace('(','').replace(')','')


#---------------------------------------------------------------------------------------------------------------------

#EXTRACT GENES IN CURATED OPERONS
#Task: from TU_table and operonset_table extract the name of the genes in each transcription unit
# with evidence "strong" or "confirmed"

#Note: if there in genes in file TU_table and Operonset_table use the file GeneProduct_table to map
# the gene name to the locus_tag which we have in our data base

#Note: we can transform all operons to their b-numbers (locus-tag) before we query the database for
# their coordinates

#Extract gene name of each transcription unit with evidence strong or confirmed
#Returns: list of gene names
def extract_genes_TU_table():
    tu_table = readfile_table(TU_table)

    #index of name in TU_file
    NAME = 3
    find_confirmed_strong = "?:Strong|Confirmed"

    gene_names = []

    for line in tu_table:
        CONFIDENCE_LEVEL = len(line)-1
        if  re.findall(line[CONFIDENCE_LEVEL],find_confirmed_strong):
            if line[NAME]:
                gene_names.append(line[NAME])
    return gene_names
#----------------------------------------------------------------------------------------------------------------------

def name_to_num():

    #Number Index
    B_NUM = 2

    #Name Index
    G_NAME = 1

    geneprod_table = readfile_table(geneproduct_table)

    name_2_num_dict = {}

    for line in geneprod_table:
        if line[B_NUM] and line[G_NAME]:
            name_2_num_dict[line[G_NAME]] = line[B_NUM]

    return name_2_num_dict


#----------------------------------------------------------------------------------------------------------------------

#Task: Extracting the data to model h1 from the tables input files
#The distances of genes inside operons will be taken directly from all
# retrieved operons with two or more genes. ** Operons with 2+ genes ONLY **
def operon_model_h1():

    genename_2_bnum = name_to_num()
    operon_table = readfile_table(operonset_table)

    # index of operon table
    OPERON_NAME = 0
    GENE_NAME = 5

    find_confirmed_strong = "?:Strong|Confirmed"

    operon_dict = {}
    gene_names = []
    operon_names = []

    for line in operon_table:
        CONFIDENCE_LEVEL = len(line) - 1
        if re.findall(line[CONFIDENCE_LEVEL], find_confirmed_strong):
            #more than 2 genes - extract
            number_genes = len(line[GENE_NAME].split(","))
            if number_genes >= 2:
                try:
                    if str(line[CONFIDENCE_LEVEL]).isalpha():
                        names = str(line[GENE_NAME])
                        operon_names.append(line[OPERON_NAME])
                        gene_names.append(line[GENE_NAME])
                        b_nums = []
                        for name in names.split(','):
                            b_num = (genename_2_bnum.get(name))
                            b_nums.append(b_num)

                        #print "\t".join([line[OPERON_NAME], format(b_nums), line[CONFIDENCE_LEVEL]])
                except IndexError:
                    return None

                operon_dict[line[OPERON_NAME]] = b_nums
    #retuns dictionary mapping operons to each b_num for calculating distances within operons
    return operon_dict




#--------------------------------------------------------------------------------------------
#Positive Control
#Tasks:
# 1.) Read in Output of SQL Query Format: b_num gene_id left    right   string'").strip("(").strip("')").strip("'(").strip("")
# 2.) Load into dictionary with b_num as key
# 3.) Calculate the size of each operon according to the size of genes contained. ** Operons with two of more genes **
def map_positive_control(operon_dict):
    #Operon Dict : [Operon Name: b_num of gene...]
    gene_length_dict = {}
    #Positive Query MySqlDB resultsa
    open_query_res = open("pos_control_res")
    pos_query = open_query_res.read()

    replace_Long_char = pos_query.replace("L","")
    #list of tuples containing B num, gene_id, left, right, sztrand
    positive_control = re.split(r',\s*(?![^()]*\))', replace_Long_char)

    BNUM = 0
    LEFT = 2
    RIGHT = 3
    operon_distances = []
    distance = 0
    #iterates over all possible gene mapped to their left and right coordinate
    for (item) in positive_control:
        tup = (item.split(","))
        str(tup)
        #access eaach item in input file
        left = str(tup[LEFT])
        right = str(tup[RIGHT])
        b_num = str(tup[BNUM].replace("'",'').replace("(",""))
        #make dictionary mappign the lenght of each genes
        gene_length_dict [b_num] = int(right) - int(left) + 1
    #iterate over each gene in operon dictionry mapping operon: gene
    for genes in operon_dict.values():
        for gene in genes:
            try:
                #distance is the sumation of all gene lengths in each operon
                distance += int(gene_length_dict[gene])
                operon_distances.append(distance)
            except KeyError:
                continue
        distance = 0
    #Returns a list of all distances of each operon
    for dit in operon_distances:
        print dit


operon_dict =  operon_model_h1()
map_positive_control(operon_dict)

#----------------------------NEGATIVE CONTROL---------------------------------------------------------------------------

#Task: Extracting the data to model h1 from the tables input files
#The distances of genes inside operons will be taken directly from all
#retrieved operons. ** Operons with all gene numbers**
def operon_model_h0():
    genename_2_bnum = name_to_num()

    operon_table = readfile_table(operonset_table)


    # index of operon table
    OPERON_NAME = 0
    GENE_NAME = 5
    #Only genes with Strong or Confirmed evidence
    find_confirmed_strong = "?:Strong|Confirmed"
    neg_operon_dict = {}

    #Iterate over whole operon table file -- all operons
    for line in operon_table:
            CONFIDENCE_LEVEL = len(line) - 1
            #find all confirmed or strong evidence ONLY
            if re.findall(line[CONFIDENCE_LEVEL], find_confirmed_strong):
                # more than 1 genes - extract
                number_genes = len(line[GENE_NAME].split(","))
                if number_genes >= 1:
                    try:
                        if str(line[CONFIDENCE_LEVEL]).isalpha():
                            names = str(line[GENE_NAME])
                            b_nums = []
                            for name in names.split(','):
                                b_num = (genename_2_bnum.get(name))
                                b_nums.append(b_num)

                                # print "\t".join([line[OPERON_NAME], format(b_nums), line[CONFIDENCE_LEVEL]])
                    except IndexError:
                        return None

                    neg_operon_dict[line[OPERON_NAME]] = b_nums
    # Returns Dictionary with genes in ALL operons
    return neg_operon_dict


def map_negative_control(neg_operon_dict):
    #NEGATIVE CONTROL: INDEX  B NUM  GENE ID  LEFT   RIGHT STRAND
    open_query_res = open("neg_control_res")
    neg_query = open_query_res.read()

    replace_Long_char = neg_query.replace("L", "")
    # list of tuples containing index, B num, gene_id, left, right, strand
    negative_control = re.split(r',\s*(?![^()]*\))', replace_Long_char)

    INDEX = 0
    B_NUM = 1
    LEFT = 3
    RIGHT = 4
    STRAND = 5

    gene_length_strand_dict = {}

    forward = "+"
    reverse = "-"

    distances = []

    #make a dictionary containing alll genes needed to be mapped with the index of gene, size of  gene, strand
    for item in negative_control:
        tup = item.split(",")

        # access eaach item in input file
        left = str(tup[LEFT])
        right = str(tup[RIGHT])
        direction = str(tup[STRAND]).replace(")", '').replace("(", "").replace("'","").replace(" ",'')
        index = str(tup[INDEX].replace("(",""))
        b_num = str(tup[B_NUM].replace("'", '').replace("(", "").replace(" ",""))
        # make dictionary mappign the lenght of each genes and strand direction
        # dictionary of tuples (index, size_gene, direction)
        gene_length_strand_dict[b_num] = (int(index), int(right) - int(left) + 1,direction)

    #use the gene/len/index/strand dictionary to map the spaces between operons in neg control list

    print "neg operons", neg_operon_dict
    print "gene_length_strand", gene_length_strand_dict

    operon_indexes = []
    sorted_indexes = []

    for genes in neg_operon_dict.values():
        try:
            #print 'genes', genes
            for gene in genes:
                the_gene = gene_length_strand_dict[gene]
                #print "curr_gene",curr_gene
                curr_gene = gene
                curr_index = the_gene[0]
                curr_direction = the_gene[2]

                #print "operon_index",operon_indexes
                #print "curr_ind", curr_index

                sorted_indexes.append((int(curr_index), curr_gene, curr_direction))
        except KeyError:
            continue
    #list of all genes in all operons sorted by the index
    all_sorted = sorted(sorted_indexes, key = lambda x:x[0])
    print all_sorted

    forward_index = 0
    reverse_index = 0

    curr_key = "b001"
    for b_num in neg_operon_dict.keys():
        if b_num == curr_key:
            forward_index = all_sorted.f








        #looking for items going in the same direction





    return distances

    # access eaach item in input file
    left = str(tup[LEFT])
    right = str(tup[RIGHT])
    b_num = str(tup[BNUM].replace("'", '').replace("(", ""))
    # make dictionary mappign the lenght of each genes
    gene_length_dict[b_num] = int(right) - int(left) + 1

#neg_operon_dict = operon_model_h0()
#print map_negative_control(neg_operon_dict)






#----------------------------------------KDE-------------------------------------------------


#Steps in Creating KDE
#1.) generate sample gaussin_kde class to feed data in
#2.) create index for our function
#3.) Plot Histogram of data

def kde():
    distances_h1 = map_positive_control()
    distances_h0 = map_negative_control()

    mu, sigma = 100, 15

    # the histogram of the data
    n, bins, patches = plt.hist(distances_h1, 50, normed=1, facecolor='green', alpha=0.50)

    # add a 'best fit' line
    y = mlab.normpdf(bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)

    plt.xlabel('Intergenic Distance')
    plt.ylabel('Frequency')


    plt.show()

#kde()







