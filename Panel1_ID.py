'''
Script to identify great vs scalloped/carolina hammerheads
Last updated: 8/23/17
A. Barker
'''

import csv
import pickle

#Open dictionary
with open('Panel1_dict.pickle', 'rb') as data:
    Great_ScallCar_dict = pickle.load(data)

# Import allele data: text file imported as a list, each individual is its own list
with open('Library_data_files/SLH3/SLH3.ScallCar_Great_SPID.formatted.txt', 'rU') as f:
    indiv_data = list(csv.reader(f, delimiter = '\t'))
# Get contig names form first line, the delete from list
contig_names = indiv_data[0]
indiv_data.remove(indiv_data[0])

num_individuals = len(indiv_data)

print indiv_data[0]
species_results = []
great_idv = []
scallCar_idv = []
ind_check = []
unknown_idv = []

scalCar_total = 0
great_total = 0

# create text files to get detailed species ID breakfdown.
unknown = open('Results/SLH3/SLH3.unknown_scallcar_great_details.txt', 'w')
scallcar = open('Results/SLH3/SLH3.scallcar_details.txt', 'w')
great = open('Results/SLH3/SLH3.Great_details.txt', 'w')

# identify individuals
for ind in range(0, len(indiv_data)):
    ind_id = []
    scallCar_allele = 0
    great_allele = 0
    for c in range(0, len(contig_names)):
        contig = contig_names[c]
        pos = c + 1 # add 1 because first position is individual name
        alleles = str(indiv_data[ind][pos]) # alleles for each contig are together (130130) and need separated (130, 130)
        allele1 = float(alleles[0:3])
        allele2 = float(alleles[3:6])
        # check each allele
        if allele1 != 0.0:
            try:
                ind_id.append(Great_ScallCar_dict[contig][allele1])
            except KeyError:
                continue
        if allele2 != 0.0:
            try:
                ind_id.append(Great_ScallCar_dict[contig][allele2])
            except KeyError:
                continue
    for i in ind_id:
        if i == 'scalloped/carolina':
            scallCar_allele += 1
        if i == 'great':
            great_allele += 1

    total_alleles = len(ind_id)
    scallCar_per = (float(scallCar_allele)/total_alleles)
    great_per = (float(great_allele)/total_alleles)
    if scallCar_per >= 0.95:
        scalCar_total += 1
        scallCar_idv.append(indiv_data[ind][0]) #to get name of individual
        scallcar.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped/carolina' + '\t' + str(
            scallCar_allele) + '\t' + str(scallCar_per) + '\n' + 'great' + '\t' + str(
            great_allele) + '\t' + str(great_per) + '\n' + '-------------------------------' + '\n')
    elif great_per >= 0.95:
        great_total += 1
        great_idv.append(indiv_data[ind][0])
        great.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped/carolina' + '\t' + str(
            scallCar_allele) + '\t' + str(scallCar_per) + '\n' + 'great' + '\t' + str(
            great_allele) + '\t' + str(great_per) + '\n' + '-------------------------------' + '\n')
    else:
        unknown_idv.append(indiv_data[ind][0])
        unknown.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped/carolina' + '\t' + str(
            scallCar_allele) + '\t' + str(scallCar_per) + '\n' + 'great' + '\t' + str(
            great_allele) + '\t' + str(great_per) + '\n' + '-------------------------------' + '\n')
    print '-------------------------------'
    print indiv_data[ind][0]
    print total_alleles
    print 'scalloped/carolina' '\t', scallCar_allele, '\t', scallCar_per
    print 'great' '\t', great_allele, '\t', great_per

unknown.close()
scallcar.close()
great.close()

print scalCar_total
print great_total
print scallCar_idv
print great_idv


# export list of ID'd individuals into respective species files
scalloped_carolina = open('Results/SLH3/SLH3.scallCar_ID.txt', 'w')
for sp in scallCar_idv:
    scalloped_carolina.write(str(sp) + '\n')
scalloped_carolina.close()
great = open('Results/SLH3/SLH3.Great_ID.txt', 'w')
for sp in great_idv:
    great.write(str(sp) + '\n')
great.close()
unknown = open('Results/SLH3/SLH3.unknown_scalcar_gr_ID.txt', 'w')
for sp in unknown_idv:
    unknown.write(str(sp) + '\n')
unknown.close()
