'''
Script to identify scalloped vs carolina hammerheads
Note: Individuals previously ID'd as great hammerheads must be removed from input file prior to running this script
Last updated: 8/31/17
A. Barker
'''


import csv
import pickle

#Open dictionary
with open('Panel2_dict.pickle', 'rb') as data:
    Scall_Car_dict = pickle.load(data)

# Import allele data: text file imported as a list, each individual is its own list
with open('Library_data_files/SLH3/SLH3.Scall_Car_SPID_EDITED.formatted.txt', 'rU') as f:
    indiv_data = list(csv.reader(f, delimiter = '\t'))
# Get contig names form first line, the delete from list
contig_names = indiv_data[0]
indiv_data.remove(indiv_data[0])

num_individuals = len(indiv_data)

species_results = []
scalloped_idv = []
carolina_idv = []
ind_check = []
unknown_idv = []

scalloped_total = 0
carolina_total = 0

# create text files to get detailed species ID breakdown.
unknown = open('Results/SLH3/SLH3.unknown_scall_car_details.txt', 'w')
scalloped = open('Results/SLH3/SLH3.scalloped_details.txt', 'w')
carolina = open('Results/SLH3/SLH3.carolina_details.txt', 'w')

# identify individuals
for ind in range(0, len(indiv_data)):
    ind_id = []
    scalloped_allele = 0
    carolina_allele = 0
    for c in range(0, len(contig_names)):
        contig = contig_names[c]
        pos = c + 1 # add 1 because first position is individual name
        alleles = str(indiv_data[ind][pos]) # alleles for each contig are together (130130) and need separated (130, 130)
        allele1 = float(alleles[0:3])
        allele2 = float(alleles[3:6])
        # check each allele
        if allele1 != 0.0:
            try:
                ind_id.append(Scall_Car_dict[contig][allele1])
            except KeyError:
                continue
        if allele2 != 0.0:
            try:
                ind_id.append(Scall_Car_dict[contig][allele2])
            except KeyError:
                continue
    for i in ind_id:
        if i == 'scalloped':
            scalloped_allele += 1
        if i == 'carolina':
            carolina_allele += 1
    total_alleles = len(ind_id)
    scalloped_per = (float(scalloped_allele)/total_alleles)
    carolina_per = (float(carolina_allele)/total_alleles)
    if scalloped_per >= 0.95:
        scalloped_total += 1
        scalloped.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped' + '\t' + str(
            scalloped_allele) + '\t' + str(scalloped_per) + '\n' + 'carolina' + '\t' + str(
            carolina_allele) + '\t' + str(carolina_per) + '\n' + '-------------------------------' + '\n')
        scalloped_idv.append(indiv_data[ind][0]) #to get name of individual
    elif carolina_per >= 0.95:
        carolina_total += 1
        carolina_idv.append(indiv_data[ind][0])
        carolina.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped' + '\t' + str(
            scalloped_allele) + '\t' + str(scalloped_per) + '\n' + 'carolina' + '\t' + str(
            carolina_allele) + '\t' + str(carolina_per) + '\n' + '-------------------------------' + '\n')
    else:
        unknown_idv.append(indiv_data[ind][0])
        unknown.write(str(indiv_data[ind][0]) + '\n' + str(total_alleles) + '\n' + 'scalloped' + '\t' + str(scalloped_allele) + '\t' + str(scalloped_per) + '\n' + 'carolina' + '\t' + str(carolina_allele) + '\t' + str(carolina_per) + '\n' + '-------------------------------' + '\n')
    print '-------------------------------'
    print indiv_data[ind][0]
    print total_alleles
    print 'scalloped' '\t', scalloped_allele, '\t', scalloped_per
    print 'carolina' '\t', carolina_allele, '\t', carolina_per

unknown.close()
scalloped.close()
carolina.close()
print scalloped_total
print carolina_total
print carolina_idv
print scalloped_idv
print unknown_idv

# export list of ID'd individuals into repsective species files
scalloped = open('Results/SLH3/SLH3.scalloped_ID.txt', 'w')
scall_write = open('Results/SLH3/SMH3_scall.bat', 'w')
for sp in scalloped_idv:
    scalloped.write(str(sp) + '\n')
    scall_write.write('ln -s ~/SLH3/' + str(sp) + '*.fq.gz .' + '\n')
scalloped.close()
scall_write.close()
carolina = open('Results/SLH3/SLH3.carolina_ID.txt', 'w')
carolina_write = open('Results/SLH3/SMH3_car.bat', 'w')
for sp in carolina_idv:
    carolina.write(str(sp) + '\n')
    carolina_write.write('ln -s ~/SLH3/' + str(sp) + '*.fq.gz .' + '\n')
carolina.close()
carolina.close()
unknown_namesonly = open('Results/SLH3/SLH3.unknown_scall_car_ID.txt', 'w')
unknown_write = open('Results/SLH3/SMH3_unknown.bat', 'w')
for sp in unknown_idv:
    unknown_namesonly.write(str(sp) + '\n')
    unknown_write.write('ln -s ~/SLH3/' + str(sp) + '*.fq.gz .' + '\n')
unknown_namesonly.close()
unknown_write.close()
