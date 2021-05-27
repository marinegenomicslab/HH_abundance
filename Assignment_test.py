'''
Script to determine how many diagnostic loci an individual must be genotyped at for accurate identification
Created: 12/14/2018
Last updated: 1/16/2019
AM Barker
'''

import csv
import pickle
import random


##### GREAT HAMMERHEAD VS SCALLOPED/CAROLINA TEST #####
## Diagnostic SNPs: 2695


# Import great vs scall/car dictionary of diagnostic SNPs
with open('Great_ScallCar_dict.pickle', 'rb') as data:
    Great_ScallCar_dict = pickle.load(data)

# Import test individuals allele data: text file imported as a list, each individual is its own list
with open('ScallCar_Great_inds.txt', 'rU') as f:
    indiv_data = list(csv.reader(f, delimiter = '\t'))

# Get contig names form first line, the delete from list
contig_names = indiv_data[0]
indiv_data.remove(indiv_data[0])

num_individuals = len(indiv_data)

# list of number of loci to test
num_loci = [5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]
#num_loci = [2, 5, 10]

num_iterations = 1000

iteration_results = []

for ni in range(0, num_iterations):
    # list to store results
    sim_results = []

    for sim in num_loci:
        # Select random loci
        rand_loci = random.sample(contig_names, sim)

        #find positions in contig list that match randomly selected loci
        contig_pos = []
        for loc in rand_loci:
            for position, item in enumerate(contig_names):
                if item == loc:
                    contig_pos.append(position)

        #print contig_pos
        #Add 1 to position in contig_pos. The first entry in indiv_data is sample name so need to add one to account for that when filtering data
        filtered_data = []
        for n in range(0,num_individuals):
            x = [indiv_data[n][i+1] for i in contig_pos]
            filtered_data.append(x)

        print len(filtered_data) #check right number of individuals
        print len(filtered_data[0]) #check right number of loci

        species_results = []
        great_idv = []
        scallCar_idv = []
        ind_check = []
        unknown_idv = []

        scalCar_total = 0
        great_total = 0

        # identify individuals
        for ind in range(0, len(filtered_data)):
            ind_id = []
            scallCar_allele = 0
            great_allele = 0
            for c in range(0, len(rand_loci)):
                contig = rand_loci[c]
                #pos = c + 1 # add 1 because first position is individual name
                pos = c
                alleles = str(filtered_data[ind][pos]) # alleles for each contig are together (130130) and need separated (130, 130)
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
                species_results.append("scalloped/carolina")
            elif great_per >= 0.95:
                great_total += 1
                great_idv.append(indiv_data[ind][0])
                species_results.append("Great")
            else:
                unknown_idv.append(indiv_data[ind][0])
                species_results.append("UND")
            print '-------------------------------'
            print indiv_data[ind][0]
            print total_alleles
            print 'scalloped/carolina' '\t', scallCar_allele, '\t', scallCar_per
            print 'great' '\t', great_allele, '\t', great_per

        print scalCar_total
        print great_total
        print scallCar_idv
        print great_idv
        print ni
        sim_results.append(species_results)

    iteration_results.append(sim_results)

# write results file
res_file = open('ScallCar_Great_sim_results.txt', 'w')

for ni in range(0, num_iterations):
    for sim in range(0, len(num_loci)):
        for ind in range(0, num_individuals):
            res_file.write(indiv_data[ind][0] + '\t' + str(ni) + '\t' + str(num_loci[sim]) + '\t' + iteration_results[ni][sim][ind] + '\n')

res_file.close()



##### SCALLOPED VS CAROLINA TEST #####
## Diagnostic SNPs; 1491


# Import scalloped vs carolina dictionary of diagnostic SNPs
with open('Scall_Car_dict.pickle', 'rb') as data:
    Scall_Car_dict = pickle.load(data)

# Import test individuals allele data: text file imported as a list, each individual is its own list
with open('Scall_Car_inds.txt', 'rU') as f:
    indiv_data = list(csv.reader(f, delimiter = '\t'))

# Get contig names form first line, the delete from list
contig_names = indiv_data[0]
indiv_data.remove(indiv_data[0])

num_individuals = len(indiv_data)

# list of number of loci to test
num_loci = [5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200] # 1 locus doesn't work, get division by zero error
#num_loci = [2, 5, 10]

num_iterations = 1000

iteration_results = []

for ni in range(0, num_iterations):
    # list to store results
    sim_results = []

    for sim in num_loci:
        # Select random loci
        rand_loci = random.sample(contig_names, sim)

        #find positions in contig list that match randomly selected loci
        contig_pos = []
        for loc in rand_loci:
            for position, item in enumerate(contig_names):
                if item == loc:
                    contig_pos.append(position)

        #print contig_pos
        #Add 1 to position in contig_pos. The first entry in indiv_data is sample name so need to add one to account for that when filtering data
        filtered_data = []
        for n in range(0,num_individuals):
            x = [indiv_data[n][i+1] for i in contig_pos]
            filtered_data.append(x)

        print len(filtered_data) #check right number of individuals
        print len(filtered_data[0]) #check right number of loci

        species_results = []
        scall_idv = []
        car_idv = []
        ind_check = []
        unknown_idv = []

        scall_total = 0
        car_total = 0

        # identify individuals
        for ind in range(0, len(filtered_data)):
            ind_id = []
            scall_allele = 0
            car_allele = 0
            for c in range(0, len(rand_loci)):
                contig = rand_loci[c]
                #pos = c + 1 # add 1 because first position is individual name
                pos = c
                alleles = str(filtered_data[ind][pos]) # alleles for each contig are together (130130) and need separated (130, 130)
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
                    scall_allele += 1
                if i == 'carolina':
                    car_allele += 1

            total_alleles = len(ind_id)
            scall_per = (float(scall_allele) / total_alleles)
            car_per = (float(car_allele) / total_alleles)
            if scall_per >= 0.95:
                scall_total += 1
                scall_idv.append(indiv_data[ind][0]) #to get name of individual
                species_results.append("Scalloped")
            elif car_per >= 0.95:
                car_total += 1
                car_idv.append(indiv_data[ind][0])
                species_results.append("Carolina")
            else:
                unknown_idv.append(indiv_data[ind][0])
                species_results.append("UND")
            print '-------------------------------'
            print indiv_data[ind][0]
            print total_alleles
            print 'scalloped' '\t', scall_allele, '\t', scall_per
            print 'carolina' '\t', car_allele, '\t', car_per

        print scall_total
        print car_total
        print car_idv
        print scall_idv
        print ni
        sim_results.append(species_results)

    iteration_results.append(sim_results)

# write results file
res_file = open('Scall_Car_sim_results.txt', 'w')

for ni in range(0, num_iterations):
    for sim in range(0, len(num_loci)):
        for ind in range(0, num_individuals):
            res_file.write(indiv_data[ind][0] + '\t' + str(ni) + '\t' + str(num_loci[sim]) + '\t' + iteration_results[ni][sim][ind] + '\n')

res_file.close()
