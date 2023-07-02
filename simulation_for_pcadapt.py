# cropwild.py
# June 2022.  
# Model drift, recombination, gene flow to create null model for PCadapt.
# Carry out 100 replicates, generating 100 data files to analyze with PCadapt.
# Datafiles are names 'pcatest100x.txt' where x ranges from 0 to 99.


import random
import copy

# Constants for model
LOCI_PER_LG = 100  # number of SNP loci per chromosome
NPOP = 100  # Sets effective size of population
YEARS = 25 # Sets number of years for model.
LG_COUNT = 17  # number of chromosomes
F1_RATE = 0.05  # fraction of F1s generated when gene flow occurs
GENEFLOW_YEAR_LIST = [0,5] # list of years when gene flow occurs
TRIALS = 100 # number of sets of simulations to run to obtain average outcome


def init_wildpop():  # create a population with Npop individuals, with chromosomes with Loci.
    population = []
    for i in range(NPOP):
        individual = init_chromos()
        population.append(individual)
    return population

def init_chromos():  # create string of SNPs for each chromosome
    individual = []
    for i in range(LG_COUNT):
        gstring = "0" * LOCI_PER_LG  #initial population has only wild alleles
        hstring = "0" * LOCI_PER_LG
        chromo = [gstring, hstring]
        individual.append(chromo)
    return individual

def recombine(individual): # create recominbinant gamete, one cross-over per chromosome
    gamete = []
    for i in range(LG_COUNT):
        gstring = individual[i][0]
        hstring = individual[i][1]
        rpoint = random.randrange(1,LOCI_PER_LG)
        g1 = gstring[0:rpoint] + hstring[rpoint:]
        g2 = hstring[0:rpoint] + gstring[rpoint:]
        x = random.random()
        if x < .25:
            gamete.append(g1)
        elif x < .5:
            gamete.append(g2)
        elif x < .75:
            gamete.append(gstring)
        else:
            gamete.append(hstring)
    return gamete

def reproduce(population): # for each individual in population, select parents at random, generate gametes (possibly with recombination)
    x = random.randrange(0,NPOP)
    y = random.randrange(0, NPOP)
    g1 = recombine(population[x])
    g2 = recombine(population[y])
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen  # returns a genotype for one individual for next generation

def crop_pollen():
    gamete = []
    for i in range(LG_COUNT):
        gstring = "1" * LOCI_PER_LG  #crop pollen has only crop alleles
        gamete.append(gstring)
    return gamete

def geneflow(population):   # create an F1 from crop pollen, wild ovule
    x = random.randrange(0,NPOP)  # choose a random wild parent
    g1 = recombine(population[x])
    g2 = crop_pollen()
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen  # returns a genotype for an F1 between crop and wild

def main():
    for repct in range(TRIALS):
        population = init_wildpop()  # create population of wild genotypes
        for year in range(YEARS):
            mctr = 0
            if year in GENEFLOW_YEAR_LIST:  # if year is one with gene flow, form F1s
                mctr = int(F1_RATE*NPOP)  # Calculate number of F1s formed from F1 rate, population size
            newpop = []
            for i in range(mctr):  # Create F1s for next year, up to the required number
                nextgen = geneflow(population)
                newpop.append(nextgen)
            for i in range(NPOP - mctr):  # Create non-F1s for next year
                nextgen = reproduce(population)
                newpop.append(nextgen)
            population = copy.deepcopy(newpop)
    
        # get ending allele freqs
        freqdict = {}
        freqlist = []
        outstring = ""
        for i in range(LG_COUNT):  # for each chromosome
            for j in range(LOCI_PER_LG):  # for each locus
                indgenos = []
                for k in range(NPOP):  # for each invidual
                    cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                    indgenos.append(cropct)
                cropfreq = sum(indgenos) / (NPOP * 2)
                #print("LG, Locus, freq", i, j, cropfreq)
                freqlist.append(cropfreq)
                first = True
                for z in range(len(indgenos)):
                    if first:
                        outstring = outstring + str(indgenos[z])
                        first = False
                    else:
                        outstring = outstring +" " + str(indgenos[z])
                outstring = outstring + "\n"
        outfname = "pcatest100" + str(repct) + ".txt"
        outfile = open(outfname, "w")
        outfile.write(outstring)
        outfile.close()
        freqlist.sort()
        bound95 = int(LG_COUNT*LOCI_PER_LG*.95)
        print("replicate, 95pct upper CI for crop freq", repct, freqlist[bound95])

main()