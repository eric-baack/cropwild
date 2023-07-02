# maxdiffpy   Eric Baack, Luther College
# July 2022.  
# Model introgression and drift in crop-wild populations.
# Populations start with all wild alleles, but experience
# Two bouts of gene flow.
# Sample 20 individuals in the population to mirror experimental design. 
# What is the maximum crop allele frequency observed at the end?


import random
import copy
import sys


# set model constants
LOCI_PER_LG = 100  #SNPs to model per LG
NPOP = 100       # effective population size for crop-wild populations
LG_COUNT = 17   # number of chromosomes per individual
TRIALS = 20000  # number of trials to run for calculation of confidence interval
YEARS = 25 # number of years to simulate in model
GENE_FLOW_LIST = [0,4]  # list of years when gene flow occurs
F1_RATE = 0.05 # fraction of population that will be F1 in years when gene flow occurs
POPULATIONS = 7 #number of independent populations pooled for frequency analysis
SAMPLE_SIZE = 20 # number of individuals genotyped per population

def init_wildpop():
    population = []
    for i in range(NPOP):
        individual = init_chromos()
        population.append(individual)
    return population  # [individual][LG][homolog][locus]

def init_chromos():
    individual = []
    for i in range(LG_COUNT):
        gstring = "0" * LOCI_PER_LG  # individuals initially have no crop SNPs, only wild
        hstring = "0" * LOCI_PER_LG
        chromo = [gstring, hstring]
        individual.append(chromo)
    return individual

def recombine(individual): 
    gamete = []
    for i in range(LG_COUNT):
        gstring = individual[i][0]
        hstring = individual[i][1]
        rpoint = random.randrange(1,LOCI_PER_LG)
        g1 = gstring[0:rpoint] + hstring[rpoint:]
        g2 = hstring[0:rpoint] + gstring[rpoint:]
        x = random.random()   # randomly select one of four gametes to use in next generation
        if x < .25:
            gamete.append(g1)
        elif x < .5:
            gamete.append(g2)
        elif x < .75:
            gamete.append(gstring)
        else:
            gamete.append(hstring)
    return gamete

def reproduce(population):
    x = random.randrange(0,NPOP)
    y = random.randrange(0, NPOP)
    g1 = recombine(population[x])  # choose a randomly parent from the population, generate gametes, including recominant
    g2 = recombine(population[y])
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen  # returns individual 

def crop_pollen():
    gamete = []
    for i in range(LG_COUNT):
        gstring = "1" * LOCI_PER_LG  # crop pollen has only crop alleles
        gamete.append(gstring)
    return gamete

def geneflow(population):
    x = random.randrange(0,NPOP)
    g1 = recombine(population[x])  # pick individual from population to receive crop pollen adn generate gamete
    g2 = crop_pollen()
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen

def make_replicate():
    population = init_wildpop()
    for year in range(YEARS):
        mctr = 0
        if year in GENE_FLOW_LIST:
            mctr = int(NPOP * F1_RATE)  # obtain number of F1 progeny
        newpop = []
        for i in range(mctr):
            nextgen = geneflow(population)
            newpop.append(nextgen)
        for i in range(NPOP - mctr):
            nextgen = reproduce(population)
            newpop.append(nextgen)
        population = copy.deepcopy(newpop)
    return population


def calc_lg_freqs(lgfreqlist, simlist):
    for i in range(LG_COUNT):
        for j in range(LOCI_PER_LG):
            croptot = 0
            for population in simlist:
                for k in range(SAMPLE_SIZE):
                    cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                    croptot = croptot + cropct
            cropfreq = croptot / (SAMPLE_SIZE * POPULATIONS * 2)
            #print("LG, cropfreq: ", i, cropfreq)
            lgfreqlist.append(cropfreq)
    return lgfreqlist


def main():
    args = sys.argv
    rep = str(args[1])
    diffmaxlist = []
    diff20list = []
    outstr = "LGintrog_" + rep + ".txt"
    outfile = open(outstr, "w")
    for trial in range(TRIALS):
        if trial % 10 == 0:
                print("replicate: ", trial)
        lgfreqlist = []
        simlist = []
        for rep in range(POPULATIONS): 
            population = make_replicate()
            #output_genotypes(population)
            simlist.append(population)
        lgfreqlist = calc_lg_freqs(lgfreqlist, simlist)
        lgfreqlist.sort()
        diffmaxlist.append(lgfreqlist[1699])
        diff20list.append(lgfreqlist[1680])

    diffmaxlist.sort()
    diff20list.sort()
    print("Upper 95% CI for difference:", diffmaxlist[int(.975*TRIALS)])
    print("upper 95% CI for 20th rank difference", diff20list[int(.975*TRIALS)])
    print("Max difference observed; ", diffmaxlist[TRIALS-1])

main()