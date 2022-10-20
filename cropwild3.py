# cropwild2.py
# July 2022.  
# Model drift, recombination, gene flow, selection to create model for figure 4.
# 10% of loci are unfavorable for crop.  Effects are additive, 
# What is the distribution of track lengths, tracks per chromosome
# Sample 20 individuals in the population to mirror our data.  
# Potentially vary introgression frequency.
# Introgression rate = % of snps that are crop
# Expectation:  little variation between LGs in introgression rate in this model given sampling.
# Interpretation:  drift alone unlikely to explain this.  Some regions are less permissive (= stronger selection)



import random
import copy
import statistics
import math

LOCI_PER_LG = 100
NPOP = 100
LG_COUNT = 17

def init_wildpop():
    population = []
    for i in range(NPOP):
        individual = init_chromos()
        population.append(individual)
    return population  # [individual][LG][homolog][locus]

def init_chromos():
    individual = []
    for i in range(LG_COUNT):
        gstring = "0" * LOCI_PER_LG
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

def reproduce(population):
    x = random.randrange(0,NPOP)
    y = random.randrange(0, NPOP)
    g1 = recombine(population[x])
    g2 = recombine(population[y])
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen

def crop_pollen():
    gamete = []
    for i in range(LG_COUNT):
        gstring = "1" * LOCI_PER_LG
        gamete.append(gstring)
    return gamete

def geneflow(population):
    x = random.randrange(0,NPOP)
    g1 = recombine(population[x])
    g2 = crop_pollen()
    nextgen = []
    for i in range(LG_COUNT):
        chromo = [g1[i], g2[i]]
        nextgen.append(chromo)
    return nextgen

def make_replicate():
    population = init_wildpop()
    for year in range(25):
        mctr = 0
        if year == 0 or year == 4 or year == 7:
            mctr = 5
        newpop = []
        for i in range(mctr):
            nextgen = geneflow(population)
            newpop.append(nextgen)
        for i in range(NPOP - mctr):
            nextgen = reproduce(population)
            newpop.append(nextgen)
        population = copy.deepcopy(newpop)
    return population

def output_genotypes(population):
    # get ending allele freqs
    freqlist = []
    outstring = ""
    for i in range(LG_COUNT):
        for j in range(LOCI_PER_LG):
            indgenos = []
            for k in range(NPOP):
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
    outfile = open("/home/baacer01/popgen/pcatest100.txt", "w")
    outfile.write(outstring)
    outfile.close()
    freqlist.sort()
    bound95 = int(LG_COUNT*LOCI_PER_LG*.95)
    print("95pct upper CI for crop freq", freqlist[bound95])

def calc_lg_freqs(lgfreqdict, population):
    lgfreqdict = {}
    for i in range(LG_COUNT):
        croptot = 0
        for j in range(LOCI_PER_LG):
            for k in range(20):
                cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                croptot = croptot + cropct
        cropfreq = croptot / (20 * 2 * LOCI_PER_LG)
        print("LG, cropfreq: ", i, cropfreq)
        lgfreqdict[i] = cropfreq
    return lgfreqdict

def main():
    lgfreqdict = {}
    outfile = open("/home/baacer01/popgen/LGintrog.txt", "w")
    for rep in range(22):
        population = make_replicate()
        #output_genotypes(population)
        lgfreqdict = calc_lg_freqs(lgfreqdict, population)
        for i in range(LG_COUNT):
            outline = str(rep) + "," + str(i) + "," + str(lgfreqdict[i]) + "\n"
            outfile.write(outline)
    outfile.close()

main()