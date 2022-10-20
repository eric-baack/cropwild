# cropwild.py
# June 2022.  
# Model drift, recombination, gene flow to create null model for PCadapt.


import random
import copy

LOCI_PER_LG = 100
NPOP = 1000
LG_COUNT = 17

def init_wildpop():
    population = []
    for i in range(NPOP):
        individual = init_chromos()
        population.append(individual)
    return population

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

def main():
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
   
    # get ending allele freqs
    freqdict = {}
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

main()