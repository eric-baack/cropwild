# cropwild4.py
# August 2022.  
# Model drift, recombination, gene flow to create null model for loci with no introgression (crop < .1)
# Sample 20 individuals in the population to mirror our data.  
# 


import random
import copy
import statistics
import math

LOCI_PER_LG = 100
NPOP = 100
LG_COUNT = 17
TRIALS = 2000
F1FRAC = .05

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
            mctr = int(F1FRAC*NPOP)
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
    cresistctlist = []
    outstring = ""
    cresistct = 0
    for i in range(LG_COUNT):
        for j in range(LOCI_PER_LG):
            indgenos = []
            for k in range(20):
                cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                #indgenos.append(cropct)
            cropfreq = sum(indgenos) / (NPOP * 2)
            #print("LG, Locus, freq", i, j, cropfreq)
            first = True
    #         for z in range(len(indgenos)):
    #             if first:
    #                 outstring = outstring + str(indgenos[z])
    #                 first = False
    #             else:
    #                 outstring = outstring +" " + str(indgenos[z])
    #         outstring = outstring + "\n"
    # outfile = open("/home/baacer01/popgen/pcatest100.txt", "w")
    # outfile.write(outstring)
    # outfile.close()
    
    # bound95 = int(LG_COUNT*LOCI_PER_LG*.95)
    # print("95pct upper CI for crop freq", cresistctlist[bound95])
    # print("crop resistant loci count, percent: ", cresistct, cresistct / (LG_COUNT*LOCI_PER_LG))

def calc_lg_freqs(lgfreqdict, population):
    for i in range(LG_COUNT):
        croptot = 0
        for j in range(LOCI_PER_LG):
            for k in range(20):
                cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                croptot = croptot + cropct
        cropfreq = croptot / (20 * 2 * LOCI_PER_LG)
        lgfreqdict[i].append(cropfreq)
    return lgfreqdict

# def anova(lgfreqdict):
#     lgmeanlist = []
#     lgvarlist = []
#     for i in range(LG_COUNT):
#         lmean = statistics.mean(lgfreqdict[i])
#         lvar = statistics.variance(lgfreqdict[i])
#         lgmeanlist.append(lmean)
#         lgvarlist.append(lvar)
#     tmean = statistics.mean(lgmeanlist)
#     ssg = 0
#     sse = 0
#     for i in range(LG_COUNT):
#         ssg = ssg + (tmean - lgmeanlist[i])**2
#         sse = sse + lgvarlist[i]*(21)
#     msg = ssg / 16
#     mse = sse / 357
#     fstat = msg / mse
#     return fstat

# def differ(lgfreqdict):
#     lgmeanlist = []
#     for i in range(LG_COUNT):
#         lmean = statistics.mean(lgfreqdict[i])
#         lgmeanlist.append(lmean)
#     maxmean = max(lgmeanlist)
#     minmean = min(lgmeanlist)
#     maxdiff = maxmean - minmean
#     return maxdiff

def init_accumulator(accumulator):
    for acc in accumulator:
        for i in range(LG_COUNT):
            indgenos = []
            for j in range(LOCI_PER_LG):
                indgenos.append(0)
            acc.append(indgenos)
    return accumulator
                

def crop_ctr(population, accumulator):
    for i in range(LG_COUNT):
        for j in range(LOCI_PER_LG):
            for k in range(20):
                cropct = int(population[k][i][0][j]) + int(population[k][i][1][j])
                accumulator[i][j] += cropct
    return accumulator

def find_cropresist(accumulator, cresistdict, reg):
    pcount = 7
    crct = 0
    if reg == 0:
        pcount = 8
    for lg in range(LG_COUNT):
        for locus in range(LOCI_PER_LG):
            cfreq = accumulator[reg][lg][locus] / (20*pcount*2)
            # if lg == 0:
            #     print("lg1, loci, cfreq:", locus, cfreq)
            if cfreq < .1:
                crct += 1
                entry = (lg, locus)
                if entry in cresistdict:
                    cresistdict[entry] += 1
                else:
                    cresistdict[entry] = 1
    # print("crct for region: ", reg, crct)
    return cresistdict # returns dictionary of lg, loci with cropfreq < .1 in at least one region, with count of regions

def main():
    share2pct = []
    share3pct = []
    for trial in range(TRIALS):
        if trial % 10 == 0:
                print("replicate: ", trial)
        # set accumulators.  For each region, LG_COUNT, LOCI_PER_LG.  regaccum[i][j] = count of crop alleles at lg i, locus j
        nregaccum = []
        mregaccum = []
        sregaccum = []
        accumulators = [nregaccum, mregaccum, sregaccum]
        accumulators = init_accumulator(accumulators)
        for rep in range(22):
            population = make_replicate()
            if rep < 8:
                nregaccum = crop_ctr(population, nregaccum)
            elif rep < 15:
                mregaccum = crop_ctr(population, mregaccum)
            else:
                sregaccum = crop_ctr(population, sregaccum)
        cresistdict = {}  # keys = lg, loci eg (0,37).  values = regions.
        for reg in range(3):
            cresistdict = find_cropresist(accumulators, cresistdict, reg)
        crloci = cresistdict.items()
        rct3 = 0
        rct2 = 0
        for i in crloci:
            if i[1] == 3:
                rct3+=1
            if i[1] == 2:
                rct2 += 1
        share2pct.append(rct2/len(crloci))  # 
        share3pct.append(rct3/len(crloci))
        if trial == 0:
            print("total number of crop-resistant loci:", len(crloci))
        # print("overall total crresist count (out of 1700):", len(crloci))
    share2pct.sort()
    share3pct.sort()
    lowci2 = share2pct[int(.025*TRIALS)]
    highci2 = share2pct[int(.975*TRIALS)]
    lowci3 = share3pct[int(.025*TRIALS)]
    highci3 = share3pct[int(.975*TRIALS)]
    print("Trials, F1 fraction:", TRIALS, F1FRAC)
    print("share2 95% CI:", lowci2, highci2)
    print("share3 95% CI:", lowci3, highci3)
        


main()