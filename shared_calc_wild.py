#shared_calc.py
# June 2022
# Given 23391 snps, outliers, 3 regions, how many expect shared across all, across 2?
OUTLIERS = [271,419, 293]
SITES = 23391
SIMS = 20000
import random

def get_share():
    regionA = set()
    regionB = set()
    regionC = set()
    for i in range(OUTLIERS[0]):
        x = random.randrange(0, SITES)
        regionA.add(x)
    for i in range(OUTLIERS[1]):
        x = random.randrange(0, SITES)
        regionB.add(x)
    for i in range(OUTLIERS[2]):
        x = random.randrange(0,SITES)
        regionC.add(x)
    ABshare = regionA.intersection(regionB)
    BCshare = regionB.intersection(regionC)
    ACshare = regionA.intersection(regionC)
    ABCshare = ABshare.intersection(regionC)
    share2 =len(ABshare) + len(BCshare) + len(ACshare)
    return share2, len(ABCshare)

def main():
    s2list = []
    s3list = []
    for i in range(SIMS):
        share2, share3 = get_share()
        s2list.append(share2)
        s3list.append(share3)
    s2list.sort()
    s3list.sort()
    print("lower, upper95% CI share 3: ", s3list[int(SIMS*.025)], s3list[int(SIMS*.975)])
    print("lower, upper 95% CI share 2 upper", s2list[int(SIMS*.025)], s2list[int(SIMS*.975)])
   

main()
