import sys

def solve():
    """
    This script models a simple hybrid zone to analyze genetic metrics.
    We consider two parental populations that are completely diverged at a single gene,
    and a hybrid population formed by mixing them.
    """
    # Define two parental populations, Pop1 (all 'A' allele) and Pop2 (all 'a' allele)
    # p is the frequency of 'A', q is the frequency of 'a'
    p1, q1 = 1.0, 0.0
    p2, q2 = 0.0, 1.0

    # The hybrid zone is formed by mixing them equally.
    # We calculate the genetic metrics for this scenario.
    print("Analysis of a simplified hybrid zone:")
    print("="*50)

    # A. Fst: Population Differentiation
    # Fst = (Ht - Hs) / Ht
    # Hs = average heterozygosity in subpopulations
    Hs1 = 2 * p1 * q1
    Hs2 = 2 * p2 * q2
    Hs = (Hs1 + Hs2) / 2
    # Ht = expected heterozygosity in the total (merged) population
    p_total = (p1 + p2) / 2
    Ht = 2 * p_total * (1 - p_total)
    Fst = (Ht - Hs) / Ht if Ht > 0 else 0
    print("A. High Fst between populations: CAN OCCUR")
    print(f"   Calculation: Fst = (Ht - Hs) / Ht = ({Ht:.2f} - {Hs:.2f}) / {Ht:.2f} = {Fst:.2f}")
    print("   Explanation: A high Fst between parental populations is a prerequisite for a hybrid zone.")
    print("-"*50)

    # B. Dxy: Absolute Divergence
    print("B. High Dxy between populations: CAN OCCUR")
    print("   Explanation: High Dxy reflects the historical time since the two populations diverged. Gene flow doesn't erase this history.")
    print("-"*50)

    # C. Fis: Inbreeding Coefficient
    # Fis = (He - Ho) / He
    # He = expected heterozygosity in the hybrid zone based on its allele frequencies
    He_hybrid = 2 * p_total * (1 - p_total)
    # Ho = observed heterozygosity. In a simple mix, we only have parental genotypes ('AA' and 'aa'), so Ho = 0.
    Ho_hybrid = 0.0
    Fis = (He_hybrid - Ho_hybrid) / He_hybrid if He_hybrid > 0 else 0
    print("C. High Fis within a population: CAN OCCUR")
    print(f"   Calculation: Fis = (He - Ho) / He = ({He_hybrid:.2f} - {Ho_hybrid:.2f}) / {He_hybrid:.2f} = {Fis:.2f}")
    print("   Explanation: The mixture of distinct populations (Wahlund effect) causes a deficit of heterozygotes, leading to high Fis.")
    print("-"*50)
    
    # E. Pi: Nucleotide Diversity
    print("E. High Pi within a population: CAN OCCUR")
    print("   Explanation: Pi is a measure of genetic diversity. The hybrid zone contains alleles from both parental populations, so its diversity is high.")
    print("-"*50)

    # D. mu: Mutation Rate
    print("D. High mu (mutation rate) within a population: CANNOT OCCUR (as a result of gene flow)")
    print("   Explanation: The mutation rate (mu) is a basic biological property. Gene flow is the movement of alleles and does not change the rate at which new mutations arise. The two processes are independent.")
    print("="*50)


solve()
<<<D>>>