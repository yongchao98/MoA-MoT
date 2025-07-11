def analyze_hybrid_zone_genetics():
    """
    Analyzes the effects of gene flow across a hybrid zone on various population genetic metrics.
    A hybrid zone is a region where two distinct populations meet and interbreed.
    Gene flow is the transfer of alleles (gene variants) from one population to another.
    The script evaluates which of the given scenarios cannot occur as a result of this process.
    """
    print("Analyzing the consequences of gene flow in a hybrid zone...")
    print("-" * 60)

    # Option A: High Fst between populations
    print("A. High Fst between populations:")
    print("   - Fst measures genetic differentiation between populations. High Fst means populations are very different.")
    print("   - Gene flow acts to make populations more similar, which *reduces* Fst.")
    print("   - However, high Fst can still exist if gene flow is counteracted by strong selection against hybrids, or if the hybrid zone is very recent.")
    print("   - Verdict: High Fst is POSSIBLE.")
    print("-" * 60)

    # Option B: High Dxy between populations
    print("B. High Dxy between populations:")
    print("   - Dxy measures the absolute average number of DNA differences between two populations.")
    print("   - It reflects the long-term history of divergence. If two populations were separated for a long time before forming a hybrid zone, they will have a high Dxy.")
    print("   - Gene flow is a recent event and does not immediately erase this historical divergence.")
    print("   - Verdict: High Dxy is POSSIBLE and often expected.")
    print("-" * 60)

    # Option C: High Fis within a population
    print("C. High Fis within a population:")
    print("   - Fis measures inbreeding. A high positive Fis indicates a deficit of heterozygotes compared to random mating (Hardy-Weinberg Equilibrium).")
    print("   - A sample from a hybrid zone contains a mix of individuals from two distinct parental populations. This mixture, known as the Wahlund effect, results in a statistical deficit of heterozygotes.")
    print("   - Verdict: High Fis is POSSIBLE and is a classic feature of hybrid zones.")
    print("-" * 60)

    # Option D: High u (mutation rate) within a population
    print("D. High u (mutation rate) within a population:")
    print("   - 'u' (or µ) is the rate at which new mutations arise spontaneously in DNA during replication.")
    print("   - This is a fundamental molecular/biochemical property of a species. It is not a demographic parameter.")
    print("   - Gene flow is a demographic process involving the movement and mating of individuals; it mixes existing alleles but does not change the rate at which new ones are created.")
    print("   - Verdict: A high mutation rate is NOT caused by gene flow. This CANNOT occur as a consequence of the process.")
    print("-" * 60)

    # Option E: High Pi within a population
    print("E. High Pi within a population:")
    print("   - Pi (π) is nucleotide diversity, a measure of genetic variation *within* a population.")
    print("   - Gene flow from a divergent population introduces new alleles. This mixing of two different allele pools significantly increases the total genetic variation within the hybrid zone.")
    print("   - Verdict: High Pi is POSSIBLE and is an expected outcome.")
    print("-" * 60)

    print("Conclusion: Gene flow is a process of mixing existing genetic material. It affects measures of diversity (Pi), differentiation (Fst), divergence (Dxy), and population structure (Fis). It does not, however, alter the fundamental rate (u) at which new mutations are generated.")

if __name__ == '__main__':
    analyze_hybrid_zone_genetics()
    print("\n<<<D>>>")