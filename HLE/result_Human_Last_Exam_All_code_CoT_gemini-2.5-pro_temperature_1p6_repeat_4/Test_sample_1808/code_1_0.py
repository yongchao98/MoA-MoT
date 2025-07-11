def analyze_hybrid_zone():
    """
    Analyzes population genetic metrics in a simplified hybrid zone scenario
    to determine which option cannot be caused by gene flow.
    """
    # --- Model Setup ---
    # We consider a single genetic locus with two alleles, '0' and '1'.
    # Pop 1 is initially fixed for allele '0'. Its allele frequency for '1' is 0.0.
    # Pop 2 is initially fixed for allele '1'. Its allele frequency for '1' is 1.0.
    p1 = 0.0  # Frequency of allele '1' in Pop 1
    p2 = 1.0  # Frequency of allele '1' in Pop 2

    print("--- Analyzing Population Metrics in a Hybrid Zone ---")
    print(f"Parental Population 1 has allele '1' frequency (p1) = {p1}")
    print(f"Parental Population 2 has allele '1' frequency (p2) = {p2}\n")

    # A. High Fst between populations
    # Fst = (Ht - Hs) / Ht, where Ht is total heterozygosity and Hs is subpopulation heterozygosity.
    p_total = (p1 + p2) / 2.0  # Average allele frequency
    Ht = 2 * p_total * (1 - p_total) # Expected heterozygosity in total pop = 2 * 0.5 * 0.5 = 0.5
    H1 = 2 * p1 * (1 - p1) # Expected heterozygosity in Pop 1 = 2 * 0.0 * 1.0 = 0.0
    H2 = 2 * p2 * (1 - p2) # Expected heterozygosity in Pop 2 = 2 * 1.0 * 0.0 = 0.0
    Hs = (H1 + H2) / 2.0
    Fst = (Ht - Hs) / Ht if Ht > 0 else 0
    print(f"A. Fst: The initial Fst between the parental populations is {Fst:.2f} (Ht={Ht:.2f}, Hs={Hs:.2f}).")
    print("   This is the maximum possible differentiation. A hybrid zone exists because of this differentiation,")
    print("   so high Fst can certainly occur. VERDICT: POSSIBLE\n")

    # B. High Dxy between populations
    # Dxy = p1*(1-p2) + (1-p1)*p2
    Dxy = p1 * (1 - p2) + (1 - p1) * p2
    print(f"B. Dxy: The absolute divergence (Dxy) is {Dxy:.2f} = ({p1} * (1-{p2})) + ((1-{p1}) * {p2}).")
    print("   This reflects the complete initial divergence between populations. VERDICT: POSSIBLE\n")

    # C. High Fis within a population
    # Imagine a hybrid population formed by an equal mix of individuals from Pop1 (genotype 0/0)
    # and Pop2 (genotype 1/1).
    # The allele frequency of '1' in this mixed population (p_mix) is 0.5.
    p_mix = 0.5
    # The expected heterozygosity (He) is 2 * p_mix * (1 - p_mix).
    He = 2 * p_mix * (1 - p_mix)
    # The observed heterozygosity (Ho) is 0, because we only have 0/0 and 1/1 individuals.
    Ho = 0.0
    # Fis = (He - Ho) / He
    Fis = (He - Ho) / He if He > 0 else 0
    print(f"C. Fis: In a zone of admixture, we can see the Wahlund effect.")
    print(f"   Expected heterozygosity (He) is {He:.2f}. Observed heterozygosity (Ho) is {Ho:.2f}.")
    print(f"   This results in a high Fis = ({He:.2f} - {Ho:.2f}) / {He:.2f} = {Fis:.2f}. VERDICT: POSSIBLE\n")

    # E. High Pi within a population
    # Pi (Nucleotide Diversity) is the same as He for a single locus.
    # The parental populations have Pi = 0.
    # The hybrid zone population has a Pi that is higher.
    Pi_parental = H1 # or H2
    Pi_hybrid = He
    print(f"E. Pi: Gene flow increases variation. Parental Pi was {Pi_parental:.2f}.")
    print(f"   The Pi in the hybrid zone is {Pi_hybrid:.2f}, which is an increase. VERDICT: POSSIBLE\n")
    
    # D. High u (mutation rate)
    print("D. Mutation Rate (u): The mutation rate is the rate of new allele formation.")
    print("   Gene flow is the *movement* of existing alleles between populations.")
    print("   One process does not cause the other. Gene flow cannot cause a high mutation rate.")
    print("   VERDICT: CANNOT OCCUR\n")

if __name__ == '__main__':
    analyze_hybrid_zone()