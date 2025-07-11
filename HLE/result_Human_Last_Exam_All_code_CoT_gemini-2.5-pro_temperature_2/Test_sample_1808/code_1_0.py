import pandas as pd

def analyze_hybrid_zone():
    """
    Simulates gene flow into a hybrid zone and calculates population genetics statistics
    to determine the most unlikely outcome.
    """

    # --- Step 1: Define Parental Populations ---
    # Parental Population 1 (P1): Fixed for allele 'A' (Genotype 'AA')
    # Parental Population 2 (P2): Fixed for allele 'B' (Genotype 'BB')
    # Let's assume a sample size of 50 from each for the gene flow event.
    p1_migrants = 50
    p2_migrants = 50
    total_migrants = p1_migrants + p2_migrants
    
    print("--- Simulating Gene Flow ---")
    print(f"A hybrid zone is formed by {p1_migrants} individuals from Population 1 (all 'AA') and {p2_migrants} individuals from Population 2 (all 'BB').")

    # --- Step 2: Calculate Allele Frequencies in the new Hybrid Zone ---
    # Total alleles in the new population = 2 * total_migrants
    # Number of 'A' alleles = 2 * p1_migrants
    # Number of 'B' alleles = 2 * p2_migrants
    p_A = (2 * p1_migrants) / (2 * total_migrants)
    p_B = (2 * p2_migrants) / (2 * total_migrants)

    print(f"\nAllele frequency of 'A' in the hybrid zone (p_A) = {p_A:.2f}")
    print(f"Allele frequency of 'B' in the hybrid zone (p_B) = {p_B:.2f}")

    # --- Step 3: Simulate Random Mating and Calculate Genotype Frequencies ---
    # After one generation of random mating, genotype frequencies will be in Hardy-Weinberg Equilibrium.
    freq_AA = p_A**2
    freq_AB = 2 * p_A * p_B
    freq_BB = p_B**2
    
    print("\n--- Analyzing the Hybrid Zone Population After Random Mating ---")
    print("Assuming random mating (the essence of gene flow), we expect the following genotype frequencies:")
    print(f"Frequency of 'AA' = p_A^2 = {freq_AA:.2f}")
    print(f"Frequency of 'AB' (Heterozygotes) = 2*p_A*p_B = {freq_AB:.2f}")
    print(f"Frequency of 'BB' = p_B^2 = {freq_BB:.2f}")
    
    # --- Step 4: Calculate Pi and Fis for the Hybrid Zone ---
    # Expected Heterozygosity (He), which is equivalent to Nucleotide Diversity (Pi) here.
    pi_or_he = 2 * p_A * p_B
    
    # Observed Heterozygosity (Ho) is the frequency of heterozygotes after mating.
    ho = freq_AB
    
    # Fis = (He - Ho) / He
    # A special case here: because mating is random, He and Ho will be identical.
    # To avoid division by zero if He is 0, we'll handle it.
    if pi_or_he > 0:
        fis = (pi_or_he - ho) / pi_or_he
    else:
        fis = 0
        
    print("\n--- Calculating Population Statistics from Answer Choices ---")
    
    # E. High Pi within a population
    print(f"\nE. Pi (Nucleotide Diversity within the population) is equal to the expected heterozygosity.")
    print(f"   Pi = 2 * p_A * p_B = {pi_or_he:.2f}")
    print("   This value is maximized because gene flow brought two different alleles together. Thus, 'High Pi' occurs.")

    # C. High Fis within a population
    print(f"\nC. Fis (Inbreeding Coefficient) measures the heterozygote deficit.")
    print(f"   Fis = (Expected_Het - Observed_Het) / Expected_Het")
    print(f"   Fis = ({pi_or_he:.2f} - {ho:.2f}) / {pi_or_he:.2f} = {fis:.2f}")
    print("   Fis is 0, indicating no deficit of heterozygotes. In fact, gene flow creates heterozygotes, which opposes a high Fis. Therefore, 'High Fis' cannot occur as a result of gene flow.")
    
    print("\n--- Analysis of Other Options ---")
    # A. High Fst
    print("A. High Fst: Fst measures differentiation between populations. While gene flow reduces Fst, a hybrid zone often exists in a state of tension where gene flow is balanced by selection, allowing high Fst between the parental populations to be maintained. This *can* occur.")
    # B. High Dxy
    print("B. High Dxy: Dxy measures the absolute genetic divergence, reflecting the time since two populations split. Gene flow does not erase this history. If the parent populations were separated for a long time, Dxy will be high. This *can* occur.")
    # D. High u
    print("D. High u (mutation rate): Mutation rate is a basic biological property of an organism's DNA replication machinery. It is independent of mating patterns like gene flow. A population could have a high mutation rate for other reasons. This *can* occur.")

    print("\n--- Conclusion ---")
    print("The process of gene flow, which involves mating between populations, actively creates heterozygotes. This drives the Fis value towards zero (or even negative if there's an excess of heterozygotes). A high Fis signifies a strong deficit of heterozygotes, which is the opposite of what gene flow accomplishes.")

analyze_hybrid_zone()