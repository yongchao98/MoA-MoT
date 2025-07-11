def analyze_hybrid_zone_genetics():
    """
    This script explains why a high mutation rate is not a consequence of gene flow,
    and demonstrates that other metrics, like Fis, can be high in a hybrid zone.
    """

    # --- Step 1: Define two distinct parent populations in HWE ---
    # Population 1 has a high frequency of allele 'A'
    p1 = 0.9  # Frequency of allele 'A' in Pop 1
    q1 = 1 - p1 # Frequency of allele 'a' in Pop 1
    N1 = 1000 # Number of individuals in Pop 1

    # Population 2 has a high frequency of allele 'a'
    p2 = 0.1  # Frequency of allele 'A' in Pop 2
    q2 = 1 - p2 # Frequency of allele 'a' in Pop 2
    N2 = 1000 # Number of individuals in Pop 2

    # --- Step 2: Create a hybrid zone by mixing the populations ---
    N_hybrid = N1 + N2

    # Calculate overall allele frequencies in the mixed population
    # Total 'A' alleles = (p1 * N1) + (p2 * N2)
    p_hybrid = (p1 * N1 + p2 * N2) / N_hybrid
    q_hybrid = 1 - p_hybrid

    # --- Step 3: Calculate Observed vs. Expected Heterozygosity (Wahlund Effect) ---

    # Observed Heterozygosity (Ho): The actual proportion of heterozygotes.
    # We assume the parent populations were in HWE, so Ho in Pop1 is 2*p1*q1.
    # The total number of heterozygotes is the sum from each source population.
    het_count_1 = (2 * p1 * q1) * N1
    het_count_2 = (2 * p2 * q2) * N2
    Ho_hybrid = (het_count_1 + het_count_2) / N_hybrid

    # Expected Heterozygosity (He): The proportion of heterozygotes expected
    # if the mixed gene pool were in Hardy-Weinberg Equilibrium.
    He_hybrid = 2 * p_hybrid * q_hybrid

    # --- Step 4: Calculate Fis and print the explanation ---
    # Fis measures the heterozygote deficit.
    Fis = (He_hybrid - Ho_hybrid) / He_hybrid

    print("--- Analyzing Population Genetics in a Hybrid Zone ---")
    print("When two distinct populations mix, as in a hybrid zone, it creates a 'Wahlund effect'.")
    print("This leads to a high Fis (inbreeding coefficient), demonstrating that C is a possible outcome.\n")

    print(f"Parent Population 1: Freq(A) = {p1}, Freq(a) = {q1}")
    print(f"Parent Population 2: Freq(A) = {p2}, Freq(a) = {q2}\n")

    print(f"Hybrid Zone Overall Freq(A) = {p_hybrid:.2f}")
    print(f"Hybrid Zone Overall Freq(a) = {q_hybrid:.2f}\n")

    print(f"Observed Heterozygosity (Ho) in Hybrid Zone: {Ho_hybrid:.4f}")
    print(f"Expected Heterozygosity (He) in Hybrid Zone: {He_hybrid:.4f}\n")

    print("The Fis statistic is calculated as (He - Ho) / He")
    print(f"Fis = ({He_hybrid:.4f} - {Ho_hybrid:.4f}) / {He_hybrid:.4f}")
    print(f"Fis = {Fis:.4f}\n")
    print("A positive Fis shows a deficit of heterozygotes, which is characteristic of the Wahlund effect.\n")

    print("--- Final Conclusion ---")
    print("Gene flow is the MOVEMENT of alleles between populations.")
    print("Mutation (rate µ) is the CREATION of new alleles.")
    print("Gene flow does not cause a change in the fundamental mutation rate.")
    print("Therefore, a high mutation rate cannot be an outcome of gene flow across a hybrid zone.")

if __name__ == "__main__":
    analyze_hybrid_zone_genetics()