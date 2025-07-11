import math

def analyze_hybrid_zone():
    """
    Simulates gene flow between two distinct populations to analyze the effect
    on the inbreeding coefficient (Fis) within the hybrid zone.
    """
    # Step 1: Define two distinct parent populations.
    # Population 1 is fixed for allele 'A'. Frequency of 'A' (p1) is 1.0.
    # Population 2 is fixed for allele 'a'. Frequency of 'a' (q2) is 1.0.
    # This scenario represents high differentiation (high Fst) between the parent populations.
    print("--- Parent Populations ---")
    print("Population 1 is fixed for allele 'A'.")
    print("Population 2 is fixed for allele 'a'.")
    print("These populations are highly differentiated (high Fst and Dxy).\n")

    # Step 2: Simulate a hybrid population from equal mixing.
    # The allele frequencies in the new hybrid population are the average of the parents.
    p_hybrid = 0.5
    q_hybrid = 0.5
    print("--- Hybrid Zone Population ---")
    print(f"Gene flow occurs, creating a hybrid zone. The new allele frequencies are:")
    print(f"Frequency of 'A' (p) = {p_hybrid}")
    print(f"Frequency of 'a' (q) = {q_hybrid}\n")

    # Step 3: Calculate Observed and Expected Heterozygosity.
    # When an individual from Pop 1 (genotype AA) mates with an individual
    # from Pop 2 (genotype aa), all F1 offspring are heterozygotes (Aa).
    # Therefore, the observed heterozygosity (Ho) in the F1 hybrid zone is 1.0.
    Ho = 1.0
    print(f"In the first generation of hybrids, all individuals are 'Aa'.")
    print(f"Observed Heterozygosity (Ho) = {Ho}\n")

    # The expected heterozygosity (He) under Hardy-Weinberg Equilibrium is 2*p*q.
    He = 2 * p_hybrid * q_hybrid
    print(f"Expected Heterozygosity (He) based on allele frequencies (2*p*q) is:")
    print(f"He = 2 * {p_hybrid} * {q_hybrid} = {He}\n")

    # Step 4: Calculate the Inbreeding Coefficient (Fis).
    # Fis = 1 - (Ho / He). A high positive Fis (>0) means fewer heterozygotes than expected.
    # A negative Fis means more heterozygotes than expected.
    Fis = 1 - (Ho / He)
    print("--- Final Calculation: Fis ---")
    print("Fis measures the level of heterozygosity relative to expectations.")
    print("Fis = 1 - (Observed Heterozygosity / Expected Heterozygosity)")
    print(f"Fis = 1 - ({Ho} / {He})")
    print(f"Fis = {Fis}\n")

    # Conclusion
    print("Conclusion:")
    print(f"The resulting Fis is {Fis}, which is highly negative. This indicates a large")
    print("excess of heterozygotes, which is the opposite of what causes a high, positive Fis.")
    print("Therefore, high Fis within a population cannot occur as a direct result of")
    print("gene flow across a hybrid zone. In fact, gene flow causes the exact opposite effect.")

# Run the analysis
analyze_hybrid_zone()
