def analyze_sex_differentiation():
    """
    Simulates and explains why Fst can be high between sexes for specific markers.
    """
    # Define parameters for a simulated population
    num_males = 100
    num_females = 100

    # --- Marker 1: An Autosomal Marker ---
    # Alleles are 'A' and 'a'. We assume they are distributed randomly
    # in the overall population gene pool that both sexes draw from.
    # Let's say the frequency of 'A' is 0.7 in the population.
    freq_A_autosomal = 0.7
    male_autosomal_freq = freq_A_autosomal
    female_autosomal_freq = freq_A_autosomal

    # --- Marker 2: A Sex-Linked Marker (on the Y chromosome in an XY system) ---
    # The marker is present on the Y chromosome.
    # Therefore, all males have it, and no females do.
    # Let's define the frequency of the 'Y_marker' allele.
    male_y_linked_freq = 1.0
    female_y_linked_freq = 0.0

    # --- Output the analysis ---
    print("=" * 65)
    print("Analysis of Genetic Differentiation Between Males and Females")
    print("=" * 65)
    print(f"Population: {num_males} Males, {num_females} Females")
    print("\n--- CASE 1: Autosomal Marker ---")
    print("This marker is on a non-sex chromosome.")
    print("Equation for Allele Frequency in Males: P(A)_males = Population Frequency")
    print(f"Resulting Allele 'A' Frequency in Males: {male_autosomal_freq:.1f}")
    print("\nEquation for Allele Frequency in Females: P(A)_females = Population Frequency")
    print(f"Resulting Allele 'A' Frequency in Females: {female_autosomal_freq:.1f}")
    print("\n>>> Conclusion: Frequencies are identical. This results in Fst ≈ 0.")

    print("\n" + "-" * 65 + "\n")

    print("--- CASE 2: Y-Linked Sex Marker (in an XY system) ---")
    print("This marker is only on the Y chromosome.")
    print("Equation for Allele Frequency in Males: P(marker)_males = 1 (by definition)")
    print(f"Resulting Marker Frequency in Males: {male_y_linked_freq:.1f}")
    print("\nEquation for Allele Frequency in Females: P(marker)_females = 0 (no Y chromosome)")
    print(f"Resulting Marker Frequency in Females: {female_y_linked_freq:.1f}")
    print("\n>>> Conclusion: Frequencies are maximally different. This results in Fst = 1.")

    print("\n" + "=" * 65)
    print("\nFinal Explanation:")
    print("The fact that *some* markers show high Fst while others do not points directly to the markers' location on sex chromosomes (X/Y or Z/W), which have different inheritance patterns and composition between sexes.")

# Run the analysis
analyze_sex_differentiation()
