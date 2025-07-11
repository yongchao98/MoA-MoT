def calculate_fst_for_sex_chromosome():
    """
    This function explains and calculates the Fst value for a genetic marker
    on a sex-specific chromosome (e.g., the Y chromosome in an XY system).

    Fst is a measure of population differentiation, calculated as:
    Fst = (Ht - Hs) / Ht
    - Ht: The expected heterozygosity in the total, pooled population.
    - Hs: The average expected heterozygosity across the subpopulations.
    """

    print("Analyzing the cause of high Fst between males and females for specific markers.")
    print("The most direct explanation is a marker located on a sex-specific chromosome (e.g., Y or W).\n")

    print("--- Scenario: A marker on the Y-chromosome in an XY-system ---")
    print("Our two 'subpopulations' are Males and Females from one species.\n")

    # --- Step 1: Calculate Hs (Average Heterozygosity in Subpopulations) ---
    # In this context, heterozygosity (H) is 1 - sum(p_i^2), where p_i is the frequency of each allele.

    # For Males (XY): The Y-marker is present. Its frequency (p) is 1.0.
    p_males = 1.0
    H_males = 1 - (p_males**2)

    # For Females (XX): The Y-marker is absent. Its frequency (p) is 0.
    # The 'absence' allele has a frequency of 1.0.
    p_females = 0.0
    H_females = 1 - (1.0**2) # Based on the frequency of the 'absence' allele.

    # Hs is the average of the two.
    Hs = (H_males + H_females) / 2.0

    print("1. Calculate Hs (Average Subpopulation Heterozygosity):")
    print(f"   - Allele frequency in Males:   {p_males}")
    print(f"   - Heterozygosity in Males (H_males) = 1 - ({p_males}^2) = {H_males}")
    print(f"   - Allele frequency in Females: {p_females}")
    print(f"   - Heterozygosity in Females (H_females) = 1 - (1.0^2) = {H_females}")
    print(f"   - Average Heterozygosity (Hs) = ({H_males} + {H_females}) / 2 = {Hs}\n")


    # --- Step 2: Calculate Ht (Heterozygosity in Total Population) ---
    # Assuming a 1:1 sex ratio, half the population has the Y-marker allele (males)
    # and half does not (females).
    p_total = 0.5 # Frequency of the Y-marker allele in the total population.
    q_total = 0.5 # Frequency of the "absence" allele in the total population.
    Ht = 1 - (p_total**2 + q_total**2)

    print("2. Calculate Ht (Total Population Heterozygosity):")
    print(f"   - Overall frequency of Y-marker allele (p_total) = {p_total}")
    print(f"   - Overall frequency of 'absence' allele (q_total) = {q_total}")
    print(f"   - Total Heterozygosity (Ht) = 1 - ({p_total}^2 + {q_total}^2) = {Ht}\n")

    # --- Step 3: Calculate Final Fst ---
    # Avoid division by zero, although Ht > 0 here.
    Fst = (Ht - Hs) / Ht if Ht > 0 else 0

    print("3. Final Fst Calculation:")
    print("   Fst = (Ht - Hs) / Ht")
    print(f"   Fst = ({Ht} - {Hs}) / {Ht}")
    print(f"   Fst = {Fst}\n")

    print(f"Conclusion: The Fst value is {Fst}, the maximum possible value, indicating complete differentiation.")
    print("This demonstrates that sex-determination systems are a key explanation for this phenomenon.")


# Run the analysis
calculate_fst_for_sex_chromosome()