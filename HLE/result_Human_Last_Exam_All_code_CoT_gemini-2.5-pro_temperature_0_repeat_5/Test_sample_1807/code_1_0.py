def calculate_fst_for_x_linked_gene(p_A):
    """
    Calculates and explains the Fst value between males and females for an X-linked gene.

    Args:
        p_A (float): The frequency of allele 'A' for the X-linked gene in the total population.
                     Must be between 0 and 1.
    """
    if not 0 <= p_A <= 1:
        print("Error: Allele frequency p_A must be between 0 and 1.")
        return

    p_a = 1 - p_A

    print(f"--- Fst Calculation for an X-linked Gene ---")
    print(f"Assume a population with an XY sex-determination system and equal numbers of males and females.")
    print(f"Let's analyze a gene on the X chromosome with two alleles, 'A' and 'a'.")
    print(f"The overall frequency of allele 'A' in the population is p = {p_A}")
    print(f"The overall frequency of allele 'a' is q = {p_a:.2f}\n")

    # Fst = (Ht - Hs) / Ht
    # Ht: Expected heterozygosity in the total population's gene pool.
    # The total gene pool for the X chromosome consists of all X chromosomes from both
    # males and females. The allele frequency is p_A.
    Ht = 2 * p_A * p_a
    print("Step 1: Calculate Ht (Expected heterozygosity in the total gene pool)")
    print(f"Ht = 2 * p * q = 2 * {p_A} * {p_a:.2f} = {Ht:.4f}\n")

    # Hs: Average expected heterozygosity within the subpopulations (males and females).
    # Hs_female: Females are XX (diploid), so their expected heterozygosity is 2*p*q.
    Hs_female = 2 * p_A * p_a
    # Hs_male: Males are XY (hemizygous), they have only one X. They cannot be heterozygous.
    Hs_male = 0
    # Hs is the average of the two subpopulations.
    Hs = (Hs_female + Hs_male) / 2

    print("Step 2: Calculate Hs (Average heterozygosity within subpopulations)")
    print(f"   - For Females (XX), they are diploid. Hs_female = 2 * p * q = {Hs_female:.4f}")
    print(f"   - For Males (XY), they are hemizygous. Hs_male = {Hs_male}")
    print(f"Hs = (Hs_female + Hs_male) / 2 = ({Hs_female:.4f} + {Hs_male}) / 2 = {Hs:.4f}\n")

    # Calculate Fst
    if Ht == 0:
        Fst = 0 # No variation in the population
        print("Step 3: Calculate Fst. Since Ht is 0, Fst is 0.")
    else:
        Fst = (Ht - Hs) / Ht
        print("Step 3: Calculate Fst = (Ht - Hs) / Ht")
        print(f"Fst = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f} = {Fst:.4f}\n")

    print(f"Result: The calculated Fst value is {Fst:.4f}.")
    print("This is a high value for differentiation, caused simply by the difference in sex chromosome systems, not by any environmental or geographic factor.")

# --- Main execution ---
# Use an example allele frequency of 0.7 for allele 'A'.
allele_frequency_A = 0.7
calculate_fst_for_x_linked_gene(allele_frequency_A)