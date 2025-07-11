import sys

def calculate_fst_sex_linked():
    """
    Calculates and explains the theoretical Fst between sexes for an X-linked marker.

    This function demonstrates that the difference in ploidy between males (XY) and
    females (XX) for sex chromosomes inherently leads to a non-zero Fst, indicating
    genetic differentiation.
    """
    print("This script demonstrates why sex-determination systems (like XY) can cause high Fst between males and females.")
    print("It calculates the theoretical Fst for a marker on the X chromosome in a population at equilibrium.\n")

    # Plan explanation
    print("Plan:")
    print("1. Assume a population with an X-linked gene having two alleles, 'A' and 'a'.")
    print("2. Set the frequency of allele 'A' to a value 'p'. The frequency of 'a' is 'q = 1 - p'.")
    print("3. Calculate Ht (Total Heterozygosity): This is the expected heterozygosity in the total gene pool of X chromosomes. Assuming the population is in equilibrium, the overall frequency of 'A' is 'p', so Ht = 2*p*q.")
    print("4. Calculate Hs (Average within-subpopulation Heterozygosity): This is the average of the heterozygosity found within females and within males.")
    print("   - Females (XX) can be heterozygous, so H_females = 2*p*q.")
    print("   - Males (XY) are hemizygous (only one X), so they cannot be heterozygous. H_males = 0.")
    print("   - Hs is the simple average: (H_females + H_males) / 2.")
    print("5. Calculate Fst using the standard formula: Fst = (Ht - Hs) / Ht.")
    print("-" * 20)

    # --- Calculation ---

    # Define the allele frequency 'p' for allele 'A'.
    # A value between 0 and 1 (exclusive) will show differentiation.
    try:
        p = 0.7
        if not (0 < p < 1):
            raise ValueError
    except (ValueError, IndexError):
        print("Using default allele frequency p = 0.7")
        p = 0.7

    q = 1 - p

    # Ht: Expected heterozygosity in the total population gene pool (all X chromosomes)
    # The overall allele frequency is p, so Ht = 2*p*q
    Ht = 2 * p * q

    # Hs: Average heterozygosity within each subpopulation (sex)
    H_females = 2 * p * q
    H_males = 0
    Hs = (H_females + H_males) / 2

    # Fst calculation
    # Handle the case where Ht is 0 (no variation in the population)
    if Ht > 0:
        Fst = (Ht - Hs) / Ht
    else:
        Fst = 0

    # Print the results, including the equation with numbers
    print("\n--- Results ---")
    print(f"Given allele frequencies p = {p:.2f} and q = {q:.2f}:\n")
    print(f"1. Total Expected Heterozygosity (Ht):")
    print(f"   Ht = 2 * p * q")
    print(f"   Ht = 2 * {p:.2f} * {q:.2f} = {Ht:.4f}\n")

    print(f"2. Average Within-Sex Heterozygosity (Hs):")
    print(f"   Hs = (H_females + H_males) / 2")
    print(f"   Hs = ( (2 * {p:.2f} * {q:.2f}) + {H_males} ) / 2 = {Hs:.4f}\n")

    print(f"3. Final Fst Calculation:")
    print(f"   Fst = (Ht - Hs) / Ht")
    print(f"   Fst = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f}")
    print(f"   Fst = {Fst:.4f}\n")
    
    print("This high Fst value (0.5) is not due to population subdivision or selection in this model,")
    print("but arises directly from the difference in sex chromosome composition between males and females.")

if __name__ == "__main__":
    calculate_fst_sex_linked()
