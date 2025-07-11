import numpy as np

def calculate_fst(p_subpops):
    """
    Calculates Fst from a list of subpopulation allele frequencies.
    Uses the formula Fst = Var(p) / (p_bar * (1 - p_bar)),
    which is equivalent to Wright's Fst.
    """
    p_array = np.array(p_subpops)
    p_bar = np.mean(p_array) # Mean allele frequency across subpopulations
    var_p = np.var(p_array, ddof=0) # Variance in allele frequencies

    # Avoid division by zero if p_bar is 0 or 1 (marker is not polymorphic in total pop)
    if p_bar == 0 or p_bar == 1:
        return 0.0

    fst = var_p / (p_bar * (1 - p_bar))
    return fst

def explain_sex_linked_fst():
    """
    Demonstrates why sex-linked markers cause high Fst between sexes.
    """
    print("This simulation will compare Fst between males and females for two types of genetic markers.\n")

    # --- Case 1: Autosomal Marker ---
    # Autosomal markers are on non-sex chromosomes. Allele frequencies should be similar
    # between males and females, with minor differences due to random chance.
    p_male_auto = 0.51
    p_female_auto = 0.49
    fst_auto = calculate_fst([p_male_auto, p_female_auto])

    print("--- Scenario 1: Autosomal Marker ---")
    print("Marker is on a non-sex chromosome.")
    print(f"Allele frequency in Males: {p_male_auto}")
    print(f"Allele frequency in Females: {p_female_auto}")
    print(f"Calculated Fst between sexes: {fst_auto:.4f}")
    print("Result: Fst is very low, indicating almost no genetic differentiation, as expected.\n")


    # --- Case 2: Sex-Linked Marker (e.g., Y-linked in an XY system) ---
    # In an XY system, a Y-linked marker is present only in males.
    # We can model this with an allele frequency of 1.0 in males and 0.0 in females.
    # The result would be identical for a W-linked marker in a ZW system (0.0 in males, 1.0 in females).
    p_male_sex_linked = 1.0
    p_female_sex_linked = 0.0
    fst_sex_linked = calculate_fst([p_male_sex_linked, p_female_sex_linked])

    print("--- Scenario 2: Sex-Linked Marker (e.g., Y-linked) ---")
    print("Marker is on a sex-specific chromosome (e.g., Y chromosome).")
    print(f"Allele frequency in Males: {p_male_sex_linked}")
    print(f"Allele frequency in Females: {p_female_sex_linked}")
    print(f"Calculated Fst between sexes: {fst_sex_linked:.4f}")
    print("Result: Fst is at its maximum value of 1, indicating complete differentiation.")
    print("\nConclusion: The presence of sex chromosomes (XY or ZW systems) is a direct cause of high Fst at specific loci when comparing males and females.")

if __name__ == '__main__':
    explain_sex_linked_fst()