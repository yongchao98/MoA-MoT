import numpy as np

def calculate_fst(p_female, p_male, weight_female, weight_male):
    """
    Calculates Fst between two subpopulations (males and females)
    using the variance of allele frequencies.
    Fst = Var(p) / (p_bar * (1 - p_bar))

    Args:
        p_female (float): Allele frequency in females.
        p_male (float): Allele frequency in males.
        weight_female (float): Proportional contribution of females to the total gene pool for that marker.
        weight_male (float): Proportional contribution of males to the total gene pool for that marker.

    Returns:
        float: The calculated Fst value.
    """
    # Calculate the average allele frequency across the total population, weighted by contribution
    p_bar = (p_female * weight_female) + (p_male * weight_male)

    # Avoid division by zero if an allele is fixed (p_bar is 0 or 1)
    if p_bar == 0 or p_bar == 1:
        return 0.0

    # Calculate the weighted variance of allele frequencies between the two groups
    var_p = weight_female * (p_female - p_bar)**2 + weight_male * (p_male - p_bar)**2

    # Calculate Fst
    fst = var_p / (p_bar * (1 - p_bar))
    return fst

def explain_sex_linked_fst():
    """
    Explains and demonstrates why sex-linkage can cause high Fst between sexes.
    """
    print("This script demonstrates why a sex-determination system (like XY vs. XX) is a potential explanation for high Fst between males and females at some genetic markers.\n")
    
    # --- Case 1: Autosomal Marker ---
    print("--- 1. Autosomal Marker (Not on a sex chromosome) ---")
    print("For autosomal markers, both males and females have two copies and contribute equally to the gene pool (assuming a 1:1 sex ratio). Allele frequencies are expected to be nearly identical.")
    
    p_female_auto = 0.50
    p_male_auto = 0.51  # Assume a small difference due to random sampling
    # Weights are 0.5 for each sex because they contribute equally (2 chromosomes each)
    w_female_auto = 0.5
    w_male_auto = 0.5
    
    fst_auto = calculate_fst(p_female_auto, p_male_auto, w_female_auto, w_male_auto)
    
    print(f"Female allele frequency (p_female): {p_female_auto}")
    print(f"Male allele frequency (p_male):     {p_male_auto}")
    print(f"Calculated Fst for the autosomal marker: {fst_auto:.5f}")
    print("This very low Fst value shows negligible differentiation, as expected.\n")

    # --- Case 2: X-linked Marker ---
    print("--- 2. X-linked Marker (On a sex chromosome) ---")
    print("For X-linked markers, females (XX) have two copies but males (XY) have only one. Therefore, females contribute 2/3 and males 1/3 of the X chromosomes to the population.")
    print("This difference, especially under sex-specific selection, can lead to large differences in allele frequencies.")

    # Assume a scenario with strong divergence in allele frequencies
    p_female_x = 0.80
    p_male_x = 0.20
    
    # Weights are 2/3 for females and 1/3 for males based on their chromosome count
    w_female_x = 2/3
    w_male_x = 1/3

    fst_x = calculate_fst(p_female_x, p_male_x, w_female_x, w_male_x)

    print(f"Female allele frequency (p_female): {p_female_x}")
    print(f"Male allele frequency (p_male):     {p_male_x}")
    print(f"Calculated Fst for the X-linked marker: {fst_x:.5f}")
    print("This high Fst value shows pronounced differentiation.\n")
    
    print("--- Conclusion ---")
    print("The presence of a sex-determining system (like XY or ZW) means markers on sex chromosomes have different inheritance patterns and ploidy between sexes. This can cause their allele frequencies to diverge, resulting in the high Fst values observed in the study, while autosomal markers would not show such a pattern.")

# Run the explanation and simulation
explain_sex_linked_fst()
