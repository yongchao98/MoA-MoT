import math

def calculate_fst_variance_based(p_sub1, p_sub2):
    """
    Calculates Fst between two subpopulations based on allele frequency variance.
    This is suitable for cases like Y-linkage where one subpopulation lacks the locus.
    Fst = Var(p) / (p_bar * (1 - p_bar))
    """
    if p_sub1 is None or p_sub2 is None:
        return None # Locus not comparable in this way

    p_bar = (p_sub1 + p_sub2) / 2
    if p_bar == 0 or p_bar == 1:
        # No variation in the total population, Fst is undefined or 0
        return 0.0

    var_p = ((p_sub1 - p_bar)**2 + (p_sub2 - p_bar)**2) / 2
    fst = var_p / (p_bar * (1 - p_bar))
    return fst

def calculate_fst_heterozygosity_based(p_female, p_male):
    """
    Calculates Fst for an X-linked locus between males (hemizygous) and females (diploid).
    Fst = (Ht - Hs) / Ht
    Assumes a 1:1 sex ratio.
    """
    # Hs = average observed heterozygosity in subpopulations
    # For X-linked locus, males are hemizygous, so their heterozygosity is 0.
    h_female = 2 * p_female * (1 - p_female)
    h_male = 0
    hs = (h_female + h_male) / 2

    # Ht = expected heterozygosity in total population
    # Females have 2/3 of X chromosomes, males have 1/3
    p_total = (2/3) * p_female + (1/3) * p_male
    ht = 2 * p_total * (1 - p_total)

    if ht == 0:
        return 0.0

    fst = (ht - hs) / ht
    return fst

print("Demonstrating how sex-linkage affects Fst between males and females.")
print("-" * 70)

# --- Scenario 1: Autosomal Locus ---
# Allele frequencies are expected to be the same in males and females.
p_male_auto = 0.5
p_female_auto = 0.5
fst_auto = calculate_fst_variance_based(p_female_auto, p_male_auto)
print(f"Scenario 1: Autosomal Locus")
print(f"Assuming allele 'A' frequency is {p_female_auto} in females and {p_male_auto} in males.")
print(f"The calculated Fst is: {fst_auto:.4f}")
print("Result: As expected, there is no genetic differentiation (Fst=0).\n")


# --- Scenario 2: Y-linked Locus (in an XY system) ---
# The marker is present in males (freq=1) and absent in females (freq=0).
p_male_y = 1.0
p_female_y = 0.0
fst_y = calculate_fst_variance_based(p_female_y, p_male_y)
print(f"Scenario 2: Y-linked Locus")
print(f"Assuming marker is present in all males (frequency={p_male_y}) and absent in females (frequency={p_female_y}).")
print(f"The calculated Fst is: {fst_y:.4f}")
print("Result: Differentiation is maximum (Fst=1). This would cause a 'pronounced' Fst value.\n")


# --- Scenario 3: X-linked Locus (in an XY system) ---
# Allele frequencies can differ due to drift or sex-specific selection.
# Let's assume a noticeable difference.
p_female_x = 0.8
p_male_x = 0.3
fst_x = calculate_fst_heterozygosity_based(p_female_x, p_male_x)
print(f"Scenario 3: X-linked Locus")
print(f"Assuming allele 'A' frequency is {p_female_x} in females and {p_male_x} in males.")
print(f"The calculated Fst is: {fst_x:.4f}")
print("Result: Significant differentiation (high Fst) occurs due to different allele frequencies on the X chromosome.\n")

print("-" * 70)
print("Conclusion: The most plausible explanation for high Fst at specific markers between sexes is that these markers are on sex chromosomes (Y, W, X, or Z), which have different inheritance patterns and copy numbers between males and females. This corresponds to Answer Choice B.")
