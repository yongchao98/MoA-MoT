import numpy as np

def calculate_fst(p1, p2):
    """
    Calculates a simple Fst for two populations with given allele frequencies.
    Assumes equal population sizes.
    """
    p_bar = (p1 + p2) / 2
    
    # Handle cases where p_bar is 0 or 1 to avoid division by zero
    if p_bar == 0 or p_bar == 1:
        return 0.0
        
    # Expected heterozygosity in the total population
    ht = 2 * p_bar * (1 - p_bar)
    
    # Expected heterozygosities in each subpopulation
    hs1 = 2 * p1 * (1 - p1)
    hs2 = 2 * p2 * (1 - p2)
    
    # Average expected heterozygosity within subpopulations
    hs = (hs1 + hs2) / 2
    
    fst = (ht - hs) / ht
    return fst

# --- Scenario 1: Autosomal Marker ---
# Both males and females inherit autosomes in the same way.
# Let's assume an allele 'A' has a frequency of 0.6 in the population.
# Due to random mating, the allele frequencies should be very similar in both sexes.
p_male_autosomal = 0.6
p_female_autosomal = 0.6
fst_autosomal = calculate_fst(p_male_autosomal, p_female_autosomal)

print("Scenario 1: Autosomal Marker")
print(f"Allele frequency in Males: {p_male_autosomal}")
print(f"Allele frequency in Females: {p_female_autosomal}")
print(f"Calculated Fst between sexes: {fst_autosomal:.4f}\n")


# --- Scenario 2: Y-linked Marker ---
# This marker is only on the Y chromosome, so only males have it.
# We will define allele 'A' as the presence of this marker.
# All males have the marker, no females have it.
p_male_Y_linked = 1.0  # All males have the marker allele
p_female_Y_linked = 0.0 # No females have the marker allele
fst_Y_linked = calculate_fst(p_male_Y_linked, p_female_Y_linked)

print("Scenario 2: Y-linked Marker")
print(f"Allele frequency in Males: {p_male_Y_linked}")
print(f"Allele frequency in Females: {p_female_Y_linked}")
print(f"Calculated Fst between sexes: {fst_Y_linked:.4f}\n")

print("Conclusion: The high Fst for certain markers is best explained by them being on sex chromosomes (like the Y chromosome), which differ between males and females.")
