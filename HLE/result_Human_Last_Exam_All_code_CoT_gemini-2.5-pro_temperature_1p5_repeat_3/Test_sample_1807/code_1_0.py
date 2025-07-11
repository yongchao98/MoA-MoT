def calculate_fst(p_group1, p_group2):
    """
    Calculates a simple Fst between two groups for a bi-allelic marker.
    Fst = (Ht - Hs) / Ht
    where Ht is the expected heterozygosity in the total population and
    Hs is the average expected heterozygosity within the subgroups.

    Args:
    p_group1 (float): Frequency of an allele in group 1 (e.g., males).
    p_group2 (float): Frequency of the same allele in group 2 (e.g., females).

    Returns:
    float: The Fst value.
    """
    # If frequencies are identical, there is no differentiation.
    if p_group1 == p_group2:
        return 0.0, "p_males and p_females are equal"

    # Calculate average allele frequency across the two groups
    p_total = (p_group1 + p_group2) / 2.0

    # Calculate total expected heterozygosity (Ht)
    # Ht = 2 * p_total * (1 - p_total)
    ht = 2 * p_total * (1 - p_total)

    # Calculate expected heterozygosity for each group
    h_group1 = 2 * p_group1 * (1 - p_group1)
    h_group2 = 2 * p_group2 * (1 - p_group2)

    # Calculate the average subpopulation heterozygosity (Hs)
    hs = (h_group1 + h_group2) / 2.0
    
    # If Ht is 0, Fst is undefined, but this case is handled by the initial check.
    if ht == 0:
      return float('nan'), "Ht is zero"

    # Calculate Fst
    fst = (ht - hs) / ht
    
    # Return Fst and the equation components for clarity
    equation = f"Fst = (Ht - Hs) / Ht = ({ht:.2f} - {hs:.2f}) / {ht:.2f}"
    return fst, equation


# --- Scenario 1: Autosomal Marker ---
# In a randomly mating population, allele frequencies at autosomal loci
# are expected to be the same between males and females.
p_auto_males = 0.5
p_auto_females = 0.5
fst_autosomal, eq_autosomal = calculate_fst(p_auto_males, p_auto_females)

print("--- Scenario 1: Autosomal Marker ---")
print(f"Assumed allele frequency in Males: {p_auto_males}")
print(f"Assumed allele frequency in Females: {p_auto_females}")
print(f"Resulting Fst = {fst_autosomal:.4f}")
print("Explanation: Because allele frequencies are identical, there is no genetic differentiation.\n")


# --- Scenario 2: Y-Linked Marker (XY System) ---
# A marker on the Y chromosome is present only in males.
# Its frequency is 1.0 in males and 0.0 in females.
p_y_linked_males = 1.0
p_y_linked_females = 0.0
fst_y_linked, eq_y_linked = calculate_fst(p_y_linked_males, p_y_linked_females)

print("--- Scenario 2: Sex-Linked Marker (e.g., Y-Linked) ---")
print(f"Allele frequency in Males: {p_y_linked_males}")
print(f"Allele frequency in Females: {p_y_linked_females}")
print("Calculation: " + eq_y_linked)
print(f"Resulting Fst = {fst_y_linked:.4f}")
print("Explanation: The allele is fixed in one sex and absent in the other, leading to maximum differentiation (Fst=1).\n")

print("Conclusion: The pronounced differentiation at some markers is best explained by them being on sex chromosomes.")
