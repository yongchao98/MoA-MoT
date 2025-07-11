import numpy as np

def calculate_fst(pop1_alleles, pop2_alleles):
    """
    A simplified Fst calculation for a single bi-allelic locus.
    Assumes diploid individuals.
    pop_alleles is a list of genotypes, where 'A' and 'a' are the alleles.
    e.g., ['AA', 'Aa', 'aa']
    """
    def get_freq(alleles):
        total_alleles = len(alleles) * 2
        p = (alleles.count('AA') * 2 + alleles.count('Aa')) / total_alleles
        q = 1 - p
        return p, q

    # Frequencies in subpopulations
    p1, q1 = get_freq(pop1_alleles)
    p2, q2 = get_freq(pop2_alleles)

    # Average frequency in subpopulations
    p_bar = (p1 + p2) / 2
    q_bar = (q1 + q2) / 2

    # Expected heterozygosity in subpopulations
    hs = (2 * p1 * q1 + 2 * p2 * q2) / 2

    # Expected heterozygosity in total population
    ht = 2 * p_bar * q_bar

    # Fst calculation
    if ht == 0:
        return 0
    fst = (ht - hs) / ht
    return fst

# --- Simulation ---

# Scenario 1: Autosomal Marker
# Alleles are distributed randomly between sexes in the same gene pool.
females_auto = ['AA', 'Aa', 'Aa', 'aa', 'AA', 'Aa']
males_auto   = ['AA', 'Aa', 'aa', 'aa', 'AA', 'Aa']
fst_autosomal = calculate_fst(females_auto, males_auto)
print("Scenario: Autosomal Marker (e.g., on Chromosome 3)")
print(f"Randomly distributed alleles between males and females from the same gene pool.")
print(f"Calculated Fst between males and females: {fst_autosomal:.4f}\n")


# Scenario 2: Y-linked Marker (XY system)
# The marker 'Y' is only present in males. We represent its absence as '0'.
# For calculation, let's treat the 'Y' allele as 'A' and absence as 'a'.
# A male is 'AY' (hemizygous), let's model this as 'Aa' for freq calculation simplicity.
# A female is 'XX' (no Y), let's model this as 'aa'.
females_Y_linked = ['aa', 'aa', 'aa', 'aa', 'aa', 'aa'] # No Y chromosome
males_Y_linked   = ['Aa', 'Aa', 'Aa', 'Aa', 'Aa', 'Aa'] # All have one Y chromosome
fst_y_linked = calculate_fst(females_Y_linked, males_Y_linked)
print("Scenario: Y-linked Marker (in an XY system)")
print(f"The marker is only present on the Y chromosome.")
print(f"Allele 'A' (presence of marker) has frequency 0 in females and 0.5 in males.")
print(f"Calculated Fst between males and females: {fst_y_linked:.4f}\n")


# Scenario 3: W-linked Marker (ZW system)
# The marker 'W' is only present in females.
# A male is 'ZZ' (no W), let's model this as 'aa'.
# A female is 'ZW' (hemizygous), let's model this as 'Aa'.
females_W_linked = ['Aa', 'Aa', 'Aa', 'Aa', 'Aa', 'Aa'] # All have one W chromosome
males_W_linked   = ['aa', 'aa', 'aa', 'aa', 'aa', 'aa'] # No W chromosome
fst_w_linked = calculate_fst(females_W_linked, males_W_linked)
print("Scenario: W-linked Marker (in a ZW system)")
print(f"The marker is only present on the W chromosome.")
print(f"Allele 'A' (presence of marker) has frequency 0.5 in females and 0 in males.")
print(f"Calculated Fst between males and females: {fst_w_linked:.4f}\n")

print("Conclusion: As the code demonstrates, markers on sex chromosomes (Y or W) lead to very high Fst values")
print("when comparing males and females, while autosomal markers show Fst values near zero.")
print("This makes sex chromosome systems the best explanation for the observation.")