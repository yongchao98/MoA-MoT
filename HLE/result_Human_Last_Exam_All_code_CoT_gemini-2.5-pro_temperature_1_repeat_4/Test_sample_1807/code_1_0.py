import numpy as np

def calculate_fst(p_males, p_females, n_males, n_females, locus_type='autosomal'):
    """
    Calculates Fst between males and females for a given locus.
    p_males, p_females: Allele frequencies in males and females.
    n_males, n_females: Number of male and female individuals.
    locus_type: 'autosomal' or 'x_linked' (for XY systems).
    """
    if locus_type == 'autosomal':
        # Males and females are both diploid
        alleles_males = 2 * n_males
        alleles_females = 2 * n_females
    elif locus_type == 'x_linked':
        # Males are hemizygous (1 copy), females are diploid (2 copies)
        alleles_males = n_males
        alleles_females = 2 * n_females
    else:
        raise ValueError("locus_type must be 'autosomal' or 'x_linked'")

    # Total allele frequency (p_bar)
    p_total = (p_males * alleles_males + p_females * alleles_females) / (alleles_males + alleles_females)
    
    # Expected heterozygosity in the total population (Ht)
    Ht = 2 * p_total * (1 - p_total)

    # Expected heterozygosities in subpopulations (Hs)
    H_males = 2 * p_males * (1 - p_males)
    H_females = 2 * p_females * (1 - p_females)
    
    # Average expected heterozygosity, weighted by number of individuals
    Hs = (n_males * H_males + n_females * H_females) / (n_males + n_females)
    
    # Fst calculation
    if Ht == 0:
        return 0.0
    fst = (Ht - Hs) / Ht
    
    return fst

# --- Simulation Parameters ---
n_males = 100
n_females = 100
# True allele frequency in the overall gene pool
p_true = 0.6

# --- 1. Autosomal Locus Simulation ---
# In an autosomal locus, both sexes draw from the same gene pool with ploidy 2.
# Due to random sampling, their allele frequencies will be very similar.
p_males_auto = np.random.binomial(2 * n_males, p_true) / (2 * n_males)
p_females_auto = np.random.binomial(2 * n_females, p_true) / (2 * n_females)
fst_autosomal = calculate_fst(p_males_auto, p_females_auto, n_males, n_females, 'autosomal')

# --- 2. X-Linked Locus Simulation (XY system) ---
# Males (XY) have 1 copy of the X chromosome, females (XX) have 2.
# They sample different numbers of alleles from the gene pool.
p_males_x = np.random.binomial(n_males, p_true) / n_males # Sample 100 alleles for 100 males
p_females_x = np.random.binomial(2 * n_females, p_true) / (2 * n_females) # Sample 200 alleles for 100 females
fst_x_linked = calculate_fst(p_males_x, p_females_x, n_males, n_females, 'x_linked')

# --- Print Results ---
print("Demonstrating Genetic Differentiation (Fst) Between Sexes\n")
print("Scenario 1: Autosomal Locus")
print(f"Random sampling resulted in:")
print(f"  - Male allele frequency: {p_males_auto:.4f}")
print(f"  - Female allele frequency: {p_females_auto:.4f}")
print(f"Calculated Fst between sexes = {fst_autosomal:.4f}")
print("This value is typically very close to zero, indicating no differentiation.\n")

print("Scenario 2: X-Linked Locus (in an XY system)")
print(f"Random sampling resulted in:")
print(f"  - Male allele frequency (1 copy/male): {p_males_x:.4f}")
print(f"  - Female allele frequency (2 copies/female): {p_females_x:.4f}")
print(f"Calculated Fst between sexes = {fst_x_linked:.4f}")
print("This value is significantly higher, indicating genetic differentiation caused by the difference in inheritance and chromosome copy number.")
