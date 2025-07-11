import numpy as np

def calculate_fst(freqs, sizes):
    """
    Calculates Fst from allele frequencies and sizes of subpopulations.
    Fst = (Ht - Hs) / Ht
    where Ht is heterozygosity of the total population and
    Hs is the average heterozygosity of the subpopulations.
    """
    # Ensure freqs and sizes are numpy arrays
    freqs = np.array(freqs)
    sizes = np.array(sizes)
    
    # Calculate the average allele frequency in the total population (p_bar)
    p_bar = np.average(freqs, weights=sizes)
    
    # Calculate Ht (expected heterozygosity in the total population)
    ht = 2 * p_bar * (1 - p_bar)
    
    # Calculate Hs (average expected heterozygosity within subpopulations)
    h_subpops = 2 * freqs * (1 - freqs)
    hs = np.average(h_subpops, weights=sizes)
    
    # Fst is the proportion of total variance found between populations
    # Handle the case where Ht is zero to avoid division by zero
    if ht == 0:
        return 0.0
    
    fst = (ht - hs) / ht
    return fst

# --- Population Setup ---
# Let's assume a population with 100 males and 100 females
pop_sizes = [100, 100] # [males, females]

# --- Case 1: Autosomal Marker (on a non-sex chromosome) ---
# For a typical autosomal marker, males and females draw from the same gene pool.
# We expect their allele frequencies to be nearly identical.
# Let's assume the frequency of allele 'A' is 0.7 in both groups.
allele_freqs_autosomal = [0.7, 0.7] # [males, females]
fst_autosomal = calculate_fst(allele_freqs_autosomal, pop_sizes)

print("--- Case 1: Autosomal Marker ---")
print(f"Allele frequency in males: {allele_freqs_autosomal[0]}")
print(f"Allele frequency in females: {allele_freqs_autosomal[1]}")
print(f"Calculated Fst between sexes: {fst_autosomal:.4f}")
print("Result: Fst is 0, indicating no genetic differentiation.\n")


# --- Case 2: Sex-Linked Marker (on the Y-chromosome in an XY system) ---
# A marker on the Y chromosome is present only in males.
# The frequency of the marker allele is 1.0 in males and 0.0 in females.
allele_freqs_sex_linked = [1.0, 0.0] # [males, females]
fst_sex_linked = calculate_fst(allele_freqs_sex_linked, pop_sizes)

print("--- Case 2: Y-Linked Marker ---")
print(f"Allele frequency in males: {allele_freqs_sex_linked[0]}")
print(f"Allele frequency in females: {allele_freqs_sex_linked[1]}")
print(f"Calculated Fst between sexes: {fst_sex_linked:.4f}")
print("Result: Fst is 1.0, indicating complete genetic differentiation.")
print("\nConclusion: The simulation shows that markers on sex chromosomes can cause pronounced genetic differentiation between males and females.")
