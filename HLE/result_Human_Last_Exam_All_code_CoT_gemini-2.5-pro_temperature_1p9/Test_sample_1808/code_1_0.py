import sys

# Suppress writing compiled bytecode files
sys.dont_write_bytecode = True

def calculate_fst(p1, p2, population_size_1=100, population_size_2=100):
    """
    Calculates Fst between two populations given their allele frequencies for a single locus.
    
    This function demonstrates the Fst calculation (Ht - Hs) / Ht.
    p1: Frequency of allele '1' in Population 1.
    p2: Frequency of allele '1' in Population 2.
    """
    # Allele frequencies are p (frequency of allele '1') and q (frequency of allele '0')
    q1 = 1 - p1
    q2 = 1 - p2

    print("--- Calculating Hs (Average heterozygosity in subpopulations) ---")
    # Calculate expected heterozygosity within each subpopulation (H = 2pq)
    h1 = 2 * p1 * q1
    h2 = 2 * p2 * q2
    # Calculate the average of the subpopulations' heterozygosities
    hs = (h1 + h2) / 2
    
    print(f"Expected H in Pop 1 (Hs1) = 2 * {p1:.2f} * {q1:.2f} = {h1:.4f}")
    print(f"Expected H in Pop 2 (Hs2) = 2 * {p2:.2f} * {q2:.2f} = {h2:.4f}")
    print(f"Average Hs = ({h1:.4f} + {h2:.4f}) / 2 = {hs:.4f}")
    

    print("\n--- Calculating Ht (Heterozygosity in total population) ---")
    # Calculate the total frequency of allele '1' across both populations
    total_individuals = population_size_1 + population_size_2
    p_total = (p1 * population_size_1 + p2 * population_size_2) / total_individuals
    q_total = 1 - p_total
    
    # Calculate the expected heterozygosity in the total combined population
    ht = 2 * p_total * q_total
    print(f"Total allele '1' frequency (p_total) = ({p1:.2f}*{population_size_1} + {p2:.2f}*{population_size_2}) / {total_individuals} = {p_total:.4f}")
    print(f"Ht = 2 * {p_total:.4f} * {q_total:.4f} = {ht:.4f}")


    print("\n--- Calculating Final Fst ---")
    # The Fst formula
    # Handle the case of no variation in the total population (ht=0)
    if ht == 0:
        fst = 0.0 # If Ht is 0, there is no variation to partition, so Fst is 0.
    else:
        fst = (ht - hs) / ht
    
    print(f"Fst = (Ht - Hs) / Ht")
    print(f"Fst = ({ht:.4f} - {hs:.4f}) / {ht:.4f} = {fst:.4f}\n")
    print("="*50)


print("Scenario 1: No Gene Flow (Complete Divergence)\n")
# In this case, Population 1 has only allele '0' (p1=0), and Population 2 has only allele '1' (p2=1).
# This results in maximum differentiation.
calculate_fst(p1=0.0, p2=1.0)


print("Scenario 2: Some Gene Flow (Partial Convergence)\n")
# Here, gene flow has started to mix the alleles. Pop 1 now has some '1's and Pop 2 has some '0's.
calculate_fst(p1=0.2, p2=0.8)


print("Scenario 3: High Gene Flow (Near Homogenization)\n")
# After significant gene flow, the allele frequencies in the two populations have become very similar.
calculate_fst(p1=0.45, p2=0.55)
<<<A>>>