import numpy as np

# This script simulates gene flow between two populations to determine the genetic consequences.

# --- Helper functions for genetic calculations ---

def get_allele_freqs(population):
    """Calculates allele frequencies (p for 'A', q for 'a') from a list of genotypes."""
    n_total = len(population)
    if n_total == 0: return 0.0, 0.0
    
    n_A = sum(genotype.count('A') for genotype in population)
    n_a = sum(genotype.count('a') for genotype in population)
    total_alleles = n_A + n_a
    if total_alleles == 0: return 0.0, 0.0
    
    p = n_A / total_alleles
    q = n_a / total_alleles
    return p, q

def get_heterozygosity(population):
    """Calculates observed (Ho) and expected (Hs) heterozygosity."""
    n_total = len(population)
    if n_total == 0: return 0.0, 0.0
    
    # Observed heterozygosity is the actual frequency of heterozygotes
    n_het = population.count('Aa')
    Ho = n_het / n_total

    # Expected heterozygosity under Hardy-Weinberg Equilibrium
    p, q = get_allele_freqs(population)
    Hs = 2 * p * q
    
    return Ho, Hs

# --- Functions to calculate metrics from the answer choices ---

def calculate_pi(population):
    """Calculate nucleotide diversity (Pi). For a single locus, this is the expected heterozygosity."""
    _, Hs = get_heterozygosity(population)
    return Hs

def calculate_fis(population):
    """Calculate the inbreeding coefficient (Fis = (Hs - Ho) / Hs)."""
    Ho, Hs = get_heterozygosity(population)
    if Hs == 0:
        # Occurs in pure parental populations where no heterozygosity is expected. Fis is undefined.
        return 0.0 
    return (Hs - Ho) / Hs

def calculate_fst(pop1, pop2):
    """Calculate Fst between two populations."""
    _, Hs1 = get_heterozygosity(pop1)
    _, Hs2 = get_heterozygosity(pop2)
    Hs_avg = np.mean([Hs1, Hs2])
    
    # Calculate Ht from the total, pooled population
    _, Ht = get_heterozygosity(pop1 + pop2)
    if Ht == 0: return 0.0
        
    return (Ht - Hs_avg) / Ht

# --- Main Simulation and Analysis ---

# 1. Define two distinct parent populations fixed for different alleles
parent_pop1 = ['AA'] * 100
parent_pop2 = ['aa'] * 100

# 2. Simulate gene flow creating an F1 hybrid population
hybrid_pop_f1 = ['Aa'] * 100 

print("--- Analysis of Gene Flow Across a Hybrid Zone ---\n")
print(f"We model two parent populations, P1 (all 'AA') and P2 (all 'aa').")
print(f"Gene flow between them creates a hybrid population of F1 offspring (all 'Aa').\n")

# A. Analysis of Fst
fst_parental = calculate_fst(parent_pop1, parent_pop2)
print("A. Can HIGH Fst between populations occur?")
print(f"   - The Fst between the distinct parent populations P1 and P2 is: {fst_parental:.2f}")
print("   -> Yes. A hybrid zone forms between populations that are already highly differentiated (High Fst).\n")

# B. Analysis of Dxy
print("B. Can HIGH Dxy between populations occur?")
print("   - Dxy measures absolute divergence reflecting historical split time. It is not erased by recent gene flow.")
print("   -> Yes. If the parent populations have been separated for a long time, their Dxy will be high.\n")

# E. Analysis of Pi
pi_hybrid = calculate_pi(hybrid_pop_f1)
print("E. Can HIGH Pi within a population occur?")
print(f"   - Pi in the parent populations is 0, but in the new hybrid population, Pi is: {pi_hybrid:.2f}")
print("   -> Yes. Gene flow introduces new alleles, increasing genetic variation (Pi) in the hybrid zone.\n")

# D. Analysis of mutation rate (μ)
print("D. Can HIGH μ (mutation rate) occur?")
print("   - The mutation rate is a biological property of DNA replication and is independent of gene flow, which only moves existing alleles.")
print("   -> Yes. A population can have a high mutation rate for reasons unrelated to gene flow.\n")

# C. Analysis of Fis
fis_hybrid = calculate_fis(hybrid_pop_f1)
# For the Fis equation: Fis = (Hs - Ho) / Hs
Ho_hybrid, Hs_hybrid = get_heterozygosity(hybrid_pop_f1)
print("C. Can HIGH Fis within a population occur?")
print("   - Gene flow (mating between P1 and P2) creates heterozygotes ('Aa') exclusively.")
print(f"   - In the F1 hybrid population, expected heterozygosity (Hs) is {Hs_hybrid:.2f} and observed heterozygosity (Ho) is {Ho_hybrid:.2f}.")
print(f"   - The resulting Fis = ({Hs_hybrid:.2f} - {Ho_hybrid:.2f}) / {Hs_hybrid:.2f} = {fis_hybrid:.2f}")
print("   - A high POSITIVE Fis means a deficit of heterozygotes (inbreeding). Gene flow does the opposite, creating an EXCESS of heterozygotes, which results in a NEGATIVE Fis.")
print("   -> No. This CANNOT occur as a direct result of gene flow.\n")
