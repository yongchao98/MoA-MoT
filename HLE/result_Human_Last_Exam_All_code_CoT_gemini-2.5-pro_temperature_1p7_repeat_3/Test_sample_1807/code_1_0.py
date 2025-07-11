import numpy as np

def calculate_allele_freq(genotypes):
    """Calculates allele frequencies from a list of genotypes."""
    # This is a simplified frequency calculator for demonstration.
    # It assumes single-letter alleles.
    alleles = [g for g in genotypes if g != '-']
    if not alleles:
        return {}
    unique_alleles = set(alleles)
    freqs = {allele: alleles.count(allele) / len(alleles) for allele in unique_alleles}
    return freqs

def print_frequencies(marker_name, male_freqs, female_freqs):
    """Prints a formatted comparison of allele frequencies."""
    print(f"\n--- {marker_name} ---")
    print(f"Male Allele Frequencies: {male_freqs}")
    print(f"Female Allele Frequencies: {female_freqs}")
    # A simple qualitative check for differentiation
    if male_freqs != female_freqs:
        print("Result: Allele frequencies differ. High Fst is expected.")
    else:
        print("Result: Allele frequencies are similar. Low Fst is expected.")

def main():
    # Number of individuals in the simulation
    n_individuals = 50

    # --- Case 1: Autosomal Marker (e.g., on Chromosome 1) ---
    # Both males and females have two copies. Allele frequencies should be similar.
    # We'll draw from a common gene pool.
    autosomal_alleles = ['A'] * 100 + ['G'] * 100
    np.random.shuffle(autosomal_alleles)
    
    # We assign two alleles to each individual, simulating diploidy
    male_autosomal_genotypes = [autosomal_alleles.pop() for _ in range(n_individuals * 2)]
    female_autosomal_genotypes = [autosomal_alleles.pop() for _ in range(n_individuals * 2)]

    male_auto_freqs = calculate_allele_freq(male_autosomal_genotypes)
    female_auto_freqs = calculate_allele_freq(female_autosomal_genotypes)
    
    print("This simulation demonstrates why sex determination systems can cause high Fst between sexes.")
    print("We will compare an autosomal marker to a sex-linked marker in an XY system.")
    print_frequencies("Autosomal Marker", male_auto_freqs, female_auto_freqs)


    # --- Case 2: Y-Linked Marker (on the Y chromosome) ---
    # Only males have this marker. We'll represent the female genotype as absent ('-').
    male_y_genotypes = ['Y1'] * n_individuals # All males have an allele on their Y chromosome
    female_y_genotypes = ['-'] * n_individuals # Females do not have a Y chromosome

    male_y_freqs = calculate_allele_freq(male_y_genotypes)
    female_y_freqs = calculate_allele_freq(female_y_genotypes)

    print_frequencies("Y-Linked Marker", male_y_freqs, female_y_freqs)
    
    print("\nExplanation:")
    print("The Autosomal Marker shows similar allele frequencies between males and females because it is inherited by both in the same way. This would lead to a low Fst value.")
    print("The Y-Linked Marker is only present in males. This complete difference in presence/absence leads to maximally different 'allele' frequencies (1.0 vs 0.0), resulting in a very high Fst.")
    print("The same logic applies to W-linked markers in a ZW system. Therefore, sex determination systems are a key explanation for high Fst between sexes.")


if __name__ == "__main__":
    main()
