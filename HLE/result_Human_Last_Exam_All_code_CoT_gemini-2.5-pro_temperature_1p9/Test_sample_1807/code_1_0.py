import numpy as np

def calculate_fst(p1, p2):
    """
    Calculates Fst between two populations with allele frequencies p1 and p2.
    Fst = (Ht - Hs) / Ht
    where Ht is expected heterozygosity in the total population and
    Hs is the average expected heterozygosity within subpopulations.
    """
    # Allele frequencies in the two subpopulations (males and females)
    freqs = np.array([p1, p2])

    # Calculate average subpopulation heterozygosity (Hs)
    hs_subpops = 2 * freqs * (1 - freqs)
    hs_average = np.mean(hs_subpops)

    # Calculate total population allele frequency (p_total)
    p_total = np.mean(freqs)
    # Calculate total heterozygosity (Ht)
    ht = 2 * p_total * (1 - p_total)

    # Calculate Fst
    # Avoid division by zero if Ht is 0
    if ht == 0:
        return 0.0
    
    fst = (ht - hs_average) / ht
    return fst

def run_simulation():
    """
    Simulates Fst calculation for an autosomal vs. a sex-linked marker.
    """
    print("--- Fst Simulation: Males vs. Females in a Population ---\n")

    # --- Case 1: Autosomal Marker ---
    # On an autosome, allele frequencies should be nearly identical between sexes
    # in a large, randomly mating population. Let's assume allele 'A' is at 60%.
    print("Scenario 1: Autosomal Marker (e.g., on chromosome 3)")
    p_male_auto = 0.60
    p_female_auto = 0.60
    fst_auto = calculate_fst(p_male_auto, p_female_auto)
    print(f"  Allele 'A' frequency in males:   {p_male_auto:.2f}")
    print(f"  Allele 'A' frequency in females: {p_female_auto:.2f}")
    print(f"  Resulting Fst between sexes:   {fst_auto:.4f}")
    print("  Conclusion: As expected, differentiation is zero (or near zero) for an autosomal marker.\n")


    # --- Case 2: Sex-Linked Marker (in a ZW system) ---
    # Let's consider a marker on the W chromosome. Females are ZW, Males are ZZ.
    # The W chromosome is ONLY found in females.
    # So, the frequency of a W-linked marker allele is 1.0 in females and 0.0 in males.
    print("Scenario 2: Sex-Linked Marker (e.g., on the W chromosome in a bird)")
    p_male_sexlinked = 0.0
    p_female_sexlinked = 1.0
    fst_sexlinked = calculate_fst(p_male_sexlinked, p_female_sexlinked)
    print(f"  Allele 'W1' frequency in males:   {p_male_sexlinked:.2f}")
    print(f"  Allele 'W1' frequency in females: {p_female_sexlinked:.2f}")
    print(f"  Resulting Fst between sexes:   {fst_sexlinked:.4f}")
    print("  Conclusion: Differentiation is maximal (1.0) for a sex-specific marker.\n")
    
    print("This simulation demonstrates how sex-chromosome systems directly cause high Fst between sexes for linked markers.")


if __name__ == "__main__":
    run_simulation()
<<<B>>>