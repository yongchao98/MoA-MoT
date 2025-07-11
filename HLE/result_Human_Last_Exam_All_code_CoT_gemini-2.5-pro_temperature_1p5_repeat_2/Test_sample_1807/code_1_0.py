import numpy as np

def calculate_and_print_fst(marker_type, p_male, p_female, n_alleles_male, n_alleles_female):
    """
    Calculates and prints Fst between a male and female subpopulation.
    
    Args:
        marker_type (str): A description of the marker (e.g., "Autosomal").
        p_male (float): Allele frequency in males.
        p_female (float): Allele frequency in females.
        n_alleles_male (int): Total number of alleles for this locus in males.
        n_alleles_female (int): Total number of alleles for this locus in females.
    """
    print(f"--- Fst Calculation for {marker_type} Marker ---")
    
    # Calculate H_s: the average expected heterozygosity within subpopulations
    # H = 2 * p * (1-p)
    h_male = 2 * p_male * (1 - p_male)
    h_female = 2 * p_female * (1 - p_female)
    
    # We use a simple unweighted average since male/female group sizes are often equal
    # and this demonstrates the concept clearly.
    hs = (h_male + h_female) / 2
    
    # Calculate H_t: the expected heterozygosity in the total population
    total_alleles = n_alleles_male + n_alleles_female
    if total_alleles == 0:
        print("Cannot calculate Fst, no alleles in the population.\n")
        return

    # p_total is the weighted average allele frequency across subpopulations
    p_total = (p_male * n_alleles_male + p_female * n_alleles_female) / total_alleles
    ht = 2 * p_total * (1 - p_total)
    
    # Calculate Fst = (H_t - H_s) / H_t
    if ht == 0:
        # If Ht is 0, it means the locus is fixed in the total population.
        # Differentiation is considered 0 unless Hs is also non-zero, which is paradoxical.
        fst = 0.0
    else:
        fst = (ht - hs) / ht
        
    print(f"Male allele frequency (p_male): {p_male:.3f}")
    print(f"Female allele frequency (p_female): {p_female:.3f}")
    print(f"Average subpopulation heterozygosity (Hs): {hs:.3f}")
    print(f"Total population heterozygosity (Ht): {ht:.3f}")

    # Output the final equation with numbers
    print("\nFst = (Ht - Hs) / Ht")
    print(f"Fst = ({ht:.3f} - {hs:.3f}) / {ht:.3f} = {fst:.3f}\n")


def main():
    # --- Population Setup ---
    # Assume an XY system with equal numbers of males and females
    N_MALES = 100
    N_FEMALES = 100
    
    # --- Scenario 1: Autosomal Marker ---
    # Allele frequencies should be nearly identical by chance.
    # We simulate this by drawing from the same binomial distribution.
    p_base_autosomal = 0.5
    # Males are diploid (2 alleles), females are diploid (2 alleles)
    n_alleles_male_auto = 2 * N_MALES
    n_alleles_female_auto = 2 * N_FEMALES
    p_male_auto = np.random.binomial(n_alleles_male_auto, p_base_autosomal) / n_alleles_male_auto
    p_female_auto = np.random.binomial(n_alleles_female_auto, p_base_autosomal) / n_alleles_female_auto
    calculate_and_print_fst(
        "Autosomal", p_male_auto, p_female_auto, n_alleles_male_auto, n_alleles_female_auto
    )

    # --- Scenario 2: X-Linked Marker with Divergent Frequencies ---
    # This simulates a case where selection or drift has caused different
    # allele frequencies on the X chromosome between sexes.
    # Males are hemizygous (1 allele), females are diploid (2 alleles)
    n_alleles_male_x = 1 * N_MALES
    n_alleles_female_x = 2 * N_FEMALES
    p_male_x = 0.8  # e.g., allele is advantageous in males
    p_female_x = 0.3 # e.g., allele is less advantageous in females
    calculate_and_print_fst(
        "X-Linked (Divergent)", p_male_x, p_female_x, n_alleles_male_x, n_alleles_female_x
    )
    
    # --- Scenario 3: Y-Linked Marker ---
    # The marker is present on the Y chromosome and absent from the X.
    # Males are hemizygous (1 Y-allele), females have none (0 Y-alleles)
    n_alleles_male_y = 1 * N_MALES
    n_alleles_female_y = 0 * N_FEMALES # No Y chromosome
    p_male_y = 1.0      # All males have the marker
    p_female_y = 0.0    # All females lack the marker
    calculate_and_print_fst(
        "Y-Linked", p_male_y, p_female_y, n_alleles_male_y, n_alleles_female_y
    )

if __name__ == "__main__":
    main()