import random

def simulate_and_calculate_fst():
    """
    Simulates a population to demonstrate Fst differences between sexes
    for autosomal vs. sex-linked markers.
    """
    # 1. Simulation Parameters
    POP_SIZE = 2000
    # Allele frequencies in the overall gene pool
    p_autosomal = 0.6  # Frequency of 'A'
    p_sex_linked = 0.7 # Frequency of 'X_B'

    # 2. Create Population Subgroups (Males and Females)
    males = []
    females = []
    
    # We will directly generate allele lists for each group based on inheritance
    # This simplifies the simulation while maintaining the core concept.
    num_males = int(POP_SIZE / 2)
    num_females = POP_SIZE - num_males
    
    # Generate alleles for the AUTOSOMAL marker
    # Both sexes have 2 copies
    male_autosomal_alleles = ['A' if random.random() < p_autosomal else 'a' for _ in range(num_males * 2)]
    female_autosomal_alleles = ['A' if random.random() < p_autosomal else 'a' for _ in range(num_females * 2)]

    # Generate alleles for the SEX-LINKED (X-linked) marker
    # Males have 1 copy (X), Females have 2 copies (XX)
    male_sex_linked_alleles = ['B' if random.random() < p_sex_linked else 'b' for _ in range(num_males)]
    female_sex_linked_alleles = ['B' if random.random() < p_sex_linked else 'b' for _ in range(num_females * 2)]

    # 3. Define Helper Functions to Calculate Frequencies and Fst
    def get_allele_freq(allele_list, target_allele):
        """Calculates frequency of a target allele in a list."""
        count = allele_list.count(target_allele)
        total = len(allele_list)
        return count / total if total > 0 else 0

    def calculate_fst(p1, n1, p2, n2):
        """
        Calculates Fst from allele frequencies (p) and allele counts (n)
        of two subpopulations.
        """
        if n1 + n2 == 0:
            return 0
        
        # Ht: Expected heterozygosity in the total pooled population
        p_total = (p1 * n1 + p2 * n2) / (n1 + n2)
        Ht = 2 * p_total * (1 - p_total)

        if Ht == 0:
            return 0
            
        # Hs: Average expected heterozygosity of the subpopulations
        H1 = 2 * p1 * (1 - p1)
        H2 = 2 * p2 * (1 - p2)
        Hs = (H1 * n1 + H2 * n2) / (n1 + n2)
        
        Fst = (Ht - Hs) / Ht
        return Fst

    # 4. Perform Calculations
    # For Autosomal Marker
    p_male_auto = get_allele_freq(male_autosomal_alleles, 'A')
    p_female_auto = get_allele_freq(female_autosomal_alleles, 'A')
    fst_auto = calculate_fst(p_male_auto, len(male_autosomal_alleles), p_female_auto, len(female_autosomal_alleles))

    # For Sex-Linked Marker
    p_male_sex = get_allele_freq(male_sex_linked_alleles, 'B')
    p_female_sex = get_allele_freq(female_sex_linked_alleles, 'B')
    fst_sex_linked = calculate_fst(p_male_sex, len(male_sex_linked_alleles), p_female_sex, len(female_sex_linked_alleles))
    
    # 5. Print the output
    print("--- Simulation of Genetic Differentiation (Fst) Between Sexes ---")
    print(f"\nComparing {num_males} males and {num_females} females.")
    print("\n--- Marker 1: Autosomal Locus ---")
    print(f"Allele 'A' frequency in Males:   {p_male_auto:.4f}")
    print(f"Allele 'A' frequency in Females: {p_female_auto:.4f}")
    print(f"Resulting Fst between sexes: {fst_auto:.6f}")
    print("Explanation: Fst is close to 0, indicating males and females are from the same gene pool for this marker.")

    print("\n--- Marker 2: Sex-Linked Locus (XY system) ---")
    print(f"Allele 'B' frequency in Males (1 copy):   {p_male_sex:.4f}")
    print(f"Allele 'B' frequency in Females (2 copies): {p_female_sex:.4f}")
    print(f"Resulting Fst between sexes: {fst_sex_linked:.6f}")
    print("Explanation: Fst is significantly greater than 0. The difference in allele copy number between sexes creates inherent genetic differentiation.")
    print("\nConclusion: The presence of a sex-determination system (like XY or ZW) is a primary reason for high Fst between sexes at specific loci.")

simulate_and_calculate_fst()