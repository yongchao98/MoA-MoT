import random

def explain_sex_linked_fst():
    """
    Simulates genetic data for autosomal and Y-linked markers to demonstrate
    why sex-determining systems can cause high Fst between males and females.
    """
    num_males = 100
    num_females = 100
    print(f"Simulating a population with {num_males} males (XY) and {num_females} females (XX).\n")

    # --- Case 1: Autosomal Marker ---
    # Both sexes have two copies. Alleles 'A' and 'a'.
    # Frequencies should be similar, leading to low Fst.
    print("--- Scenario 1: Autosomal Marker ---")
    print("Markers on non-sex chromosomes are inherited similarly by males and females.")

    # Generate random allele frequencies for demonstration
    p_A = random.uniform(0.4, 0.6) # Allele 'A' frequency
    q_A = 1 - p_A               # Allele 'a' frequency

    # In a large, randomly mating population, frequencies between sexes should be nearly identical.
    # We simulate small sampling differences.
    male_alleles = random.choices(['A', 'a'], weights=[p_A, q_A], k=2 * num_males)
    female_alleles = random.choices(['A', 'a'], weights=[p_A, q_A], k=2 * num_females)

    p_male = male_alleles.count('A') / (2 * num_males)
    q_male = 1 - p_male
    p_female = female_alleles.count('A') / (2 * num_females)
    q_female = 1 - p_female
    
    total_alleles = male_alleles + female_alleles
    p_total = total_alleles.count('A') / len(total_alleles)
    q_total = 1 - p_total

    # Calculate Fst = (Ht - Hs) / Ht
    Ht = 2 * p_total * q_total
    Hs_male = 2 * p_male * q_male
    Hs_female = 2 * p_female * q_female
    Hs = (Hs_male + Hs_female) / 2
    
    fst_autosomal = (Ht - Hs) / Ht if Ht > 0 else 0

    print("Calculating Fst, which measures differentiation:")
    print(f"Ht (total heterozygosity) = 2 * p_total * q_total = 2 * {p_total:.4f} * {q_total:.4f} = {Ht:.4f}")
    print(f"Hs (avg. subpopulation heterozygosity) = (H_male + H_female) / 2 = ({Hs_male:.4f} + {Hs_female:.4f}) / 2 = {Hs:.4f}")
    print("\nFinal Fst Equation:")
    print(f"Fst = (Ht - Hs) / Ht = ({Ht:.4f} - {Hs:.4f}) / {Ht:.4f} = {fst_autosomal:.4f}")
    print("\nResult: For an autosomal marker, Fst is very close to 0, indicating no significant differentiation.\n")


    # --- Case 2: Y-linked Marker ---
    # Only males have a Y chromosome.
    # This leads to maximum Fst.
    print("--- Scenario 2: Y-linked Marker ---")
    print("Markers on the Y chromosome exist only in males.")
    
    male_y_alleles = ['Y_allele_1', 'Y_allele_2']
    female_y_alleles = ['None']

    print(f"Alleles present in Males: {male_y_alleles}")
    print(f"Alleles present in Females: {female_y_alleles}")
    print("\nResult: The allele pools for males and females are completely different (disjoint).")
    print("This represents the maximum possible genetic differentiation.")
    print("\nFinal Fst Equation:")
    print("Fst = 1.0")
    print("\nConclusion: A marker on a sex chromosome (like the Y) is a strong potential explanation for pronounced differentiation.")

explain_sex_linked_fst()