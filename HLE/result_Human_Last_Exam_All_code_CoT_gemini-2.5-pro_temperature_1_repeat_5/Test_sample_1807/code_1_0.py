import numpy as np

def calculate_fst(p1, p2, n1, n2):
    """
    Calculates Fst between two subpopulations for a bi-allelic locus.
    Fst = (Ht - Hs) / Ht
    
    Args:
        p1 (float): Allele 'A' frequency in subpopulation 1 (e.g., males).
        p2 (float): Allele 'A' frequency in subpopulation 2 (e.g., females).
        n1 (int): Number of individuals in subpopulation 1.
        n2 (int): Number of individuals in subpopulation 2.
        
    Returns:
        float: The Fst value. Returns -1 if Ht is zero.
    """
    # Calculate average allele frequency in the total population
    p_total = (p1 * n1 + p2 * n2) / (n1 + n2)
    
    # Expected heterozygosity in the total population (Ht)
    # Ht = 2 * p * q
    ht = 2 * p_total * (1 - p_total)
    
    # Expected heterozygosity within each subpopulation
    h1 = 2 * p1 * (1 - p1)
    h2 = 2 * p2 * (1 - p2)
    
    # Average heterozygosity across subpopulations (Hs)
    hs = (h1 * n1 + h2 * n2) / (n1 + n2)
    
    if ht == 0:
        # Avoid division by zero; if Ht is 0, Fst is typically considered 0 or undefined.
        # In our Y-chromosome case, it signifies complete differentiation.
        if hs == 0 and p1 != p2: # Case of fixed but different alleles
             return 1.0
        return 0.0

    fst = (ht - hs) / ht
    return fst

def run_simulation():
    """
    Simulates genetic data and calculates Fst to explain differentiation between sexes.
    """
    print("### Analyzing Genetic Differentiation Between Males and Females ###\n")
    print("The question asks for a potential explanation for why *some* genetic markers show")
    print("pronounced differentiation (high Fst) between males and females in a population.")
    print("The most direct cause is the presence of sex chromosomes (XY or ZW systems).\n")
    print("Let's demonstrate this with a simulation of an XY system.\n")

    # --- Simulation Parameters ---
    num_males = 100
    num_females = 100

    # --- Case 1: Autosomal Marker (present in both sexes) ---
    # Allele frequencies are expected to be similar between males and females.
    # Let's assume allele 'A' has a frequency of 0.7
    p_autosomal_males = 0.7
    p_autosomal_females = 0.7
    
    fst_autosomal = calculate_fst(p_autosomal_males, p_autosomal_females, num_males, num_females)
    
    print("--- Case 1: Autosomal Marker ---")
    print(f"This marker is on a non-sex chromosome, inherited by both sexes.")
    print(f"Allele frequency in Males:  {p_autosomal_males}")
    print(f"Allele frequency in Females: {p_autosomal_females}")
    print(f"Calculated Fst between males and females: {fst_autosomal:.4f}")
    print("Result: As expected, Fst is 0, indicating no genetic differentiation.\n")
    
    # --- Case 2: Y-linked Marker (present only in males) ---
    # We can model this by saying males have allele 'Y1' (freq=1.0) and 
    # females have a 'null' allele (freq=0.0 for 'Y1').
    p_y_linked_males = 1.0  # All males have the marker
    p_y_linked_females = 0.0 # All females lack the marker
    
    fst_y_linked = calculate_fst(p_y_linked_males, p_y_linked_females, num_males, num_females)

    print("--- Case 2: Y-linked Marker ---")
    print("This marker is on the Y chromosome, inherited only by males.")
    print(f"Frequency of the marker allele in Males:   {p_y_linked_males}")
    print(f"Frequency of the marker allele in Females:  {p_y_linked_females}")
    print(f"Calculated Fst between males and females: {fst_y_linked:.4f}")
    print("Result: Fst is 1, indicating complete genetic differentiation.\n")
    
    # --- Conclusion ---
    print("### Conclusion ###")
    print("The simulation shows that markers on sex chromosomes (like the Y-linked marker)")
    print("will exhibit maximum genetic differentiation between sexes, while autosomal markers will not.")
    print("This perfectly matches the observation that 'some markers exhibit pronounced differentiation'.")
    print("Therefore, the underlying biological mechanism is the sex-determination system.")

run_simulation()

<<<B>>>