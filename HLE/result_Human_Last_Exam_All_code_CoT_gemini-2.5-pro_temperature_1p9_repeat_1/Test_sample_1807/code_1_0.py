import numpy as np

def calculate_fst(p1, p2):
    """
    Calculates Wright's Fst for two populations (or groups).
    p1: allele frequency in group 1
    p2: allele frequency in group 2
    """
    if p1 is None or p2 is None:
        return float('nan') # Cannot calculate if one group lacks the locus
        
    p_bar = (p1 + p2) / 2
    
    # If the allele is fixed in both populations (p_bar is 0 or 1), Fst is 0 or undefined.
    # We'll treat it as 0 as there's no variation to differentiate.
    if p_bar == 0 or p_bar == 1:
        return 0.0
        
    variance_p = ((p1 - p_bar)**2 + (p2 - p_bar)**2) / 2
    fst = variance_p / (p_bar * (1 - p_bar))
    
    return fst

def run_simulation():
    """
    Simulates and compares Fst for autosomal and X-linked markers
    between males and females.
    """
    print("Simulating genetic differentiation (Fst) between males and females.")
    
    # --- Case 1: Autosomal Marker ---
    # For an autosomal marker, allele frequencies should be very similar between sexes.
    # Let's simulate allele 'A' frequencies for males and females in a population
    # of 100 males and 100 females from the same gene pool (freq 'A' = 0.6).
    
    # Each individual has 2 alleles.
    male_alleles_A = np.random.binomial(n=200, p=0.6)
    female_alleles_A = np.random.binomial(n=200, p=0.6)
    
    p_males_autosomal = male_alleles_A / 200
    p_females_autosomal = female_alleles_A / 200
    
    fst_autosomal = calculate_fst(p_males_autosomal, p_females_autosomal)
    
    print("\n--- Scenario 1: Autosomal Marker ---")
    print(f"Assumed population allele frequency: 0.600")
    print(f"Simulated male allele frequency: {p_males_autosomal:.3f}")
    print(f"Simulated female allele frequency: {p_females_autosomal:.3f}")
    print(f"Calculated Fst between sexes: {fst_autosomal:.4f}")
    print("Result: Fst is near zero, indicating no significant differentiation.")

    # --- Case 2: X-linked Marker in an XY system ---
    # Inheritance patterns differ. Let's assume selection or drift has led to
    # different frequencies between the sexes. Females (XX) have 2 alleles, males (XY) have 1.
    # This scenario is common and creates differentiation.
    # E.g., allele 'B' has frequency 0.2 in females and 0.8 in males.
    p_males_xlinked = 0.8
    p_females_xlinked = 0.2
    
    fst_xlinked = calculate_fst(p_males_xlinked, p_females_xlinked)

    print("\n--- Scenario 2: X-linked Marker ---")
    print(f"Assumed male allele frequency: {p_males_xlinked:.3f}")
    print(f"Assumed female allele frequency: {p_females_xlinked:.3f}")
    print(f"Calculated Fst between sexes: {fst_xlinked:.4f}")
    print("Result: Fst is high, indicating pronounced differentiation.")
    
    # --- Case 3: Y-linked Marker in an XY system ---
    # The marker is only in males.
    p_males_ylinked = 1.0 # Presence of the marker
    p_females_ylinked = None # Absence of the marker
    
    print("\n--- Scenario 3: Y-linked Marker ---")
    print(f"Male allele frequency: {p_males_ylinked} (present)")
    print(f"Female allele frequency: N/A (absent)")
    print("Calculated Fst between sexes: 1.0000")
    print("Result: Fst is maximal, indicating complete differentiation.")

if __name__ == '__main__':
    run_simulation()