import random

def calculate_fst(p1, p2):
    """
    Calculates Fst and intermediate values for two subpopulations.
    p1: Allele frequency in subpopulation 1 (e.g., males)
    p2: Allele frequency in subpopulation 2 (e.g., females)
    Returns: (fst, p_bar, var_p)
    """
    p_bar = (p1 + p2) / 2
    
    # Handle cases where the allele is fixed across both populations
    if p_bar == 0 or p_bar == 1:
        # A special case is p1=0, p2=1 (or vice versa), where Fst is 1.
        if (p1 == 0 and p2 == 1) or (p1 == 1 and p2 == 0):
             var_p = ((p1 - p_bar)**2 + (p2 - p_bar)**2) / 2
             return 1.0, p_bar, var_p
        return 0.0, p_bar, 0.0

    var_p = ((p1 - p_bar)**2 + (p2 - p_bar)**2) / 2
    denominator = p_bar * (1 - p_bar)
    fst = var_p / denominator if denominator > 0 else 0.0
    
    return fst, p_bar, var_p

def run_simulation(pop_size=2000, p_autosomal_allele=0.6):
    """
    Simulates a population and calculates Fst between males and females
    for an autosomal and a Y-linked marker, showing the calculation steps.
    """
    num_males = pop_size // 2
    num_females = pop_size - num_males

    # --- 1. Autosomal Marker Simulation ---
    # Frequencies will differ slightly due to random sampling.
    male_autosomal_A_count = sum(1 for _ in range(2 * num_males) if random.random() < p_autosomal_allele)
    female_autosomal_A_count = sum(1 for _ in range(2 * num_females) if random.random() < p_autosomal_allele)
    
    p_male_auto = male_autosomal_A_count / (2 * num_males)
    p_female_auto = female_autosomal_A_count / (2 * num_females)

    # --- 2. Y-linked Marker Simulation ---
    # Frequency is 1 in males (all have a Y) and 0 in females (no Y).
    p_male_y = 1.0
    p_female_y = 0.0

    print("--- Fst Calculation: Males vs. Females ---")
    print(f"Formula: Fst = Variance(p) / (p_bar * (1 - p_bar))")

    # --- Calculation and Output for Autosomal Marker ---
    fst_auto, p_bar_auto, var_p_auto = calculate_fst(p_male_auto, p_female_auto)
    denominator_auto = p_bar_auto * (1 - p_bar_auto)

    print("\n1. Autosomal Marker:")
    print(f"   p(males)   = {p_male_auto:.4f}")
    print(f"   p(females) = {p_female_auto:.4f}")
    print(f"   p_bar = ({p_male_auto:.4f} + {p_female_auto:.4f}) / 2 = {p_bar_auto:.4f}")
    print(f"   Variance(p) = (({p_male_auto:.4f} - {p_bar_auto:.4f})^2 + ({p_female_auto:.4f} - {p_bar_auto:.4f})^2) / 2 = {var_p_auto:.8f}")
    print(f"   Fst = {var_p_auto:.8f} / ({p_bar_auto:.4f} * (1 - {p_bar_auto:.4f})) = {fst_auto:.6f}")
    print("   Result: Fst is near 0, indicating negligible differentiation.")

    # --- Calculation and Output for Y-linked Marker ---
    fst_y, p_bar_y, var_p_y = calculate_fst(p_male_y, p_female_y)
    denominator_y = p_bar_y * (1 - p_bar_y)
    
    print("\n2. Y-linked Marker:")
    print(f"   p(males)   = {p_male_y:.4f}")
    print(f"   p(females) = {p_female_y:.4f}")
    print(f"   p_bar = ({p_male_y:.4f} + {p_female_y:.4f}) / 2 = {p_bar_y:.4f}")
    print(f"   Variance(p) = (({p_male_y:.4f} - {p_bar_y:.4f})^2 + ({p_female_y:.4f} - {p_bar_y:.4f})^2) / 2 = {var_p_y:.8f}")
    print(f"   Fst = {var_p_y:.8f} / ({p_bar_y:.4f} * (1 - {p_bar_y:.4f})) = {fst_y:.6f}")
    print("   Result: Fst is 1, indicating complete differentiation.")

# Run the simulation
run_simulation()