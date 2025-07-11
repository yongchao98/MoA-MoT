import math

def solve_biology_problem():
    """
    This function calculates and prints the step-by-step solution to the problem.
    """

    # --- Step 1: Calculate the mass increase from the insertion ---
    print("--- Step 1: Mass Increase Calculation ---")
    insertion_nucleotides = 105
    # Every 3 nucleotides code for 1 amino acid.
    amino_acid_increase = insertion_nucleotides / 3
    # The average molecular weight of an amino acid is ~110 Daltons (Da).
    avg_aa_mass_da = 110
    # Mass increase in Daltons
    mass_increase_da = amino_acid_increase * avg_aa_mass_da
    # Convert Daltons to kiloDaltons (kDa)
    mass_increase_kda = mass_increase_da / 1000
    
    print(f"The insertion is {insertion_nucleotides} nucleotides long.")
    print(f"This translates to {insertion_nucleotides} / 3 = {amino_acid_increase:.0f} amino acids.")
    print(f"The estimated mass increase is {amino_acid_increase:.0f} amino acids * {avg_aa_mass_da} Da/amino acid = {mass_increase_da:.0f} Da.")
    print(f"This is equal to {mass_increase_kda:.2f} kDa, which is approximately 4.0 kDa.\n")
    
    # --- Step 2: Interpret the experimental data ---
    print("--- Step 2: Interpretation of Experimental Data ---")
    print("Co-expression results show Par22 levels drop with E3ub-wt, indicating it's an active ubiquitin ligase that degrades Par22.")
    print("Par22 levels increase with E3ub-insert105, indicating it is NOT an active ubiquitin ligase.")
    print("Native Mass Spectrometry shows a 110kDa complex for E3ub-wt + Par22 (50kDa + 60kDa), indicating they interact.")
    print("For E3ub-insert105, separate peaks are observed (including a 50kDa peak for free Par22), indicating Par22 cannot interact with E3ub-insert105.\n")

    # --- Step 3: Calculate the percentage of resistant offspring ---
    print("--- Step 3: Offspring Resistance Calculation ---")
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    
    # In self-pollination (wt/insert x wt/insert), offspring are 1/4 wt/wt, 2/4 wt/insert, 1/4 insert/insert.
    # Resistance is dominant, so 3/4 of offspring are resistant.
    resistance_rate_selfing = 0.75
    
    # In cross-pollination (wt/insert x wt/wt), offspring are 1/2 wt/wt, 1/2 wt/insert.
    # 1/2 of offspring are resistant.
    resistance_rate_crossing = 0.50
    
    total_resistance_rate = (self_pollination_rate * resistance_rate_selfing) + (cross_pollination_rate * resistance_rate_crossing)
    total_resistance_percentage = total_resistance_rate * 100

    print("The final proportion of resistant offspring is the sum of contributions from self-pollination and cross-pollination.")
    print(f"Final Equation: ({self_pollination_rate} * {resistance_rate_selfing}) + ({cross_pollination_rate} * {resistance_rate_crossing}) = {total_resistance_rate}")
    print(f"Therefore, the theoretical percentage of resistant offspring is {total_resistance_rate} * 100 = {total_resistance_percentage:.2f}%.\n")
    
    # --- Step 4: Synthesize the final conclusion ---
    print("--- Step 4: Final Conclusion ---")
    final_sentence = (
        f"Theoretically, {total_resistance_percentage:.2f}% of the offspring should be drought-resistant. "
        "Only E3ub-wt is an active ubiquitin ligase. "
        "Par22 cannot interact with E3ub-insert105. "
        f"The insertion increases the mass of E3ub by approximately {round(mass_increase_kda, 1)}kDA."
    )
    print(final_sentence)

# Execute the function to solve the problem
solve_biology_problem()