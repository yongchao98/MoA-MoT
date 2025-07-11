def solve_rye_genetics():
    """
    This function calculates the outcomes based on the provided biological data
    and identifies the correct conclusion.
    """

    # Part 1: Calculate the percentage of resistant offspring
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # In self-pollination (Wt/Ins x Wt/Ins), offspring are:
    # 1/4 Wt/Wt (not resistant), 2/4 Wt/Ins (resistant), 1/4 Ins/Ins (resistant)
    # Probability of resistant offspring from selfing is 3/4 = 0.75
    prob_resistant_from_selfing = 0.75

    # In cross-pollination (Wt/Ins x Wt/Wt), offspring are:
    # 1/2 Wt/Wt (not resistant), 1/2 Wt/Ins (resistant)
    # Probability of resistant offspring from crossing is 1/2 = 0.50
    prob_resistant_from_crossing = 0.50

    # Total probability of resistant offspring
    total_resistant_prob = (self_pollination_rate * prob_resistant_from_selfing) + \
                           (cross_pollination_rate * prob_resistant_from_crossing)
    
    total_resistant_percentage = total_resistant_prob * 100

    print("--- Offspring Resistance Calculation ---")
    print(f"The final equation for the percentage of resistant offspring is:")
    print(f"({self_pollination_rate} * {prob_resistant_from_selfing}) + ({cross_pollination_rate} * {prob_resistant_from_crossing}) = {total_resistant_prob}")
    print(f"Calculated percentage of resistant offspring: {total_resistant_percentage:.2f}%\n")

    # Part 2: Calculate the protein mass increase
    insertion_length_bp = 105
    bp_per_amino_acid = 3
    avg_mass_per_amino_acid_da = 110

    num_amino_acids = insertion_length_bp / bp_per_amino_acid
    mass_increase_da = num_amino_acids * avg_mass_per_amino_acid_da
    mass_increase_kda = mass_increase_da / 1000

    print("--- Protein Mass Increase Calculation ---")
    print(f"The insertion of {insertion_length_bp} bp codes for {int(num_amino_acids)} amino acids.")
    print(f"Estimated mass increase: {int(num_amino_acids)} aa * {avg_mass_per_amino_acid_da} Da/aa = {int(mass_increase_da)} Da")
    print(f"This is approximately {mass_increase_kda:.1f} kDa.\n")

    # Part 3: Summarize experimental findings and state the final conclusion
    print("--- Summary of Findings ---")
    print("1. Ligase Activity: The co-expression data shows Par22 levels drop with E3ub-wt (700 -> 200) but not with E3ub-insert105 (700 -> 3000). Therefore, only E3ub-wt is an active ubiquitin ligase.")
    print("2. Protein Interaction: The mass spectrometry data shows a complex (110kDa) for E3ub-wt and Par22, but free Par22 (50kDa) for E3ub-insert105. Therefore, Par22 cannot interact with E3ub-insert105.")
    print("\n--- Final Conclusion ---")
    print("Based on all calculations and experimental evidence, the correct statement is:")
    print("Theoretically, 51.25% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately 4.0kDA.")

solve_rye_genetics()
<<<J>>>