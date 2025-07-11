def solve_biology_problem():
    """
    This function calculates the mass increase of a protein from a nucleotide insertion
    and the theoretical percentage of resistant offspring based on pollination rates
    and genetic principles. It then prints a step-by-step analysis to determine the
    correct multiple-choice answer.
    """

    # --- Step 1: Calculate Mass Increase ---
    insertion_length_nt = 105
    nt_per_codon = 3
    # The average molecular weight of an amino acid is ~110 Daltons.
    avg_aa_mass_da = 110

    added_amino_acids = insertion_length_nt / nt_per_codon
    mass_increase_da = added_amino_acids * avg_aa_mass_da
    mass_increase_kda = mass_increase_da / 1000

    print("--- Analysis ---")
    print("\nStep 1: Calculate the mass increase of the E3ub protein.")
    print(f"The insertion adds {added_amino_acids} amino acids ({insertion_length_nt} nt / {nt_per_codon} nt/codon).")
    print(f"The calculated mass increase is {mass_increase_da} Da, or approximately {mass_increase_kda:.1f} kDa.")
    print("This confirms the mass increase is approximately 4.0 kDa.")

    # --- Step 2: Interpret Experimental Data ---
    print("\nStep 2: Interpret the protein activity and interaction data.")
    print("- Co-expression data shows only E3ub-wt degrades Par22, so only E3ub-wt is an active ubiquitin ligase.")
    print("- Mass spectrometry data shows only E3ub-wt forms a complex with Par22, meaning Par22 cannot interact with E3ub-insert105.")

    # --- Step 3: Calculate Offspring Resistance Percentage ---
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    # For a dominant trait from a heterozygous parent (Wt/Insert):
    # Self-pollination (Wt/Insert x Wt/Insert) yields 75% resistant offspring.
    resistance_from_selfing = 0.75
    # Cross-pollination (Wt/Insert x Wt/Wt) yields 50% resistant offspring.
    resistance_from_crossing = 0.50

    total_resistance = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
    total_resistance_percent = total_resistance * 100

    print("\nStep 3: Calculate the theoretical percentage of resistant offspring.")
    print("The resistance allele is dominant. The total percentage of resistant offspring is calculated based on pollination methods.")
    print("Final Equation: Total Resistance = (Self-pollination Rate * Resistance from Selfing) + (Cross-pollination Rate * Resistance from Crossing)")
    # Remember to output each number in the final equation!
    print(f"Calculation: Total Resistance = ({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing})")
    print(f"Result: Total Resistance = {total_resistance}")
    print(f"This equals {total_resistance_percent:.2f}% of offspring being resistant.")

    print("\n--- Conclusion ---")
    print("The correct answer choice must state that:")
    print(f"1. The theoretical resistance is {total_resistance_percent:.2f}%.")
    print("2. Only E3ub-wt is an active ubiquitin ligase.")
    print("3. Par22 cannot interact with E3ub-insert105.")
    print(f"4. The insertion increases E3ub mass by ~4.0 kDa.")
    print("Answer choice 'J' matches all these points.")

solve_biology_problem()
<<<J>>>