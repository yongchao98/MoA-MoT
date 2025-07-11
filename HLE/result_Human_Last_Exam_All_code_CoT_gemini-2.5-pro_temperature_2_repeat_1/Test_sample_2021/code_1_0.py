def solve_rye_genetics_problem():
    """
    Solves the provided biological problem by calculating genetic probabilities
    and interpreting experimental data.
    """

    # Part 1: Calculate the percentage of resistant offspring
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Parent is Ee (resistant), general population is EE (susceptible)
    # Self-pollination (Ee x Ee) -> offspring are 1/4 EE, 1/2 Ee, 1/4 ee.
    # Resistant offspring (Ee, ee) are 3/4 or 0.75
    resistance_from_selfing = 0.75

    # Cross-pollination (Ee x EE) -> offspring are 1/2 EE, 1/2 Ee.
    # Resistant offspring (Ee) are 1/2 or 0.50
    resistance_from_crossing = 0.50

    total_resistant_offspring_proportion = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
    total_resistant_offspring_percentage = total_resistant_offspring_proportion * 100

    print("Step 1: Inheritance Calculation")
    print(f"Contribution from self-pollination: {self_pollination_rate} * {resistance_from_selfing} = {self_pollination_rate * resistance_from_selfing}")
    print(f"Contribution from cross-pollination: {cross_pollination_rate} * {resistance_from_crossing} = {cross_pollination_rate * resistance_from_crossing}")
    print(f"Total proportion of resistant offspring = {total_resistant_offspring_proportion}")
    print(f"Therefore, theoretically, {total_resistant_offspring_percentage:.2f}% of the offspring should be drought-resistant.")
    print("-" * 20)

    # Part 2: Interpret E3 ligase activity
    # Control (Par22 alone): 700 units
    # Par22 + E3ub-wt: 200 units. Level dropped, indicating Par22 was degraded.
    # Par22 + E3ub-insert105: 3000 units. Level increased, indicating no degradation.
    print("Step 2: E3 Ligase Activity Analysis")
    print("The level of Par22 decreases with E3ub-wt, indicating it's an active E3 ubiquitin ligase that degrades Par22.")
    print("The level of Par22 is not decreased with E3ub-insert105, indicating it is not an active ligase.")
    print("Conclusion: Only E3ub-wt is an active ubiquitin ligase.")
    print("-" * 20)

    # Part 3: Interpret Protein Interaction
    # Par22 (50kDa) + E3ub-wt (60kDa) -> 110kDa complex. This shows interaction.
    # Par22 (50kDa) + E3ub-insert105 (~69kDa) -> separate peaks at 50kDa and 69kDa. This shows no interaction.
    print("Step 3: Protein Interaction Analysis")
    print("Native mass spectrometry shows a complex of 110 kDa for Par22 and E3ub-wt (50 kDa + 60 kDa), indicating they interact.")
    print("Separate peaks for Par22 and E3ub-insert105 show they do not interact.")
    print("Conclusion: Par22 cannot interact with E3ub-insert105.")
    print("-" * 20)

    # Part 4: Calculate the mass increase of the protein
    nucleotide_insertion_length = 105
    # 3 nucleotides code for 1 amino acid
    amino_acid_insertion_length = nucleotide_insertion_length / 3
    # Average molecular weight of an amino acid is ~110 Daltons
    avg_amino_acid_mass_da = 110
    mass_increase_da = amino_acid_insertion_length * avg_amino_acid_mass_da
    # Convert Daltons to kiloDaltons
    mass_increase_kda = mass_increase_da / 1000

    print("Step 4: Mass Increase Calculation")
    print(f"An insertion of {nucleotide_insertion_length} nucleotides codes for {int(amino_acid_insertion_length)} amino acids.")
    print(f"The estimated mass increase is {int(amino_acid_insertion_length)} amino acids * ~{avg_amino_acid_mass_da} Da/amino acid = {mass_increase_da} Da.")
    print(f"This is approximately {mass_increase_kda:.1f} kDa.")
    print("-" * 20)

    # Final Answer Assembly
    # Combining the conclusions:
    # 1. 51.25% resistant offspring.
    # 2. Only E3ub-wt is active.
    # 3. E3ub-insert105 does not interact with Par22.
    # 4. Mass increase is ~4.0 kDa.
    # This matches choice J.
    print("Final Conclusion:")
    final_statement = (
        f"Theoretically, {total_resistant_offspring_percentage:.2f}% of the offspring should be drought-resistant. "
        "Only E3ub-wt is an active ubiquitin ligase. "
        "Par22 cannot interact with E3ub-insert105. "
        f"The insertion increases the mass of E3ub by approximately {mass_increase_kda:.1f}kDA."
    )
    print(final_statement)


solve_rye_genetics_problem()
<<<J>>>