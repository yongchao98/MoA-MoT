def solve_rye_genetics():
    """
    This function calculates the expected percentage of drought-resistant offspring
    and the mass increase of the E3ub protein based on the provided data.
    """

    # Part 1: Calculate the percentage of resistant offspring

    # Given probabilities for pollination types
    prob_self_pollination = 0.05
    prob_cross_pollination = 0.95

    # In self-pollination (wt/insert x wt/insert), resistant offspring are
    # genotypes (wt/insert) and (insert/insert).
    # This corresponds to 3/4 of the progeny.
    percent_resistant_from_selfing = 0.75

    # In cross-pollination with the wild-type population (wt/insert x wt/wt),
    # resistant offspring are genotype (wt/insert).
    # This corresponds to 1/2 of the progeny.
    percent_resistant_from_crossing = 0.50

    # Calculate the total expected percentage of resistant offspring
    total_resistant_percentage = (prob_self_pollination * percent_resistant_from_selfing) + \
                                 (prob_cross_pollination * percent_resistant_from_crossing)

    # Convert to percentage format
    total_resistant_percentage *= 100

    print("Step 1: Calculating the theoretical percentage of resistant offspring.")
    print(f"The parent plant is heterozygous (wt/insert). Resistance is conferred by the 'insert' allele.")
    print("Scenario A: Self-pollination (occurs 5% of the time)")
    print(f"   - Offspring from (wt/insert) x (wt/insert) are 75% resistant.")
    print("Scenario B: Cross-pollination with wild-type population (occurs 95% of the time)")
    print(f"   - Offspring from (wt/insert) x (wt/wt) are 50% resistant.")
    print("\nCalculating the total probability:")
    print(f"({prob_self_pollination} * {percent_resistant_from_selfing*100}%) + ({prob_cross_pollination} * {percent_resistant_from_crossing*100}%) = {total_resistant_percentage:.2f}%")
    print("-" * 20)

    # Part 2: Calculate the mass increase of the protein
    insertion_length_nt = 105
    nucleotides_per_amino_acid = 3
    avg_mass_per_amino_acid_da = 110 # Daltons
    daltons_per_kilodalton = 1000

    # Calculate number of amino acids added
    num_added_amino_acids = insertion_length_nt / nucleotides_per_amino_acid

    # Calculate mass increase in Daltons
    mass_increase_da = num_added_amino_acids * avg_mass_per_amino_acid_da

    # Convert mass to kiloDaltons (kDa)
    mass_increase_kda = mass_increase_da / daltons_per_kilodalton

    print("Step 2: Calculating the protein mass increase from the nucleotide insertion.")
    print(f"The insertion has a length of {insertion_length_nt} nucleotides.")
    print("Each amino acid is coded by 3 nucleotides.")
    print(f"Number of added amino acids = {insertion_length_nt} / {nucleotides_per_amino_acid} = {int(num_added_amino_acids)}")
    print(f"Mass increase = {int(num_added_amino_acids)} amino acids * {avg_mass_per_amino_acid_da} Da/amino acid = {mass_increase_da} Da")
    print(f"Converting to kiloDaltons (kDa): {mass_increase_da} Da / {daltons_per_kilodalton} = {mass_increase_kda:.2f} kDa")
    print(f"This value is approximately 4.0 kDa.")
    print("-" * 20)
    
    print("Summary of conclusions based on experimental data:")
    print("1. Protein Activity: Only E3ub-wt is an active ubiquitin ligase (degrades Par22).")
    print("2. Protein Interaction: Only E3ub-wt interacts with Par22 (forms a complex). Par22 cannot interact with E3ub-insert105.")
    print("3. Offspring Resistance: 51.25% of offspring are expected to be resistant.")
    print("4. Mass Increase: The insertion adds approximately 4.0 kDa to the protein mass.")
    print("\nThese conclusions match option J.")

solve_rye_genetics()