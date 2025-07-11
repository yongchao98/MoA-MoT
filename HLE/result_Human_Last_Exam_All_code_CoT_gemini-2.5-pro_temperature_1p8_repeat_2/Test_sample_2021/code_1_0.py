def solve_biology_problem():
    """
    This function calculates the genetic and molecular properties described in the problem
    and identifies the correct answer choice.
    """

    # 1. Calculate the percentage of resistant offspring
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # In a self-cross of a heterozygous plant (Ee x Ee), offspring are 1/4 EE, 2/4 Ee, 1/4 ee.
    # The 'e' allele confers resistance, so Ee and ee are resistant.
    resistance_from_selfing = 0.75  # 3/4

    # In a cross with the wild-type population (Ee x EE), offspring are 1/2 Ee, 1/2 EE.
    # Only Ee is resistant.
    resistance_from_crossing = 0.50  # 1/2

    # Calculate the weighted average for the total resistant offspring
    total_resistance_decimal = (self_pollination_rate * resistance_from_selfing) + (cross_pollination_rate * resistance_from_crossing)
    total_resistance_percentage = total_resistance_decimal * 100

    print("--- Offspring Resistance Calculation ---")
    print(f"The parent plant is heterozygous for the resistance allele.")
    print(f"Self-pollination occurs at a rate of {self_pollination_rate * 100}%.")
    print(f"A self-cross results in {resistance_from_selfing * 100}% resistant offspring.")
    print(f"Cross-pollination occurs at a rate of {cross_pollination_rate * 100}%.")
    print(f"A cross with the wild-type population results in {resistance_from_crossing * 100}% resistant offspring.")
    print("\nThe final equation for the total percentage of resistant offspring is:")
    print(f"({self_pollination_rate} * {resistance_from_selfing}) + ({cross_pollination_rate} * {resistance_from_crossing}) = {total_resistance_decimal:.4f}")
    print(f"Theoretically, {total_resistance_percentage:.2f}% of the offspring should be drought-resistant.")
    print("\n" + "="*40 + "\n")

    # 2. Analyze protein function from experimental data
    print("--- Functional Analysis ---")
    print("1. E3 Ligase Activity:")
    print("   - Par22 + E3ub-wt: Par22 levels decrease (700 -> 200), so E3ub-wt IS an active E3 ubiquitin ligase.")
    print("   - Par22 + E3ub-insert105: Par22 levels do not decrease (700 -> 3000), so E3ub-insert105 is NOT active.")
    print("\n2. Protein Interaction:")
    print("   - Par22 (50kDa) + E3ub-wt (60kDa): A complex of 110kDa is formed. They DO interact.")
    print("   - Par22 (50kDa) + E3ub-insert105: Separate peaks observed. They DO NOT interact.")
    print("\n" + "="*40 + "\n")

    # 3. Calculate the mass increase of the protein
    nucleotide_insertion_length = 105
    nucleotides_per_amino_acid = 3
    avg_mass_per_amino_acid_da = 110

    # Calculate number of amino acids added
    amino_acids_added = nucleotide_insertion_length / nucleotides_per_amino_acid
    # Calculate mass increase in Daltons
    mass_increase_da = amino_acids_added * avg_mass_per_amino_acid_da
    # Convert to kiloDaltons (kDa)
    mass_increase_kda = mass_increase_da / 1000

    print("--- Mass Increase Calculation ---")
    print(f"Insertion length: {nucleotide_insertion_length} nucleotides")
    print(f"Amino acids added: {nucleotide_insertion_length} / {nucleotides_per_amino_acid} = {int(amino_acids_added)}")
    print(f"Estimated mass increase: {int(amino_acids_added)} aa * {avg_mass_per_amino_acid_da} Da/aa = {int(mass_increase_da)} Da")
    print(f"This is approximately {mass_increase_kda:.2f} kDa (~4.0 kDa).")
    print("\n" + "="*40 + "\n")
    
    # 4. Final Conclusion
    print("--- Conclusion ---")
    print("Combining the findings:")
    print(f"- Offspring resistance: {total_resistance_percentage:.2f}%")
    print("- E3 Ligase Activity: Only E3ub-wt is active.")
    print("- Interaction: Only E3ub-wt interacts with Par22 (or Par22 cannot interact with E3ub-insert105).")
    print(f"- Mass Increase: ~{mass_increase_kda:.1f} kDa")
    print("\nThese results match answer choice J.")


solve_biology_problem()