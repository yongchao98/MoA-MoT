def solve_biology_problem():
    """
    This function calculates and explains the solution to the given biology problem.
    """

    # Part 1: Calculating the percentage of resistant offspring
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Proportions from self-pollination (Rr x Rr)
    prop_RR_from_self = 0.25
    prop_Rr_from_self = 0.50

    # Proportions from cross-pollination (Rr x rr)
    prop_Rr_from_cross = 0.50

    # Total proportion of RR genotype
    # P(RR) = P(RR from self) * P(self)
    total_prop_RR = prop_RR_from_self * self_pollination_rate

    # Total proportion of Rr genotype
    # P(Rr) = (P(Rr from self) * P(self)) + (P(Rr from cross) * P(cross))
    total_prop_Rr = (prop_Rr_from_self * self_pollination_rate) + (prop_Rr_from_cross * cross_pollination_rate)

    # The resistant genotypes are RR and Rr (dominant inheritance)
    total_resistant_prop = total_prop_RR + total_prop_Rr
    total_resistant_percentage = total_resistant_prop * 100

    print("--- Calculation of Resistant Offspring Percentage ---")
    print("The final percentage is calculated by summing the proportions of resistant genotypes (RR and Rr) from both self- and cross-pollination scenarios.")
    print(f"Resistant Proportion = P(RR) + P(Rr)")
    print(f"Resistant Proportion = ({prop_RR_from_self} * {self_pollination_rate}) + (({prop_Rr_from_self} * {self_pollination_rate}) + ({prop_Rr_from_cross} * {cross_pollination_rate}))")
    print(f"Resistant Proportion = {total_prop_RR} + {total_prop_Rr} = {total_resistant_prop}")
    print(f"Final percentage of resistant offspring = {total_resistant_prop} * 100 = {total_resistant_percentage:.2f}%\n")


    # Part 2: Calculating the mass increase from the insertion
    nucleotide_insertion_length = 105
    nucleotides_per_codon = 3
    avg_mass_amino_acid_Da = 110

    num_amino_acids = nucleotide_insertion_length / nucleotides_per_codon
    mass_increase_Da = num_amino_acids * avg_mass_amino_acid_Da
    mass_increase_kDa = mass_increase_Da / 1000

    print("--- Calculation of Protein Mass Increase ---")
    print("The mass increase is calculated based on the number of amino acids added by the nucleotide insertion.")
    print(f"Mass Increase (kDa) = ({nucleotide_insertion_length} nucleotides / {nucleotides_per_codon} nt_per_aa) * {avg_mass_amino_acid_Da} Da_per_aa / 1000")
    print(f"Mass Increase = {mass_increase_kDa:.2f} kDa")
    print("This value is approximately 4.0 kDa.\n")

    # Part 3: Stating the final conclusion
    print("--- Summary of Findings ---")
    print(f"1. Offspring Resistance: Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant.")
    print("2. Ligase Activity: Co-expression data shows only E3ub-wt is an active ubiquitin ligase that degrades Par22.")
    print("3. Protein Interaction: Mass spectrometry shows that Par22 cannot interact with E3ub-insert105.")
    print(f"4. Mass Increase: The {nucleotide_insertion_length} nucleotide insertion increases the mass of E3ub by approximately 4.0 kDa.")

solve_biology_problem()
<<<J>>>