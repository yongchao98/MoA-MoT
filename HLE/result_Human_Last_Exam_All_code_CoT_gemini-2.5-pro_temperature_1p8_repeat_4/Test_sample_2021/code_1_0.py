import math

def solve_rye_genetics_problem():
    """
    This function analyzes the provided biological data and calculates the theoretical
    percentage of resistant offspring, then prints a comprehensive explanation.
    """

    # Part 1: Molecular Biology Analysis

    # Calculate the mass increase from the insertion
    insertion_nucleotides = 105
    codons_per_amino_acid = 3
    avg_mass_per_amino_acid_Da = 110 # Approximate average mass in Daltons
    Da_per_kDa = 1000

    amino_acids_added = insertion_nucleotides / codons_per_amino_acid
    mass_increase_Da = amino_acids_added * avg_mass_per_amino_acid_Da
    mass_increase_kDa = mass_increase_Da / Da_per_kDa

    # Analyze protein activity and interaction based on the text
    # Densitometry: Par22 (control 700), with E3ub-wt (200), with E3ub-insert105 (3000)
    # Conclusion: E3ub-wt degrades Par22, so it's an active ligase. E3ub-insert105 does not.
    activity_conclusion = "Only E3ub-wt is an active ubiquitin ligase that degrades Par22."

    # Mass Spec: Par22 (50kDa) + E3ub-wt (60kDa) -> 110kDa complex. Par22 + E3ub-insert105 -> 50kDa peak (free Par22)
    # Conclusion: Par22 interacts with E3ub-wt, but not with E3ub-insert105.
    interaction_conclusion = "Par22 cannot interact with E3ub-insert105, but it can interact with E3ub-wt."


    # Part 2: Genetic Calculation

    # Define pollination rates
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Calculate the probability of resistant offspring for each pollination type.
    # The parent is heterozygous (WI) and resistant. The population is homozygous wild-type (WW).
    # Resistance is conferred by having at least one 'I' allele.
    # Case 1: Self-pollination (WI x WI -> 1/4 WW, 2/4 WI, 1/4 II)
    prob_resistant_self = 0.75  # (WI + II = 3/4)
    # Case 2: Cross-pollination (WI x WW -> 1/2 WI, 1/2 WW)
    prob_resistant_cross = 0.50 # (WI = 1/2)

    # Calculate the total probability using the law of total probability
    total_prob_resistant = (prob_resistant_self * self_pollination_rate) + (prob_resistant_cross * cross_pollination_rate)
    percentage_resistant = total_prob_resistant * 100

    # Part 3: Print the full analysis
    print("--- Molecular Biology Analysis ---")
    print(f"1. Mass Increase Calculation:")
    print(f"   - The 105 nucleotide insertion codes for {int(amino_acids_added)} amino acids.")
    print(f"   - The estimated mass increase is {int(mass_increase_Da)} Da, or approximately {mass_increase_kDa:.1f} kDa.")
    print(f"2. Activity and Interaction Summary:")
    print(f"   - Ligase Activity: {activity_conclusion}")
    print(f"   - Protein Interaction: {interaction_conclusion}")
    print("\n--- Genetic Offspring Calculation ---")
    print("The theoretical percentage of resistant offspring is calculated as a weighted average of self-pollination and cross-pollination outcomes.")
    print("\nFinal Equation:")
    print(f"({prob_resistant_self} * {self_pollination_rate}) + ({prob_resistant_cross} * {cross_pollination_rate}) = {total_prob_resistant:.4f}")
    print(f"\nResult: The theoretical percentage of offspring that are drought-resistant is {percentage_resistant:.2f}%.")
    print("\n--- Final Conclusion ---")
    print("The correct answer must state that 51.25% of offspring are resistant, only E3ub-wt is an active ligase, Par22 cannot interact with E3ub-insert105, and the insertion adds ~4.0 kDa. This matches choice J.")

# Execute the analysis
solve_rye_genetics_problem()