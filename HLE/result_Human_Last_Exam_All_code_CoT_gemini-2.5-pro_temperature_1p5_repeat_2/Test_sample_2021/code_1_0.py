import math

def solve_biology_problem():
    """
    Analyzes the provided biological data to determine the correct conclusions
    and calculates the required numerical values.
    """

    # Step 1: Define initial parameters from the problem description
    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95
    insertion_nt_length = 105
    avg_aa_mass_da = 110  # Average mass of an amino acid in Daltons

    # The original resistant plant is heterozygous (E3ub-wt/E3ub-insert105).
    # Resistance is associated with the 'insert' allele, which prevents Par22 degradation.
    # Therefore, genotypes wt/insert and insert/insert are resistant.

    # For self-pollination (wt/insert x wt/insert), offspring are:
    # 1/4 wt/wt (non-resistant), 2/4 wt/insert (resistant), 1/4 insert/insert (resistant)
    resistance_rate_selfing = 0.75

    # For cross-pollination with wild-type (wt/insert x wt/wt), offspring are:
    # 1/2 wt/wt (non-resistant), 1/2 wt/insert (resistant)
    resistance_rate_crossing = 0.50

    # Step 2: Analyze qualitative experimental data
    # Co-expression: E3ub-wt reduces Par22 levels, while E3ub-insert105 does not.
    # This implies only E3ub-wt is an active ubiquitin ligase targeting Par22.
    activity_conclusion = "Only E3ub-wt is an active ubiquitin ligase."

    # Mass Spectrometry: E3ub-wt (60kDa) + Par22 (50kDa) forms a 110kDa complex.
    # E3ub-insert105 + Par22 results in separate peaks (no complex).
    # This implies Par22 interacts with E3ub-wt but not with E3ub-insert105.
    interaction_conclusion = "Par22 cannot interact with E3ub-insert105."

    # Step 3: Calculate the mass increase of the protein
    # 3 nucleotides code for 1 amino acid
    added_amino_acids = insertion_nt_length / 3
    mass_increase_da = added_amino_acids * avg_aa_mass_da
    mass_increase_kda = mass_increase_da / 1000

    # Step 4: Calculate the total percentage of resistant offspring
    total_resistant_offspring_rate = (self_pollination_rate * resistance_rate_selfing) + (cross_pollination_rate * resistance_rate_crossing)
    percentage_resistant_offspring = total_resistant_offspring_rate * 100

    # Step 5: Print the comprehensive analysis and results
    print("Summary of Analysis:")
    print("-" * 30)

    print("1. Protein Activity & Interaction Analysis:")
    print(f"- Conclusion on activity: {activity_conclusion}")
    print(f"- Conclusion on interaction: {interaction_conclusion}")
    print("-" * 30)

    print("2. Protein Mass Increase Calculation:")
    print(f"The insertion of {insertion_nt_length} nucleotides adds {int(added_amino_acids)} amino acids.")
    print(f"The estimated mass increase is {int(mass_increase_da)} Da, which is approximately {mass_increase_kda:.1f} kDa.")
    print("-" * 30)

    print("3. Resistant Offspring Percentage Calculation:")
    print("The final equation for the total percentage of resistant offspring is:")
    print(f"   (Self-pollination Rate * Resistance from Selfing) + (Cross-pollination Rate * Resistance from Crossing)")
    # The final instruction requests to output each number in the final equation.
    print(f"   ({self_pollination_rate} * {resistance_rate_selfing}) + ({cross_pollination_rate} * {resistance_rate_crossing}) = {total_resistant_offspring_rate}")
    print(f"The theoretical percentage of drought-resistant offspring is {percentage_resistant_offspring:.2f}%.")
    print("-" * 30)

    print("\nFinal Combined Conclusion:")
    print("Theoretically, 51.25% of the offspring should be drought-resistant. Only E3ub-wt is an active ubiquitin ligase. Par22 cannot interact with E3ub-insert105. The insertion increases the mass of E3ub by approximately 4.0 kDa.")


if __name__ == "__main__":
    solve_biology_problem()