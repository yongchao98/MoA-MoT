import math

def solve_biology_problem():
    """
    This function calculates the quantitative parts of the problem to find the correct answer.
    It calculates:
    1. The expected frequency of resistant offspring.
    2. The increase in protein mass due to the nucleotide insertion.
    """

    # --- Part 1: Calculate the frequency of resistant offspring ---

    # Pollination rates
    self_pollination_rate = 0.05  # 5%
    cross_pollination_rate = 0.95 # 95%

    # Probabilities of resistant offspring for each pollination type
    # Self-pollination (Rr x Rr -> RR, Rr, rR) -> 3/4 are resistant
    prob_resistant_from_selfing = 0.75
    # Cross-pollination (Rr x rr -> Rr, Rr) -> 1/2 are resistant
    prob_resistant_from_crossing = 0.50

    # Total probability calculation
    total_prob_resistant = (self_pollination_rate * prob_resistant_from_selfing) + (cross_pollination_rate * prob_resistant_from_crossing)
    total_prob_resistant_percent = total_prob_resistant * 100

    print("--- Offspring Resistance Frequency Calculation ---")
    print("The parent plant is heterozygous (Rr).")
    print("Probability of resistant offspring from self-pollination (Rr x Rr) = 3/4 = 0.75")
    print("Probability of resistant offspring from cross-pollination (Rr x rr) = 1/2 = 0.50")
    print("\nFinal Equation for Total Probability:")
    print(f"P(Resistant) = (Self-Pollination Rate * P(Resistant | Self)) + (Cross-Pollination Rate * P(Resistant | Cross))")
    print(f"P(Resistant) = ({self_pollination_rate} * {prob_resistant_from_selfing}) + ({cross_pollination_rate} * {prob_resistant_from_crossing})")
    print(f"P(Resistant) = {self_pollination_rate * prob_resistant_from_selfing} + {cross_pollination_rate * prob_resistant_from_crossing}")
    print(f"P(Resistant) = {total_prob_resistant}")
    print(f"Theoretically, {total_prob_resistant_percent:.2f}% of the offspring should be drought-resistant.")


    # --- Part 2: Calculate the protein mass increase ---

    insertion_sequence = "gatgaaaaagtgtttaccaaagaactggatcagaccattgaacagctgaacgaatgcggccagctgagcgaaagccaggtgaaaagcctgatgtgcgaagcgacc"
    nucleotide_count = len(insertion_sequence)
    nucleotides_per_codon = 3
    avg_amino_acid_mass_da = 110  # Daltons

    # Calculation
    num_amino_acids = nucleotide_count / nucleotides_per_codon
    mass_increase_da = num_amino_acids * avg_amino_acid_mass_da
    mass_increase_kda = mass_increase_da / 1000

    print("\n--- Protein Mass Increase Calculation ---")
    print(f"The insertion has {nucleotide_count} nucleotides.")
    print("Each codon consists of 3 nucleotides.")
    print("The average mass of an amino acid is ~110 Daltons.")
    print("\nFinal Equation for Mass Increase:")
    print(f"Amino Acids Added = Nucleotide Count / Nucleotides per Codon")
    print(f"Amino Acids Added = {nucleotide_count} / {nucleotides_per_codon} = {int(num_amino_acids)}")
    print(f"\nMass Increase (Da) = Amino Acids Added * Average Mass per Amino Acid")
    print(f"Mass Increase (Da) = {int(num_amino_acids)} * {avg_amino_acid_mass_da} = {int(mass_increase_da)}")
    print(f"\nMass Increase (kDa) = {int(mass_increase_da)} / 1000 = {mass_increase_kda:.2f} kDa")
    print(f"The insertion increases the mass of the protein by approximately {round(mass_increase_kda)} kDa.")

solve_biology_problem()