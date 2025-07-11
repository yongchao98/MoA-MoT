import math

def solve_rye_genetics_problem():
    """
    This function analyzes the provided biological data to determine the correct statement among the choices.
    It calculates offspring genetics, interprets biochemical experiments, and calculates protein mass change.
    """

    # --- Step 1: Calculate the percentage of drought-resistant offspring ---

    self_pollination_rate = 0.05
    cross_pollination_rate = 0.95

    # Genotype of resistant parent: Ee (E=wild type, e=insertion/resistance allele)
    # Genotype of wild population for cross-pollination: EE

    # Offspring proportions from self-pollination (Ee x Ee) -> 0.25 EE, 0.50 Ee, 0.25 ee
    p_Ee_self = 0.50
    p_ee_self = 0.25

    # Offspring proportions from cross-pollination (Ee x EE) -> 0.50 EE, 0.50 Ee, 0.00 ee
    p_Ee_cross = 0.50
    p_ee_cross = 0.0

    # Calculate total probability for resistant genotypes (Ee and ee)
    # Total P(Ee) = P(Ee from self) + P(Ee from cross)
    total_p_Ee = (p_Ee_self * self_pollination_rate) + (p_Ee_cross * cross_pollination_rate)
    # Total P(ee) = P(ee from self) + P(ee from cross)
    total_p_ee = (p_ee_self * self_pollination_rate) + (p_ee_cross * cross_pollination_rate)

    # Total resistant offspring are those with at least one 'e' allele (Ee or ee)
    total_p_resistant = total_p_Ee + total_p_ee

    print("--- Step 1: Genetic Calculation ---")
    print("The final equation for the proportion of resistant offspring is:")
    print(f"P(Resistant) = (P(Ee_self) * {self_pollination_rate}) + (P(Ee_cross) * {cross_pollination_rate}) + (P(ee_self) * {self_pollination_rate}) + (P(ee_cross) * {cross_pollination_rate})")
    print(f"P(Resistant) = ({p_Ee_self} * {self_pollination_rate}) + ({p_Ee_cross} * {cross_pollination_rate}) + ({p_ee_self} * {self_pollination_rate}) + ({p_ee_cross} * {cross_pollination_rate})")
    print(f"P(Resistant) = {total_p_Ee:.4f} + {total_p_ee:.4f} = {total_p_resistant:.4f}")
    print(f"Theoretically, {total_p_resistant * 100:.2f}% of the offspring should be drought-resistant.\n")

    # --- Step 2: Interpretation of Biochemical Data ---

    print("--- Step 2: Biochemical Analysis ---")
    print("Co-expression results:")
    print("  - E3ub-wt reduces Par22 levels, indicating it is an active E3 ubiquitin ligase that targets Par22 for degradation.")
    print("  - E3ub-insert105 does not reduce Par22 levels, indicating the insertion inactivates its ligase function.")
    print("Conclusion 1: Only E3ub-wt is an active ubiquitin ligase.\n")

    print("Native mass spectrometry results:")
    print("  - Par22 (50 kDa) + E3ub-wt (60 kDa) form a 110 kDa complex, showing they interact.")
    print("  - Par22 + E3ub-insert105 show separate peaks, indicating they do NOT interact.")
    print("Conclusion 2: Par22 cannot interact with E3ub-insert105.\n")

    # --- Step 3: Calculation of Mass Increase ---

    insertion_length_nt = 105
    # Each amino acid is coded by 3 nucleotides
    added_amino_acids = insertion_length_nt / 3
    # Average molecular weight of an amino acid is ~110 Daltons
    avg_mass_aa_da = 110
    mass_increase_da = added_amino_acids * avg_mass_aa_da
    mass_increase_kda = mass_increase_da / 1000

    print("--- Step 3: Protein Mass Increase Calculation ---")
    print("The insertion adds a new sequence of 105 nucleotides.")
    print(f"Number of added amino acids = {insertion_length_nt} / 3 = {added_amino_acids:.0f}")
    print(f"Approximate mass increase = {added_amino_acids:.0f} amino acids * {avg_mass_aa_da} Da/aa = {mass_increase_da:.0f} Da")
    print(f"This is equal to {mass_increase_kda:.1f} kDa, which is approximately 4.0 kDa.\n")

    print("--- Step 4: Final Conclusion ---")
    print("Combining the findings:")
    print(f"1. Resistance in offspring: {total_p_resistant * 100:.2f}%")
    print("2. Ligase Activity: Only E3ub-wt is active.")
    print("3. Protein Interaction: Only E3ub-wt interacts with Par22.")
    print(f"4. Mass Increase: Approximately 4.0 kDa.")
    print("\nOption J is the only statement that correctly integrates all these points.")

# Run the analysis
solve_rye_genetics_problem()