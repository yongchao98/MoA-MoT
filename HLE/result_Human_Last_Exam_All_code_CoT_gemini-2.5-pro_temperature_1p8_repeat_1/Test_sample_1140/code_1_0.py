def calculate_synthesis_yields():
    """
    Calculates and compares the theoretical yields for synthesizing a 100aa peptide
    using direct Solid-Phase Peptide Synthesis (SPPS) versus Native Chemical Ligation (NCL).
    """

    # --- Parameters ---
    peptide_length = 100
    # A single coupling efficiency of 99% is optimistic for SPPS, especially for long peptides.
    coupling_efficiency_per_step = 0.99
    # Ligation efficiency is typically high for the NCL reaction.
    ligation_efficiency = 0.90 # 90%

    print("--- Comparing Peptide Synthesis Strategies ---\n")
    print(f"Peptide Length: {peptide_length} amino acids")
    print(f"Assumed SPPS Coupling Efficiency per Step: {coupling_efficiency_per_step * 100}%\n")

    # --- 1. Direct SPPS Method ---
    print("Method 1: Direct Solid-Phase Peptide Synthesis (SPPS)")
    num_spps_couplings = peptide_length - 1
    # The final yield is the efficiency raised to the power of the number of steps.
    yield_spps = coupling_efficiency_per_step ** num_spps_couplings

    print(f"Number of coupling steps = {peptide_length} - 1 = {num_spps_couplings}")
    print(f"Final Yield Equation: {coupling_efficiency_per_step}^{num_spps_couplings}")
    print(f"Calculated Direct SPPS Yield: {yield_spps:.2%}\n")


    # --- 2. Native Chemical Ligation (NCL) Method ---
    print("Method 2: Native Chemical Ligation (NCL)")
    # We split the 100aa peptide into two 50aa fragments.
    fragment1_len = 50
    fragment2_len = 50

    # Calculate yield for Fragment 1
    num_frag1_couplings = fragment1_len - 1
    yield_frag1 = coupling_efficiency_per_step ** num_frag1_couplings
    print(f"Fragment 1: {fragment1_len} amino acids")
    print(f"  - Number of coupling steps = {fragment1_len} - 1 = {num_frag1_couplings}")
    print(f"  - Synthesis Yield Equation: {coupling_efficiency_per_step}^{num_frag1_couplings}")
    print(f"  - Calculated Fragment 1 Yield: {yield_frag1:.2%}\n")

    # The synthesis yield for Fragment 2 is the same.
    # The amount of peptide available for ligation is determined by the yield of fragment synthesis.
    yield_of_fragments = yield_frag1

    # Final NCL yield is the product of the fragment synthesis yield and the ligation yield.
    yield_ncl = yield_of_fragments * ligation_efficiency

    print(f"Ligation Step:")
    print(f"  - Ligation Reaction Efficiency: {ligation_efficiency * 100}%")
    print(f"Final Yield Equation: (Fragment Synthesis Yield) * (Ligation Efficiency)")
    print(f"Final Yield Equation with numbers: {yield_of_fragments:.4f} * {ligation_efficiency}")
    print(f"Calculated Final NCL Yield: {yield_ncl:.2%}\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    if yield_ncl > yield_spps:
        print(f"The NCL approach ({yield_ncl:.2%}) provides a significantly higher theoretical yield than direct SPPS ({yield_spps:.2%}).")
        print("Therefore, Native Chemical Ligation is the most helpful technique for this synthesis.")
    else:
        print("Based on these parameters, direct SPPS appears more efficient.")

if __name__ == "__main__":
    calculate_synthesis_yields()