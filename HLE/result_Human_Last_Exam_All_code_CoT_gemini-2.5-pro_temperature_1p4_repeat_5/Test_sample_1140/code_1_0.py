import sys

def calculate_synthesis_yields():
    """
    Calculates and compares the theoretical yields for synthesizing a long peptide
    using Direct Solid-Phase Peptide Synthesis (SPPS) vs. Native Chemical Ligation (NCL).
    """
    # --- Parameters ---
    total_length = 100
    # A realistic, yet optimistic, coupling efficiency for each amino acid addition in SPPS.
    coupling_efficiency = 0.99
    # A typical yield for the NCL reaction step after fragments are purified.
    ligation_yield = 0.70

    print("--- Analysis of Synthesis Strategies for a 100-Residue Peptide ---")
    print(f"Peptide Length: {total_length} amino acids")
    print(f"Assumed SPPS Coupling Efficiency per Step: {coupling_efficiency * 100:.1f}%")
    print("-" * 60)

    # --- Strategy 1: Direct SPPS ---
    print("\nStrategy 1: Direct Solid-Phase Peptide Synthesis (SPPS)")
    num_couplings_direct = total_length - 1
    direct_spps_yield = coupling_efficiency ** num_couplings_direct

    print(f"This strategy requires {num_couplings_direct} consecutive coupling steps.")
    print("The final yield is calculated by the equation: (coupling_efficiency) ^ (number_of_couplings)")
    print(f"Equation with numbers: {coupling_efficiency}^{num_couplings_direct}")
    print(f"Calculated Theoretical Yield: {direct_spps_yield:.4f} (or {direct_spps_yield * 100:.2f}%)")
    print("This low yield results in a difficult-to-purify mixture.")
    print("-" * 60)

    # --- Strategy 2: Native Chemical Ligation (NCL) ---
    print("\nStrategy 2: Native Chemical Ligation (NCL)")
    # For NCL, we can split the 100aa peptide into two 50aa fragments.
    fragment_length = 50
    num_couplings_fragment = fragment_length - 1
    fragment_synthesis_yield = coupling_efficiency ** num_couplings_fragment

    # The overall yield of the NCL strategy is limited by the synthesis of the
    # peptide fragment and the efficiency of the subsequent ligation step.
    # We assume one fragment is the limiting reagent.
    overall_ncl_yield = fragment_synthesis_yield * ligation_yield

    print("This strategy involves synthesizing two shorter fragments (e.g., 50aa each) and then joining them.")
    print(f"The yield of synthesizing a single {fragment_length}aa fragment is calculated first.")
    print(f"Fragment Synthesis Yield Equation: {coupling_efficiency}^{num_couplings_fragment}")
    print(f"Calculated Fragment Yield: {fragment_synthesis_yield:.4f} (or {fragment_synthesis_yield * 100:.2f}%)")
    print("\nNext, we calculate the overall process yield after ligation.")
    print(f"Assumed Ligation Step Yield: {ligation_yield * 100:.1f}%")
    print("Overall NCL Yield Equation: (Fragment Synthesis Yield) * (Ligation Yield)")
    print(f"Equation with numbers: {fragment_synthesis_yield:.4f} * {ligation_yield:.2f}")
    print(f"Calculated Overall Yield: {overall_ncl_yield:.4f} (or {overall_ncl_yield * 100:.2f}%)")
    print("-" * 60)

    # --- Conclusion ---
    print("\n--- Conclusion ---")
    print(f"Comparing the final theoretical yields:")
    print(f"  - Direct SPPS: {direct_spps_yield * 100:.2f}%")
    print(f"  - NCL Strategy: {overall_ncl_yield * 100:.2f}%")
    print("\nThe NCL strategy is significantly more efficient. It produces a higher yield of the target peptide and simplifies purification by dealing with smaller, purer fragments.")
    print("Therefore, Native Chemical Ligation is the most helpful technique for this task.")

if __name__ == "__main__":
    calculate_synthesis_yields()
