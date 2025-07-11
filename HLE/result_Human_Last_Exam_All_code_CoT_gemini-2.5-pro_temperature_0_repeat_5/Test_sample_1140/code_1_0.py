def calculate_spps_yield():
    """
    Calculates the theoretical maximum yield for direct Solid-Phase Peptide Synthesis (SPPS)
    of a long peptide, demonstrating why this method is inefficient for sequences of ~100aa.
    """
    peptide_length = 100
    # For a peptide of N amino acids, there are N-1 peptide bonds to be formed (couplings).
    number_of_couplings = peptide_length - 1
    # A high, yet realistic, coupling efficiency for each step in SPPS.
    coupling_efficiency = 0.995  # 99.5%

    # The final yield is the efficiency raised to the power of the number of coupling steps.
    final_yield = coupling_efficiency ** number_of_couplings

    print("This script calculates the theoretical maximum yield for direct Solid-Phase Peptide Synthesis (SPPS).")
    print("The result shows why synthesizing a 100aa peptide this way is very challenging.")
    print("-" * 30)
    print("Calculation Steps:")
    # The prompt asks to output each number in the final equation.
    print(f"Final Yield = (Coupling Efficiency) ^ (Number of Couplings)")
    print(f"Final Yield = ({coupling_efficiency}) ^ ({number_of_couplings})")
    print(f"Final Yield = {final_yield:.4f}")
    print("-" * 30)
    print(f"Result: Even with a high {coupling_efficiency:.1%} efficiency at each step, the maximum possible yield is only {final_yield:.2%}.")

calculate_spps_yield()