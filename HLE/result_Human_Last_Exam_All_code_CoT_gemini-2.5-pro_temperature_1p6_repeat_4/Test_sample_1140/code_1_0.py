def calculate_spps_yield(peptide_length, coupling_efficiency):
    """
    Calculates the theoretical maximum yield for solid-phase peptide synthesis.

    Args:
        peptide_length (int): The number of amino acids in the peptide.
        coupling_efficiency (float): The efficiency of a single coupling step (e.g., 0.99 for 99%).
    """
    if not (0 < coupling_efficiency <= 1.0):
        print("Error: Coupling efficiency must be between 0 and 1.0.")
        return

    # The number of coupling reactions is the length of the peptide minus one.
    num_couplings = peptide_length - 1

    # The final yield is the efficiency raised to the power of the number of couplings.
    final_yield = coupling_efficiency ** num_couplings

    print("--- SPPS Yield Calculation ---")
    print(f"Peptide Length: {peptide_length} aa")
    print(f"Assumed per-step coupling efficiency: {coupling_efficiency * 100:.1f}%")
    print("\nThe theoretical maximum yield is calculated as (efficiency) ^ (number of couplings).")
    
    # Printing each number in the final equation
    print("\nFinal Equation:")
    print(f"Yield = {coupling_efficiency} ^ {num_couplings}")
    
    print(f"\nResult:")
    print(f"Theoretical Maximum Yield = {final_yield * 100:.4f}%")
    print("\nThis low theoretical yield demonstrates why continuous SPPS is not feasible for a 100aa peptide, making Native Chemical Ligation the preferred method.")


# --- Parameters for the 100aa peptide ---
# A very optimistic coupling efficiency for an automated synthesizer is 99.0%
length_of_peptide = 100
efficiency_per_step = 0.99

calculate_spps_yield(length_of_peptide, efficiency_per_step)