def calculate_spps_yield(peptide_length, coupling_efficiency):
    """
    Calculates the theoretical maximum yield for Solid-Phase Peptide Synthesis (SPPS).

    Args:
        peptide_length (int): The number of amino acids in the peptide.
        coupling_efficiency (float): The efficiency of a single coupling step (e.g., 0.99 for 99%).
    """
    if not (0 < coupling_efficiency <= 1):
        raise ValueError("Coupling efficiency must be between 0 and 1.")
    if peptide_length < 2:
        # No coupling steps needed for a single amino acid
        return 1.0

    # The number of coupling reactions is one less than the peptide length
    num_couplings = peptide_length - 1

    # Overall yield is the coupling efficiency raised to the power of the number of couplings
    total_yield = coupling_efficiency ** num_couplings

    return total_yield, num_couplings

# --- Parameters ---
peptide_length_full = 100
peptide_length_fragment = 50
efficiency = 0.99

# --- Calculations ---
yield_full, couplings_full = calculate_spps_yield(peptide_length_full, efficiency)
yield_fragment, couplings_fragment = calculate_spps_yield(peptide_length_fragment, efficiency)

# --- Output ---
print("--- Comparing Theoretical SPPS Yields ---")
print(f"Assuming a per-step coupling efficiency of: {efficiency * 100}%\n")

print(f"1. Direct synthesis of a {peptide_length_full}aa peptide:")
print(f"   Equation: Yield = {efficiency} ^ {couplings_full}")
print(f"   Theoretical Yield: {yield_full:.4f} (or {yield_full:.2%})")
print("\n")

print(f"2. Synthesis of a {peptide_length_fragment}aa fragment (for NCL):")
print(f"   Equation: Yield = {efficiency} ^ {couplings_fragment}")
print(f"   Theoretical Yield: {yield_fragment:.4f} (or {yield_fragment:.2%})")
print("\n")

print("Conclusion: The yield for a 50aa fragment is significantly higher than for a 100aa peptide,")
print("making Native Chemical Ligation (NCL) a much more viable strategy.")
