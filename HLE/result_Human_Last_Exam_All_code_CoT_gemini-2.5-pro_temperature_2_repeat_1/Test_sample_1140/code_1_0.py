def calculate_spps_yield(peptide_length, coupling_efficiency):
    """
    Calculates and prints the theoretical maximum yield for Solid-Phase Peptide Synthesis (SPPS).

    Args:
        peptide_length (int): The number of amino acids in the peptide.
        coupling_efficiency (float): The efficiency of each coupling step (e.g., 0.99 for 99%).
    """
    if peptide_length < 2:
        print("Peptide length must be at least 2 for a coupling reaction.")
        return

    # The number of coupling reactions is the peptide length minus one.
    num_couplings = peptide_length - 1

    # Calculate the overall theoretical yield.
    overall_yield = coupling_efficiency ** num_couplings

    print(f"For a {peptide_length}aa peptide with a {coupling_efficiency*100:.1f}% per-step efficiency:")
    print("The final yield is calculated as: (Step Yield) ^ (Number of Couplings)")
    print(f"Overall Yield = {coupling_efficiency} ^ {num_couplings}")
    print(f"Result: {overall_yield:.4f} (or {overall_yield*100:.2f}%)\n")

# --- Analysis ---
print("--- Yield Analysis for Synthesizing a 100aa Peptide ---\n")

print("Strategy 1: Direct SPPS synthesis of the full 100aa peptide.")
# Assuming a very good, but realistic, 99.0% coupling efficiency per step.
calculate_spps_yield(peptide_length=100, coupling_efficiency=0.99)

print("-----------------------------------------------------------\n")

print("Strategy 2: SPPS synthesis of a 50aa fragment for Native Chemical Ligation.")
print("This represents the yield for *each* of the two fragments needed for NCL.")
# Using the same 99.0% efficiency for comparison.
calculate_spps_yield(peptide_length=50, coupling_efficiency=0.99)

print("The yield for the 50aa fragment is significantly higher, making the NCL approach far more practical.")
