import math

def calculate_spps_yield(peptide_length, coupling_efficiency):
    """
    Calculates the theoretical maximum yield for Solid-Phase Peptide Synthesis (SPPS).

    Args:
        peptide_length (int): The number of amino acids in the peptide.
        coupling_efficiency (float): The efficiency of each coupling step (e.g., 0.99 for 99%).
    """
    if not (0 < coupling_efficiency <= 1.0):
        print("Error: Coupling efficiency must be between 0 and 1.0.")
        return

    # For a peptide of N amino acids, there are N-1 coupling steps.
    num_couplings = peptide_length - 1
    
    # The overall yield is the efficiency raised to the power of the number of steps.
    overall_yield = math.pow(coupling_efficiency, num_couplings)
    
    print(f"For a {peptide_length} amino acid peptide:")
    print(f"With a per-step coupling efficiency of {coupling_efficiency * 100:.1f}%, there are {num_couplings} coupling steps.")
    print("The theoretical maximum yield is calculated by the equation:")
    # The prompt requires printing each number in the final equation.
    print(f"    Yield = (Coupling Efficiency) ^ (Number of Couplings)")
    print(f"    {overall_yield:.4f} = {coupling_efficiency} ^ {num_couplings}")
    print(f"This results in a final theoretical yield of {overall_yield:.2%}.\n")

# --- Main execution ---
# Demonstrate with a typical high efficiency
calculate_spps_yield(peptide_length=100, coupling_efficiency=0.99)

# Demonstrate with a slightly lower, more realistic efficiency
calculate_spps_yield(peptide_length=100, coupling_efficiency=0.98)

# Demonstrate the dramatic effect on an even longer peptide
calculate_spps_yield(peptide_length=150, coupling_efficiency=0.99)
