def calculate_spps_yield(peptide_length, coupling_efficiency_percent):
    """
    Calculates and prints the theoretical maximum yield for Solid-Phase Peptide Synthesis (SPPS).
    
    This function demonstrates how cumulative efficiency losses make it difficult
    to synthesize long peptides directly.
    """
    # The number of coupling reactions is one less than the number of amino acids.
    num_couplings = peptide_length - 1
    
    # Convert the efficiency from a percentage to a decimal fraction.
    efficiency_decimal = coupling_efficiency_percent / 100.0
    
    # The overall theoretical yield is the efficiency raised to the power of the number of couplings.
    overall_yield_fraction = efficiency_decimal ** num_couplings
    
    # Convert the final yield to a percentage for easier interpretation.
    overall_yield_percent = overall_yield_fraction * 100
    
    # Print the parameters and the equation used for the calculation.
    print(f"Calculating for a {peptide_length}aa peptide with {coupling_efficiency_percent}% per-step efficiency:")
    print(f"The equation for overall yield is: (Per-Step Efficiency) ^ (Number of Couplings)")
    print(f"Calculation: {overall_yield_percent:.2f}% = ({efficiency_decimal}) ^ ({num_couplings}) * 100")
    print("-" * 40)

# --- Main Script ---
# Define the length of the peptide from the problem.
peptide_length = 100

print(f"### SPPS Yield Calculation for a {peptide_length}-Amino-Acid Peptide ###\n")
print("This script shows why direct synthesis of long peptides is inefficient,")
print("making a method like Native Chemical Ligation (NCL) a better choice.\n")

# Case 1: An extremely optimistic, near-perfect coupling efficiency.
calculate_spps_yield(peptide_length, 99.5)

# Case 2: A more standard, excellent coupling efficiency.
calculate_spps_yield(peptide_length, 99.0)

# Case 3: A good, but more realistic, coupling efficiency.
calculate_spps_yield(peptide_length, 98.0)

print("\nAs shown, even with excellent per-step efficiency, the final yield for a 100aa peptide is very low.")
print("NCL overcomes this by ligating shorter, high-yield fragments.")
