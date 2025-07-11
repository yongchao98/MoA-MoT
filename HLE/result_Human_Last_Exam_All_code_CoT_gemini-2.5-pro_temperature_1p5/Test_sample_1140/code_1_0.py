import math

def calculate_spps_yield():
    """
    Calculates the theoretical maximum yield for synthesizing a peptide
    via Solid-Phase Peptide Synthesis (SPPS) at various coupling efficiencies.
    """
    # Define the length of the peptide
    peptide_length = 100
    
    # The number of coupling reactions is peptide_length - 1
    num_couplings = peptide_length - 1
    
    # List of different per-step coupling efficiencies to test
    coupling_efficiencies = [98.0, 99.0, 99.5, 99.9]

    print(f"Calculating theoretical yield for a {peptide_length}-amino-acid peptide after {num_couplings} coupling steps:\n")

    for efficiency in coupling_efficiencies:
        # Convert percentage to decimal for calculation
        efficiency_decimal = efficiency / 100.0
        
        # Calculate the total theoretical yield
        total_yield = math.pow(efficiency_decimal, num_couplings)
        
        # Print the final equation and the result
        # The numbers in the equation are explicitly shown in the output string
        print(f"With a {efficiency:.1f}% per-step efficiency:")
        print(f"The final equation is: ({efficiency_decimal}) ^ {num_couplings}")
        print(f"Resulting maximum yield: {total_yield:.2%}\n")

if __name__ == '__main__':
    calculate_spps_yield()