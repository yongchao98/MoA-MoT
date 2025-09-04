import sympy
import itertools

def check_partition_function():
    """
    This function calculates the partition function for the given 3-spin system from first principles
    and compares it to the provided answer.
    """
    # Define symbolic variables for the expression
    J, beta = sympy.symbols('J beta')

    # --- Step 1: Derive the correct partition function from first principles ---

    # Define the possible spin values
    spin_values = [+1, -1]

    # Generate all 2^3 = 8 possible states for the three spins (S1, S2, S3)
    all_states = list(itertools.product(spin_values, repeat=3))

    # Dictionary to store the calculated energy levels and their degeneracies
    energy_degeneracy_map = {}

    # Iterate through all possible states to calculate their energy
    for state in all_states:
        s1, s2, s3 = state
        
        # Calculate the energy for the current state using the given formula
        # E = -J[ S1S2 + S1S3 + S2S3 ]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        
        # Count the degeneracy for each energy level
        energy_degeneracy_map[energy] = energy_degeneracy_map.get(energy, 0) + 1

    # Construct the theoretical partition function Z from the derived energy levels and degeneracies
    # Z = Σ_levels g * e^(-βE)
    Z_theoretical = 0
    for energy, degeneracy in energy_degeneracy_map.items():
        Z_theoretical += degeneracy * sympy.exp(-beta * energy)

    # Simplify the resulting symbolic expression
    Z_theoretical_simplified = sympy.simplify(Z_theoretical)

    # --- Step 2: Define the candidate answer to be checked ---
    # The provided answer is A: Z = 2*e^(3*J*beta) + 6*e^(-J*beta)
    Z_candidate = 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta)

    # --- Step 3: Compare the theoretical result with the candidate answer ---
    # The most reliable way to check for equality of symbolic expressions is to see if their difference simplifies to zero.
    if sympy.simplify(Z_theoretical_simplified - Z_candidate) == 0:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        reason = (
            "Incorrect. The provided answer does not match the partition function derived from first principles.\n\n"
            "Reasoning:\n"
            "1. The system has two distinct energy levels with the following degeneracies:\n"
        )
        for energy, degeneracy in energy_degeneracy_map.items():
            reason += f"   - Energy Level: {energy}, Degeneracy (number of states): {degeneracy}\n"
        
        reason += (
            f"\n2. Based on these levels, the correct partition function is Z = Σ g * e^(-βE), which calculates to:\n"
            f"   Z_correct = {Z_theoretical_simplified}\n\n"
            f"3. The provided answer was:\n"
            f"   Z_provided = {Z_candidate}\n\n"
            "These two expressions are not equivalent."
        )
        return reason

# Run the check and print the result
result = check_partition_function()
print(result)