import numpy as np

def check_partition_function():
    """
    Checks the correctness of the proposed partition function formula by comparing it
    to a direct calculation from first principles.
    """

    # The formula from the proposed answer (Option C)
    # Z = 2 * e^(3Jβ) + 6 * e^(-Jβ)
    def z_formula_C(J, beta):
        return 2 * np.exp(3 * J * beta) + 6 * np.exp(-J * beta)

    # Function to calculate Z by summing over all 8 microstates
    def z_first_principles(J, beta):
        z_sum = 0
        spins = [-1, 1]
        # Iterate through all 2^3 = 8 possible spin configurations
        for s1 in spins:
            for s2 in spins:
                for s3 in spins:
                    # Calculate the energy for the current state
                    energy = -J * (s1 * s2 + s1 * s3 + s2 * s3)
                    # Add the Boltzmann factor for this state to the sum
                    z_sum += np.exp(-beta * energy)
        return z_sum

    # We will test for a few arbitrary non-trivial values of J and beta
    test_cases = [
        {'J': 1.0, 'beta': 1.0},
        {'J': -0.5, 'beta': 2.0},
        {'J': 0.0, 'beta': 1.0},  # Special case: J=0, all energies are 0, Z should be 8
        {'J': 1.0, 'beta': 0.0}   # Special case: beta=0 (infinite T), Z should be 8
    ]

    for params in test_cases:
        J = params['J']
        beta = params['beta']
        
        # Calculate Z using both methods
        z_from_formula = z_formula_C(J, beta)
        z_from_summation = z_first_principles(J, beta)

        # Compare the results using a small tolerance for floating point numbers
        if not np.isclose(z_from_formula, z_from_summation):
            return (f"Incorrect. The provided formula from option C is wrong.\n"
                    f"For test case J={J}, beta={beta}:\n"
                    f" - The formula gives Z = {z_from_formula}\n"
                    f" - The correct value from summing all states is Z = {z_from_summation}")

    # If all test cases pass, the formula is correct.
    # Let's double-check the derivation of the energy levels and their degeneracies.
    # E = -J * (s1*s2 + s1*s3 + s2*s3)
    # Case 1: All spins aligned (+1,+1,+1) or (-1,-1,-1)
    # s1*s2+s1*s3+s2*s3 = 1+1+1 = 3. Energy E = -3J. Degeneracy = 2.
    # Term: 2 * exp(-beta * -3J) = 2 * exp(3*J*beta)
    # Case 2: One spin flipped, e.g., (+1,+1,-1) or (-1,-1,+1)
    # s1*s2+s1*s3+s2*s3 = 1-1-1 = -1. Energy E = J. Degeneracy = 6 (3 with one spin down, 3 with one spin up).
    # Term: 6 * exp(-beta * J)
    # Total Z = 2*exp(3*J*beta) + 6*exp(-J*beta).
    # The derivation matches the formula and the numerical check.
    
    return "Correct"

# Execute the check and print the result
result = check_partition_function()
print(result)