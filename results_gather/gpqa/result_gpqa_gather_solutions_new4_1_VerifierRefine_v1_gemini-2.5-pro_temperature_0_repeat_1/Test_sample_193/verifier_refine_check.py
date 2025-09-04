import sympy

def check_partition_function_correctness():
    """
    Calculates the partition function for the three-spin system from first principles
    and compares it to the expression from the provided answer.
    """
    # Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # The system has three spins, each can be +1 or -1
    spins = [+1, -1]
    
    # Initialize the partition function Z to zero
    Z_calculated = 0
    
    # A dictionary to store energy levels and their degeneracies for verification
    energy_levels = {}

    # Iterate through all 2^3 = 8 possible states (s1, s2, s3)
    for s1 in spins:
        for s2 in spins:
            for s3 in spins:
                # Calculate the energy for the current state using the Hamiltonian
                # E = -J[ S1*S2 + S1*S3 + S2*S3 ]
                energy = -J * (s1*s2 + s1*s3 + s2*s3)
                
                # Add the Boltzmann factor for this state to the total partition function
                Z_calculated += sympy.exp(-beta * energy)

                # Count the degeneracy for each energy level
                if energy in energy_levels:
                    energy_levels[energy] += 1
                else:
                    energy_levels[energy] = 1

    # The calculated partition function is the sum of all terms.
    # Sympy automatically groups identical terms, but we can simplify for a canonical form.
    Z_calculated = sympy.simplify(Z_calculated)

    # The provided answer is C, which corresponds to the expression:
    # Z = 2 * e^(3*J*beta) + 6 * e^(-J*beta)
    Z_answer_C = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # Check if the calculated partition function matches the expression from answer C.
    # A robust way to check for symbolic equality is to see if their difference simplifies to zero.
    if sympy.simplify(Z_calculated - Z_answer_C) == 0:
        # As an additional check, verify the intermediate steps (degeneracies)
        # The correct derivation should find two energy levels: -3J (degeneracy 2) and +J (degeneracy 6)
        expected_degeneracies = {-3*J: 2, J: 6}
        if energy_levels == expected_degeneracies:
            return "Correct"
        else:
            # This case is unlikely but would indicate a correct final answer derived from flawed intermediate steps.
            return f"Incorrect. The final expression for Z is correct, but the derivation of energy level degeneracies in the provided answer is flawed. The code calculated degeneracies as: {energy_levels}"
    else:
        return f"Incorrect. The provided answer is C, which evaluates to Z = {Z_answer_C}. However, the correct partition function calculated from first principles is Z = {Z_calculated}."

# Run the check
result = check_partition_function_correctness()
print(result)