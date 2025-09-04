import sympy
from itertools import product

def check_partition_function():
    """
    This function verifies the partition function for the given 3-spin system.
    It calculates the partition function from first principles by enumerating all
    possible states, calculating their energies, and summing their Boltzmann factors.
    It then compares this derived result with the expression from the selected answer.
    """
    # Define symbolic variables for the calculation
    J, beta = sympy.symbols('J beta')

    # Define the possible spin values
    spin_values = [+1, -1]

    # Generate all 2^3 = 8 possible spin configurations
    all_states = list(product(spin_values, repeat=3))

    # Dictionary to store energy levels and their degeneracies
    energy_degeneracy = {}

    # Calculate the energy for each state and count degeneracies
    for state in all_states:
        s1, s2, s3 = state
        # Energy function E = -J[S1S2 + S1S3 + S2S3]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        
        # Update the degeneracy count for this energy level
        energy_degeneracy[energy] = energy_degeneracy.get(energy, 0) + 1

    # Construct the partition function Z from the derived energy levels and degeneracies
    # Z = sum(g_i * exp(-beta * E_i))
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # The final answer provided by the LLM is <<<C>>>.
    # Let's define the expression for option C from the question.
    # C) Z = 2 e^(3Jβ) + 6e^(-Jβ)
    answer_C_expression = 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta)

    # Check if the calculated partition function is symbolically equal to the answer's expression.
    # sympy.simplify(expr1 - expr2) will be 0 if they are equivalent.
    if sympy.simplify(calculated_Z - answer_C_expression) == 0:
        return "Correct"
    else:
        # If they don't match, provide the reason.
        reason = (
            "The answer is incorrect.\n"
            f"The provided answer corresponds to the expression: {answer_C_expression}\n"
            f"However, the correct partition function calculated from first principles is: {sympy.simplify(calculated_Z)}\n"
            "The discrepancy arises from an incorrect expression for the partition function."
        )
        return reason

# Run the check
result = check_partition_function()
print(result)