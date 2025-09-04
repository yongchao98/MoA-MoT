import sympy
import itertools
from collections import Counter

def check_partition_function():
    """
    Calculates the partition function for the 3-spin system and checks it against the given answer D.
    """
    # 1. Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # 2. Generate all possible microstates
    # Each of the 3 spins can be +1 or -1.
    spin_values = [-1, 1]
    all_states = list(itertools.product(spin_values, repeat=3))
    
    # Check if the total number of states is correct (2^3 = 8)
    if len(all_states) != 8:
        return f"Constraint violated: The number of microstates should be 8, but was {len(all_states)}."

    # 3. Calculate the energy for each state and find the degeneracies
    energies = []
    for state in all_states:
        s1, s2, s3 = state
        # Energy E = -J[ S1S2 + S1S3 + S2S3 ]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        energies.append(energy)
    
    # The Counter will group identical energies and count their occurrences (degeneracy)
    energy_degeneracy = Counter(energies)

    # 4. Construct the partition function Z from the energy levels and degeneracies
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # 5. Define the expression from the given answer (D)
    # Answer D is Z = 2*e^(3Jβ) + 6*e^(-Jβ)
    answer_D_expr = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 6. Compare the calculated Z with the answer's expression
    # sympy.simplify() ensures that the expressions are compared in their simplest form.
    # If the expressions are identical, their difference simplifies to 0.
    if sympy.simplify(calculated_Z - answer_D_expr) == 0:
        return "Correct"
    else:
        # If they don't match, provide a clear reason.
        # This helps diagnose if the calculation or the provided answer is wrong.
        expected_degeneracies = {
            -3*J: 2,
            J: 6
        }
        if dict(energy_degeneracy) != dict(expected_degeneracies):
             return f"The calculated energy levels or degeneracies are incorrect. Calculated: {dict(energy_degeneracy)}. Expected for answer D: {dict(expected_degeneracies)}."

        return (f"The final expression for the partition function is incorrect.\n"
                f"Calculated Z = {calculated_Z}\n"
                f"Answer D's Z = {answer_D_expr}")

# Execute the check
result = check_partition_function()
print(result)