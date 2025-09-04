import itertools
from collections import Counter
import sympy

def check_partition_function_correctness():
    """
    This function calculates the partition function for a system of three spins
    and checks if it matches the provided answer.
    """
    # Define J and beta as symbolic variables for the calculation
    J, beta = sympy.symbols('J beta')

    # 1. Generate all possible microstates for the three spins (+1 or -1)
    spin_values = [1, -1]
    num_spins = 3
    all_states = list(itertools.product(spin_values, repeat=num_spins))

    # Check if the total number of states is correct (2^3 = 8)
    if len(all_states) != 8:
        return f"Constraint check failed: The total number of microstates should be 2^{num_spins} = 8, but the code generated {len(all_states)}."

    # 2. Calculate the energy for each state
    energies = []
    for state in all_states:
        s1, s2, s3 = state
        # Energy function E = -J[S1*S2 + S1*S3 + S2*S3]
        energy = -J * (s1 * s2 + s1 * s3 + s2 * s3)
        energies.append(energy)

    # 3. Find the degeneracy of each energy level
    # Counter will create a dictionary of {energy: count}
    energy_degeneracy = Counter(energies)

    # 4. Construct the partition function Z from the calculated degeneracies
    Z_calculated = 0
    for energy, degeneracy in energy_degeneracy.items():
        Z_calculated += degeneracy * sympy.exp(-beta * energy)

    # 5. Define the partition function from the provided answer (Option B)
    # Z_B = 2 * e^(3Jβ) + 6 * e^(-Jβ)
    Z_answer_B = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 6. Compare the calculated Z with the answer's Z
    # The most reliable way to check for symbolic equality is to simplify their difference.
    # If the difference simplifies to zero, they are equivalent.
    if sympy.simplify(Z_calculated - Z_answer_B) == 0:
        return "Correct"
    else:
        # If they don't match, construct a detailed error message.
        # Let's check the calculated degeneracies against the logic in the answer.
        # The answer states there are 2 states with energy -3J and 6 states with energy J.
        expected_degeneracy = {-3*J: 2, J: 6}
        if energy_degeneracy == expected_degeneracy:
             # This means the calculation is correct, but the final expression in the answer is wrong.
             # (This case shouldn't happen if the answer is B, but it's good practice to check).
             reason = (f"The provided answer is incorrect. "
                       f"The calculated partition function is Z = {Z_calculated}, "
                       f"which does not simplify to the expression in option B: {Z_answer_B}.")
        else:
             # This means the energy calculation or degeneracy counting in the provided answer is flawed.
             reason = (f"The provided answer's reasoning is incorrect. "
                       f"The calculated energy levels and degeneracies are {dict(energy_degeneracy)}, "
                       f"which does not match the degeneracies described in the step-by-step analysis. "
                       f"The correct partition function is Z = {Z_calculated}.")
        return reason

# Run the check and print the result
result = check_partition_function_correctness()
print(result)