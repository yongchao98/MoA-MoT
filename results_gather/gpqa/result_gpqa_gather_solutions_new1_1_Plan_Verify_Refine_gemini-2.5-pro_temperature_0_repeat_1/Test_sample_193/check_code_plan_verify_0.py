import sympy
from collections import defaultdict
import itertools

def check_correctness():
    """
    Checks the correctness of the partition function calculation for a 3-spin system.
    """
    # 1. Define symbolic variables
    J, beta = sympy.symbols('J beta')

    # 2. Generate all 8 possible microstates
    spin_values = [1, -1]
    all_states = list(itertools.product(spin_values, repeat=3))

    # Dictionary to store energy levels and their degeneracies
    energy_degeneracy = defaultdict(int)

    # 3. & 4. Calculate energy for each state and find degeneracies
    for s1, s2, s3 in all_states:
        # Energy formula: E = -J[S1*S2 + S1*S3 + S2*S3]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        energy_degeneracy[energy] += 1

    # Verify the intermediate step: energy levels and degeneracies
    # Expected: E1 = -3J (g=2), E2 = +J (g=6)
    expected_energies = {-3*J: 2, J: 6}
    if energy_degeneracy != expected_energies:
        return (f"Incorrect intermediate calculation. The energy levels and/or degeneracies are wrong.\n"
                f"Calculated: {dict(energy_degeneracy)}\n"
                f"Expected: {expected_energies}")

    # 5. Construct the partition function from the calculated energies
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        # Z = Σ gᵢ * e^(-βEᵢ)
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # 6. Create a symbolic expression for the given answer (Option B)
    # The final answer provided is <<<B>>>, which corresponds to Z = 2e^(3Jβ) + 6e^(-Jβ)
    answer_Z_symbolic = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 7. Compare the derived partition function with the answer's expression
    # sympy.simplify() ensures that two mathematically equivalent but differently
    # structured expressions are compared correctly.
    if sympy.simplify(calculated_Z - answer_Z_symbolic) == 0:
        return "Correct"
    else:
        return (f"The final expression is incorrect.\n"
                f"Derived partition function: Z = {sympy.simplify(calculated_Z)}\n"
                f"Expression from the answer (Option B): Z = {sympy.simplify(answer_Z_symbolic)}")

# Run the check
result = check_correctness()
print(result)