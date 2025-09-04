import sympy
from collections import Counter
import itertools

def check_partition_function():
    """
    Calculates the partition function for the three-spin system and checks it against the given answer.
    """
    # 1. Define symbolic variables for J and beta
    J, beta = sympy.symbols('J beta')

    # 2. Define the possible spin values and generate all 8 microstates
    spin_values = [-1, 1]
    # All possible combinations of (S1, S2, S3)
    states = list(itertools.product(spin_values, repeat=3))

    # Check if the number of states is correct (2^3 = 8)
    if len(states) != 8:
        return f"Constraint check failed: The total number of states should be 8, but was {len(states)}."

    # 3. Calculate the energy for each state and count degeneracies
    energies = []
    for s1, s2, s3 in states:
        # Energy formula: E = -J[S1S2 + S1S3 + S2S3]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        energies.append(energy)
    
    # The Counter object will automatically group identical energies and count them
    energy_degeneracy = Counter(energies)

    # 4. Construct the partition function Z from the calculated energies and degeneracies
    calculated_Z = 0
    for energy, degeneracy in energy_degeneracy.items():
        # Z = sum(g(E) * exp(-beta * E))
        calculated_Z += degeneracy * sympy.exp(-beta * energy)

    # 5. Define the candidate answer as a symbolic expression.
    # The given answer is A) Z = 2e^(3Jβ) + 6e^(-Jβ)
    candidate_answer_Z = 2 * sympy.exp(3 * J * beta) + 6 * sympy.exp(-J * beta)

    # 6. Compare the calculated result with the candidate answer.
    # sympy.simplify() is used to ensure the comparison is robust.
    # If the expressions are equivalent, their difference simplifies to 0.
    if sympy.simplify(calculated_Z - candidate_answer_Z) == 0:
        return "Correct"
    else:
        # If they don't match, provide a detailed reason.
        reason = (
            "The calculated partition function does not match the provided answer.\n"
            f"Calculated Z = {calculated_Z}\n"
            f"Provided Answer Z = {candidate_answer_Z}\n"
            "Let's re-check the energy levels and degeneracies:\n"
            f"Calculated energy levels and their degeneracies: {dict(energy_degeneracy)}\n"
        )
        # Expected degeneracies based on the correct answer
        expected_energies = {-3*J: 2, J: 6}
        if dict(energy_degeneracy) != expected_energies:
            reason += f"The calculated degeneracies are incorrect. They should be {expected_energies}."
        else:
            reason += "The degeneracies were calculated correctly, but the final expression is still wrong, which indicates a potential error in constructing Z from the degeneracies."
        return reason

# Run the check
result = check_partition_function()
print(result)