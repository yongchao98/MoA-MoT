import sympy
from itertools import product
from collections import Counter

def check_correctness():
    """
    This function verifies the correctness of the provided answer for a statistical mechanics problem.
    It calculates the partition function for a system of three spins from first principles
    and compares it to the given answer options.

    Question:
    Consider a system of three spins S1, S2 and S3. Each of which can take spin +1 and -1.
    The energy of the system is given by, E = -J[ S1S2 + S1S3 + S2S3 ].
    Find the partition function Z of the system. (β = 1/kT)

    Options:
    A) Z = 6*e^(2*J*β) + 2*e^(-2*J*β)
    B) Z = 2*e^(2*J*β) + 6*e^(-2*J*β)
    C) Z = 2*e^(-3*J*β) + 6*e^(J*β)
    D) Z = 2*e^(3*J*β) + 6*e^(-J*β)

    The provided answer to check is D.
    """

    # Define symbolic variables to handle the math abstractly
    J, beta = sympy.symbols('J beta')

    # The answer given by the LLM is D
    answer_to_check = "D"
    options = {
        "A": 6 * sympy.exp(2*J*beta) + 2 * sympy.exp(-2*J*beta),
        "B": 2 * sympy.exp(2*J*beta) + 6 * sympy.exp(-2*J*beta),
        "C": 2 * sympy.exp(-3*J*beta) + 6 * sympy.exp(J*beta),
        "D": 2 * sympy.exp(3*J*beta) + 6 * sympy.exp(-J*beta)
    }
    answer_expression = options[answer_to_check]

    # --- Step 1: Enumerate all possible microstates ---
    # Each of the 3 spins can be +1 or -1, so there are 2^3 = 8 states.
    spin_values = [1, -1]
    all_states = list(product(spin_values, repeat=3))

    # --- Step 2: Calculate the energy for each state and find degeneracies ---
    # We will store the energy of each state in a list.
    energy_levels = []
    for s1, s2, s3 in all_states:
        # Energy formula: E = -J[S1*S2 + S1*S3 + S2*S3]
        energy = -J * (s1*s2 + s1*s3 + s2*s3)
        energy_levels.append(energy)

    # Use a Counter to find the degeneracy (g_i) of each unique energy level (E_i)
    energy_degeneracy = Counter(energy_levels)

    # --- Step 3: Construct the partition function Z from the calculated energies ---
    # The partition function is the sum of Boltzmann factors over all states:
    # Z = Σ_i g_i * exp(-β * E_i)
    Z_calculated = sum(g * sympy.exp(-beta * E) for E, g in energy_degeneracy.items())

    # --- Step 4: Compare the calculated Z with the expression from the given answer ---
    # We use sympy.simplify to ensure that two mathematically equivalent expressions
    # are recognized as equal, regardless of their form.
    if sympy.simplify(Z_calculated - answer_expression) == 0:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed explanation.
        explanation_parts = []
        for E, g in sorted(energy_degeneracy.items(), key=lambda item: str(item[0])):
             explanation_parts.append(f"  - Energy E = {E} has a degeneracy of {g}.")
        degeneracy_explanation = "\n".join(explanation_parts)

        reason = (
            f"The provided answer {answer_to_check} is incorrect.\n"
            f"The expression for answer {answer_to_check} is: Z = {answer_expression}\n\n"
            f"The correct derivation is as follows:\n"
            f"The system has 8 microstates. By calculating the energy E = -J[S1*S2 + S1*S3 + S2*S3] for each state, we find two distinct energy levels:\n"
            f"{degeneracy_explanation}\n\n"
            f"The partition function is Z = Σ g_i * exp(-β*E_i), which gives:\n"
            f"Z = {next(g for E, g in energy_degeneracy.items() if E == -3*J)}*exp(-β*({-3*J})) + {next(g for E, g in energy_degeneracy.items() if E == J)}*exp(-β*({J}))\n"
            f"This simplifies to the correct expression:\n"
            f"Z = {Z_calculated}\n\n"
            f"The provided answer does not match this result."
        )
        return reason

# Run the check and print the result
result = check_correctness()
print(result)