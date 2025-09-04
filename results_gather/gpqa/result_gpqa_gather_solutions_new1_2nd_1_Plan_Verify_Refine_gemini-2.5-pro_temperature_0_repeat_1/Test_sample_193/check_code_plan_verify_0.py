import sympy

def check_answer():
    """
    This function checks the correctness of the provided answer for the partition function Z.
    It calculates Z from first principles and compares it to the given answer.
    """
    # Define symbols for the expression
    J, beta = sympy.symbols('J beta')
    exp = sympy.exp

    # The final answer from the LLM is 'D', which corresponds to the expression:
    # Z = 2 * e^(3Jβ) + 6 * e^(-Jβ)
    llm_answer_expression = 2 * exp(3 * J * beta) + 6 * exp(-J * beta)

    # --- Step 1: Calculate the partition function from first principles ---

    # Define the possible spin values
    spins = [+1, -1]

    # Dictionary to store energies and their degeneracies
    energy_levels = {}

    # Iterate through all 2^3 = 8 possible states
    for s1 in spins:
        for s2 in spins:
            for s3 in spins:
                # Calculate the energy for the current state
                # E = -J[ S1S2 + S1S3 + S2S3 ]
                energy = -J * (s1 * s2 + s1 * s3 + s2 * s3)

                # Count the degeneracy for each energy level
                energy_levels[energy] = energy_levels.get(energy, 0) + 1

    # --- Step 2: Verify the intermediate constraints (energy levels and degeneracies) ---

    # Constraint 1: Total number of states must be 8
    total_degeneracy = sum(energy_levels.values())
    if total_degeneracy != 8:
        return f"Incorrect: The derivation failed to account for all 8 states. Found {total_degeneracy} states."

    # Constraint 2: Check the energy levels and their degeneracies
    # Expected from manual derivation: E1 = -3J (g=2), E2 = +J (g=6)
    expected_energies = {-3 * J: 2, J: 6}
    if energy_levels != expected_energies:
        # Format the dictionaries for clear comparison in the error message
        found_str = ", ".join([f"E={k}: g={v}" for k, v in energy_levels.items()])
        expected_str = ", ".join([f"E={k}: g={v}" for k, v in expected_energies.items()])
        return (f"Incorrect: The energy levels or their degeneracies are wrong.\n"
                f"Found: {found_str}\n"
                f"Expected: {expected_str}")

    # --- Step 3: Construct the partition function from the derived energy levels ---

    calculated_Z = 0
    for energy, degeneracy in energy_levels.items():
        # Z = sum(g_i * exp(-beta * E_i))
        calculated_Z += degeneracy * exp(-beta * energy)

    # --- Step 4: Compare the calculated result with the LLM's answer ---

    # Simplify both expressions to ensure a canonical form for comparison
    simplified_calculated_Z = sympy.simplify(calculated_Z)
    simplified_llm_answer = sympy.simplify(llm_answer_expression)

    # The check: is the difference zero?
    if sympy.simplify(simplified_calculated_Z - simplified_llm_answer) == 0:
        return "Correct"
    else:
        return (f"Incorrect: The provided answer expression does not match the correctly calculated partition function.\n"
                f"Provided answer expression: {simplified_llm_answer}\n"
                f"Correctly calculated expression: {simplified_calculated_Z}")

# Print the result of the check
print(check_answer())