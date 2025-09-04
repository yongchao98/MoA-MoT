import itertools
from collections import Counter

def check_partition_function():
    """
    Calculates the exact partition function for the three-spin system and
    compares it to the provided answer.
    """
    # Define the possible spin values
    spins = [+1, -1]
    num_spins = 3

    # Generate all 2^3 = 8 possible microstates
    all_states = list(itertools.product(spins, repeat=num_spins))

    # Calculate the energy for each state and store the energy coefficients (E/J)
    energy_coeffs = []
    for state in all_states:
        S1, S2, S3 = state
        # Energy E = -J * (S1*S2 + S1*S3 + S2*S3)
        # We calculate the coefficient E/J
        energy_coeff = -(S1 * S2 + S1 * S3 + S2 * S3)
        energy_coeffs.append(energy_coeff)

    # Count the degeneracy of each energy level
    # The result is a dictionary mapping {energy_coeff: degeneracy}
    calculated_degeneracies = Counter(energy_coeffs)

    # The provided answer is D) Z = 2 * e^(3Jβ) + 6 * e^(-Jβ)
    # Let's analyze the terms in the answer to find the expected degeneracies.
    # The general form of a term in the partition function is g * exp(-βE),
    # where g is the degeneracy.
    # This can be written as g * exp(-(E/J) * Jβ).

    # Term 1: 2 * e^(3Jβ)
    # Degeneracy g1 = 2.
    # Exponent is 3Jβ, so -(E1/J) * Jβ = 3Jβ => E1/J = -3.
    
    # Term 2: 6 * e^(-Jβ)
    # Degeneracy g2 = 6.
    # Exponent is -Jβ, so -(E2/J) * Jβ = -Jβ => E2/J = 1.

    # So, the expected degeneracies from answer D are:
    # Energy E = -3J (coeff = -3) should have degeneracy 2.
    # Energy E = +J  (coeff =  1) should have degeneracy 6.
    expected_degeneracies = {-3: 2, 1: 6}

    # Check if the calculated degeneracies match the expected ones.
    # The Counter object behaves like a dictionary for comparison.
    if calculated_degeneracies == expected_degeneracies:
        return "Correct"
    else:
        # If they don't match, provide a detailed reason.
        reason = (
            "The answer is incorrect.\n"
            f"The provided answer implies the following energy levels and degeneracies: {expected_degeneracies}.\n"
            f"The exact calculation for all 8 states yields these energy levels and degeneracies: {dict(calculated_degeneracies)}.\n"
            "The calculated results do not match the degeneracies derived from the answer."
        )
        # Check if the total number of states is correct
        if sum(calculated_degeneracies.values()) != 8:
            reason += f"\nAdditionally, the calculated total number of states is {sum(calculated_degeneracies.values())}, which should be 8."
        if sum(expected_degeneracies.values()) != 8:
             reason += f"\nAdditionally, the total number of states implied by the answer is {sum(expected_degeneracies.values())}, which should be 8."
        return reason

# Run the check
result = check_partition_function()
print(result)