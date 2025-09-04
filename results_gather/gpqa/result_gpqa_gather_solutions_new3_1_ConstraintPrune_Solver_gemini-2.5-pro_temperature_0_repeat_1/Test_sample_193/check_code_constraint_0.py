import itertools
import collections

def check_partition_function():
    """
    Calculates the partition function for the 3-spin system from first principles
    and checks it against the provided answer D.
    """
    # The final answer from the LLM is D.
    # D) Z = 2e^(3Jβ) + 6e^(-Jβ)
    # This mathematical form implies two distinct energy levels and their degeneracies:
    # 1. A term with coefficient 2 and exponent 3Jβ -> degeneracy=2, energy E=-3J.
    # 2. A term with coefficient 6 and exponent -Jβ -> degeneracy=6, energy E=+J.
    # We can represent this as a dictionary mapping {Energy_Coefficient: Degeneracy}
    expected_degeneracies = {-3: 2, 1: 6}

    # Step 1: Generate all 2^3 = 8 possible microstates for the three spins.
    spins = [-1, 1]
    all_states = list(itertools.product(spins, repeat=3))

    # Sanity check for the number of states.
    if len(all_states) != 8:
        return "Constraint check failed: The total number of states should be 2^3 = 8."

    # Step 2 & 3: Calculate the energy for each state and count the degeneracies.
    # We use a dictionary to store the count of states for each energy level.
    # The key will be the energy coefficient (E/J).
    energy_counts = collections.defaultdict(int)
    for s1, s2, s3 in all_states:
        # The energy is E = -J * (s1*s2 + s1*s3 + s2*s3).
        # We calculate the energy coefficient, E/J.
        energy_coeff = -(s1*s2 + s1*s3 + s2*s3)
        energy_counts[energy_coeff] += 1
    
    calculated_degeneracies = dict(energy_counts)

    # Step 4: Compare the calculated result with the expected result from answer D.
    if calculated_degeneracies == expected_degeneracies:
        return "Correct"
    else:
        reason = (
            "Incorrect. The calculated energy levels and/or degeneracies do not match "
            "those implied by the provided answer D.\n"
            f"Answer D implies: {expected_degeneracies}\n"
            f"Calculation from first principles gives: {calculated_degeneracies}"
        )
        return reason

# Execute the check and print the result.
result = check_partition_function()
print(result)