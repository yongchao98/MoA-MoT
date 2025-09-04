import itertools
from collections import Counter

def check_correctness():
    """
    This function checks the correctness of the provided answer for the partition function
    of a 3-spin system. It does this by recalculating the result from first principles.

    The energy of the system is E = -J[S1*S2 + S1*S3 + S2*S3].
    The partition function is Z = sum(exp(-beta*E)) over all states.

    The provided answer is A: Z = 2*e^(3*J*beta) + 6*e^(-J*beta).
    This implies:
    - An energy level E1 with degeneracy g1=2, where -beta*E1 = 3*J*beta => E1 = -3J.
    - An energy level E2 with degeneracy g2=6, where -beta*E2 = -J*beta => E2 = J.
    
    The code will verify these energy levels and degeneracies.
    """
    
    # 1. Define the system and generate all possible states.
    spin_values = [-1, 1]
    num_spins = 3
    all_states = list(itertools.product(spin_values, repeat=num_spins))

    # Constraint Check: The total number of states must be 2^3 = 8.
    if len(all_states) != 8:
        return f"Constraint check failed: The number of generated states is {len(all_states)}, but it should be 8."

    # 2. Calculate the energy for each state.
    # We will calculate the energy coefficient (E/J) for simplicity.
    # E/J = -(S1*S2 + S1*S3 + S2*S3)
    energy_coeffs = []
    for state in all_states:
        s1, s2, s3 = state
        coeff = -(s1 * s2 + s1 * s3 + s2 * s3)
        energy_coeffs.append(coeff)

    # 3. Determine the unique energy levels and their degeneracies.
    calculated_degeneracies = Counter(energy_coeffs)
    
    # 4. Define the expected energy levels and degeneracies based on the provided answer 'A'.
    # The keys are the energy coefficients (E/J) and values are the degeneracies.
    expected_degeneracies = {
        -3: 2,  # Corresponds to E = -3J, degeneracy = 2
         1: 6   # Corresponds to E = +J,  degeneracy = 6
    }

    # 5. Compare the calculated results with the expected results.
    # The `Counter` objects will be equal if they have the same items and counts.
    if calculated_degeneracies == expected_degeneracies:
        return "Correct"
    else:
        # If they don't match, provide a detailed reason.
        # Check if the number of unique energy levels is wrong.
        if len(calculated_degeneracies) != len(expected_degeneracies):
            return (f"Incorrect. The calculation resulted in {len(calculated_degeneracies)} unique energy levels, "
                    f"but the answer implies {len(expected_degeneracies)} levels. "
                    f"Calculated energy coefficients (E/J): {list(calculated_degeneracies.keys())}.")
        
        # If the number of levels is correct, the energies or degeneracies must be wrong.
        return (f"Incorrect. The calculated energy levels or their degeneracies do not match the answer. "
                f"Calculated (Energy_Coeff: Degeneracy): {dict(calculated_degeneracies)}. "
                f"Expected based on answer A: {expected_degeneracies}.")

# Run the check and print the result.
result = check_correctness()
print(result)