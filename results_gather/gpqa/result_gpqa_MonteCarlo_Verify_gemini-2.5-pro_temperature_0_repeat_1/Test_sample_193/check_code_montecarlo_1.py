import itertools
from collections import Counter
import numpy as np

def check_correctness():
    """
    This function checks the correctness of the given answer for the partition function.

    The problem involves a system of three spins S1, S2, S3 (+1 or -1).
    The energy is E = -J[S1*S2 + S1*S3 + S2*S3].
    The partition function is Z = sum(exp(-beta * E)) over all states.

    The function will:
    1. Enumerate all 2^3 = 8 possible states of the system.
    2. Calculate the energy for each state.
    3. Group states by energy level to find the degeneracy of each level.
    4. Construct the theoretical partition function from these degeneracies.
    5. Compare the derived function with the provided answer D.
    """
    # Define the possible values for a single spin
    spin_values = [-1, 1]

    # 1. Generate all 8 possible states for the three-spin system
    all_states = list(itertools.product(spin_values, repeat=3))

    # 2. Calculate the energy for each state.
    # For simplicity, we calculate the energy term E/(-J) = S1*S2 + S1*S3 + S2*S3
    energy_terms = []
    for state in all_states:
        s1, s2, s3 = state
        term = s1 * s2 + s1 * s3 + s2 * s3
        energy_terms.append(term)

    # 3. Count the number of states (degeneracy) for each energy level.
    # The keys of the counter are the energy terms (E / -J), and values are the degeneracies.
    degeneracy = Counter(energy_terms)

    # The partition function Z is sum_levels(g(E) * exp(-beta * E)).
    # With E = -J * term, this becomes Z = sum_terms(g(term) * exp(J * beta * term)).
    # The calculated degeneracies are expected to be:
    # - term = 3 (for E = -3J): 2 states ((+1,+1,+1) and (-1,-1,-1))
    # - term = -1 (for E = J): 6 states (all other states)
    
    # 4. Define the expected degeneracies based on our calculation and the provided answer D.
    # Answer D: Z = 2 * e^(3*J*beta) + 6 * e^(-J*beta)
    # This implies:
    # - A term exp(3*J*beta), which means energy term = 3, with a coefficient (degeneracy) of 2.
    # - A term exp(-J*beta), which means energy term = -1, with a coefficient (degeneracy) of 6.
    expected_degeneracy = {3: 2, -1: 6}

    # 5. Compare the calculated degeneracy with the one from the answer.
    if degeneracy == expected_degeneracy:
        return "Correct"
    else:
        # Construct the calculated partition function string for a clear error message
        calculated_z_parts = []
        # Sort by energy term descending to match the usual representation
        for term, count in sorted(degeneracy.items(), key=lambda x: x[0], reverse=True):
            # The exponent in Z is J*beta*term
            calculated_z_parts.append(f"{count}*exp({term}*J*beta)")
        calculated_z_str = " + ".join(calculated_z_parts)

        reason = (
            f"The answer is incorrect. The calculated degeneracies of the energy levels do not match the answer.\n"
            f"The system has 2^3 = 8 states. The energy is E = -J*(S1*S2 + S1*S3 + S2*S3).\n"
            f"The calculated energy levels and their degeneracies are:\n"
            f" - For states (+1,+1,+1) and (-1,-1,-1), the term (S1*S2+...) is 3. Energy = -3J. Degeneracy = 2.\n"
            f" - For the other 6 states, the term (S1*S2+...) is -1. Energy = J. Degeneracy = 6.\n"
            f"This leads to a partition function Z = 2*exp(-beta*(-3J)) + 6*exp(-beta*J) = 2*exp(3*J*beta) + 6*exp(-J*beta).\n"
            f"The calculated degeneracies are: {dict(degeneracy)}.\n"
            f"The degeneracies implied by answer D are: {expected_degeneracy}.\n"
            f"The calculated result ({calculated_z_str}) does not match the provided answer."
        )
        return reason

# Run the check
result = check_correctness()
print(result)