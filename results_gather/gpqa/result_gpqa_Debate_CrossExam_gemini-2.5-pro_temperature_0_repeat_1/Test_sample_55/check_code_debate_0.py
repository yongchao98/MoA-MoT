import numpy as np

def check_answer_correctness():
    """
    Checks the correctness of the provided answer for the H-atom decay problem.
    The code verifies the transition path against dipole selection rules and
    re-calculates the probability based on physical principles.
    """

    # --- Define the problem based on the selected answer (D) ---
    # Path: |3,0,0> -> |2,1,0> -> |1,0,0>
    # Probability: 1/3
    initial_state = {'n': 3, 'l': 0, 'm': 0}
    intermediate_state = {'n': 2, 'l': 1, 'm': 0}
    final_state = {'n': 1, 'l': 0, 'm': 0}
    given_probability = 1/3

    # --- Constraint 1: Check if the transition path is valid ---
    # A transition is allowed if Δl = ±1 and Δm = 0, ±1, and Δn < 0 for decay.

    # Check first transition: |3,0,0> -> |2,1,0>
    delta_n1 = intermediate_state['n'] - initial_state['n']
    delta_l1 = intermediate_state['l'] - initial_state['l']
    delta_m1 = intermediate_state['m'] - initial_state['m']

    if not (delta_n1 < 0 and abs(delta_l1) == 1 and abs(delta_m1) <= 1):
        return (f"Incorrect. The first transition |3,0,0> -> |2,1,0> violates the dipole selection rules. "
                f"Found Δn={delta_n1}, Δl={delta_l1}, Δm={delta_m1}.")

    # Check second transition: |2,1,0> -> |1,0,0>
    delta_n2 = final_state['n'] - intermediate_state['n']
    delta_l2 = final_state['l'] - intermediate_state['l']
    delta_m2 = final_state['m'] - intermediate_state['m']

    if not (delta_n2 < 0 and abs(delta_l2) == 1 and abs(delta_m2) <= 1):
        return (f"Incorrect. The second transition |2,1,0> -> |1,0,0> violates the dipole selection rules. "
                f"Found Δn={delta_n2}, Δl={delta_l2}, Δm={delta_m2}.")

    # --- Constraint 2: Check if the probability is correct ---
    # The probability of the route is P(step1) * P(step2).

    # Calculate P(step1): Branching ratio for |3,0,0> -> |2,1,0>
    # From an l=0 state, transitions to the l'=1 manifold are equally probable for all allowed m' values.
    # Allowed m' values are m + { -1, 0, 1 } = 0 + { -1, 0, 1 } = { -1, 0, 1 }.
    num_possible_intermediate_paths = 3  # to |2,1,-1>, |2,1,0>, |2,1,1>
    prob_step1 = 1.0 / num_possible_intermediate_paths

    # Calculate P(step2): Probability of |2,1,0> -> |1,0,0>
    # From the n=2 shell, the only possible decay is to the n=1 ground state.
    # Therefore, the probability of this decay, given the atom is in the |2,1,0> state, is 1.
    prob_step2 = 1.0

    # Calculate the total probability for the specific route
    calculated_probability = prob_step1 * prob_step2

    # Compare calculated probability with the given probability
    if not np.isclose(calculated_probability, given_probability):
        return (f"Incorrect. The transition path is valid, but the probability is wrong. "
                f"The initial state |3,0,0> (l=0) can decay to {num_possible_intermediate_paths} equally likely intermediate states in the n=2, l=1 manifold. "
                f"The probability of the first step is therefore 1/{num_possible_intermediate_paths}. "
                f"The probability of the second step is 1. "
                f"The total calculated probability is {calculated_probability:.4f}, but the answer provides {given_probability:.4f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_answer_correctness()
print(result)