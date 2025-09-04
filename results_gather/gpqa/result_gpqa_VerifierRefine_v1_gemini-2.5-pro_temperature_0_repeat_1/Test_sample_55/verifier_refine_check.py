import re
from fractions import Fraction

def check_correctness_of_llm_answer():
    """
    Checks the correctness of the LLM's answer for the H-atom decay problem.

    The logic is based on:
    1. Verifying the transition routes using electric dipole selection rules.
    2. Calculating the branching ratio for the specified cascade decay.
    """

    # Define a helper function to parse quantum state strings like '|n,l,m>'
    def parse_state(state_str):
        match = re.search(r'\|(\d+),(\d+),(-?\d+)', state_str)
        if match:
            return tuple(map(int, match.groups()))
        return None

    # Define a helper function to check E1 selection rules
    def check_dipole_selection_rules(initial_state, final_state):
        n_i, l_i, m_i = initial_state
        n_f, l_f, m_f = final_state
        
        # For decay, n must decrease
        if n_f >= n_i:
            return False, f"Transition fails: n must decrease (n_f={n_f} is not < n_i={n_i})."
        
        # Delta l must be +/- 1
        if abs(l_f - l_i) != 1:
            return False, f"Transition fails: |Δl| must be 1 (is {abs(l_f - l_i)})."
            
        # Delta m must be 0 or +/- 1
        if abs(m_f - m_i) > 1:
            return False, f"Transition fails: |Δm| must be <= 1 (is {abs(m_f - m_i)})."
            
        return True, "Allowed"

    # --- Physics Analysis ---

    # Step 1: Probability of the first transition: |3,0,0> -> |2,1,m'>
    # The initial state |3,0,0> is spherically symmetric (l=0).
    # The possible final states |2,1,m'> for m' in {-1, 0, 1} are degenerate.
    # Due to symmetry, the decay probability is equal for all three m' channels.
    # Assuming decay to n=2 is the only significant channel from n=3, l=0.
    num_channels_step1 = 3
    prob_step1 = Fraction(1, num_channels_step1)

    # Step 2: Probability of the second transition: |2,1,m'> -> |1,0,0>
    # From a |2,1,m'> state, the only allowed dipole decay to a lower energy level
    # is to the n=1 shell. The n=1 shell only has l=0.
    # Therefore, the transition to |1,0,0> is the only possible decay route.
    # The probability of this step, given the atom is in the intermediate state, is 1.
    prob_step2 = Fraction(1, 1)

    # Total probability for any valid cascade route
    correct_probability = prob_step1 * prob_step2

    # --- Verification of Options ---
    options = {
        'A': {'route': ('|3,0,0>', '|2,1,-1>', '|1,0,0>'), 'prob': Fraction(1, 4)},
        'B': {'route': ('|3,0,0>', '|2,1,0>', '|1,0,0>'), 'prob': Fraction(2, 3)},
        'C': {'route': ('|3,0,0>', '|2,1,1>', '|1,0,0>'), 'prob': Fraction(1, 4)},
        'D': {'route': ('|3,0,0>', '|2,1,0>', '|1,0,0>'), 'prob': Fraction(1, 3)}
    }
    
    llm_choice = 'D' # Based on the provided response " <<<D>>> "
    chosen_option = options[llm_choice]

    # 1. Check if the route in the chosen option is valid
    initial_s_str, inter_s_str, final_s_str = chosen_option['route']
    initial_s = parse_state(initial_s_str)
    inter_s = parse_state(inter_s_str)
    final_s = parse_state(final_s_str)

    is_allowed_1, reason_1 = check_dipole_selection_rules(initial_s, inter_s)
    if not is_allowed_1:
        return f"The chosen answer {llm_choice} is incorrect. The first transition {initial_s_str} -> {inter_s_str} is forbidden: {reason_1}"

    is_allowed_2, reason_2 = check_dipole_selection_rules(inter_s, final_s)
    if not is_allowed_2:
        return f"The chosen answer {llm_choice} is incorrect. The second transition {inter_s_str} -> {final_s_str} is forbidden: {reason_2}"

    # 2. Check if the probability is correct
    if chosen_option['prob'] != correct_probability:
        return (f"The chosen answer {llm_choice} is incorrect. The route is valid, but the probability is wrong. "
                f"The calculated probability is {float(correct_probability):.3f} ({correct_probability}), "
                f"but the option states {float(chosen_option['prob']):.3f} ({chosen_option['prob']}).")

    # 3. If all checks pass for the chosen answer, it is correct.
    return "Correct"

# Run the check
result = check_correctness_of_llm_answer()
print(result)