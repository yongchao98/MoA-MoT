import re
from fractions import Fraction

def check_correctness():
    """
    Checks the correctness of the LLM's answer about a two-step dipole transition.
    
    The function verifies two main aspects:
    1. The validity of the transition route based on electric dipole selection rules.
    2. The correctness of the associated probability (branching ratio).
    """

    # --- Problem Definition & LLM's Answer ---
    # The final answer provided by the LLM is 'D', which corresponds to:
    # Route: |3,0,0> -> |2,1,0> -> |1,0,0>
    # Probability: 1/3
    
    llm_answer_route_str = "|3,0,0>->|2,1,0>->|1,0,0>"
    llm_answer_prob_str = "1/3"

    # --- Helper Functions ---
    def parse_state(state_str):
        """Parses a state string like '|n,l,m>' into a tuple (n, l, m)."""
        match = re.search(r'\|\s*(\d+)\s*,\s*(\d+)\s*,\s*([+-]?\d+)\s*>', state_str)
        if match:
            n, l, m = map(int, match.groups())
            return (n, l, m)
        return None

    def check_dipole_transition(state1, state2):
        """Checks if an electric dipole transition is allowed between state1 and state2."""
        n1, l1, m1 = state1
        n2, l2, m2 = state2

        # Rule 1: Principal quantum number must decrease for decay.
        if n2 >= n1:
            return False, f"Invalid transition: n must decrease for decay (from n={n1} to n={n2})."

        # Rule 2: Change in orbital angular momentum quantum number must be ±1.
        if abs(l2 - l1) != 1:
            return False, f"Invalid transition: Δl must be ±1, but it was {l2 - l1} (from l={l1} to l={l2})."

        # Rule 3: Change in magnetic quantum number must be 0 or ±1.
        if abs(m2 - m1) > 1:
            return False, f"Invalid transition: Δm must be 0 or ±1, but it was {m2 - m1} (from m={m1} to m={m2})."

        return True, "Allowed"

    # --- Verification Logic ---

    # 1. Parse the states from the answer's route string
    try:
        initial_state_str, intermediate_state_str, final_state_str = llm_answer_route_str.split("->")
        initial_state = parse_state(initial_state_str)
        intermediate_state = parse_state(intermediate_state_str)
        final_state = parse_state(final_state_str)
    except (ValueError, AttributeError):
        return "Error: Could not parse the transition route string in the answer."

    # 2. Check if the transition route is valid according to selection rules
    is_allowed1, reason1 = check_dipole_transition(initial_state, intermediate_state)
    if not is_allowed1:
        return f"The answer is incorrect. The first transition in the route {llm_answer_route_str} is forbidden. Reason: {reason1}"
    
    is_allowed2, reason2 = check_dipole_transition(intermediate_state, final_state)
    if not is_allowed2:
        return f"The answer is incorrect. The second transition in the route {llm_answer_route_str} is forbidden. Reason: {reason2}"

    # 3. Calculate the correct probability for this route
    # The initial state |3,0,0> is spherically symmetric (l=0). It decays with equal
    # probability into the three possible intermediate m' states (-1, 0, 1) of the l'=1 level.
    # Number of possible intermediate m' states for l'=1 is 2*l' + 1 = 3.
    prob_step1 = Fraction(1, 3)
    
    # The intermediate state |2,1,m'> can only decay to |1,0,0> via a dipole transition.
    # Therefore, the probability of the second step is 1.
    prob_step2 = 1
    
    correct_probability = prob_step1 * prob_step2
    
    # 4. Compare with the probability given in the answer
    try:
        answer_probability = Fraction(llm_answer_prob_str)
    except ValueError:
        return f"Error: Could not parse the probability '{llm_answer_prob_str}' from the answer."
    
    if correct_probability == answer_probability:
        return "Correct"
    else:
        return (f"The answer is incorrect. The transition route is valid, but the probability is wrong. "
                f"The answer claims the probability is {answer_probability}, but the correct probability "
                f"for any single path from a spherically symmetric state is {correct_probability}.")

# Execute the check and print the result
result = check_correctness()
print(result)