import numpy as np
import re

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by:
    1. Defining the problem's physical constraints (selection rules).
    2. Calculating the theoretical probability for the decay path.
    3. Parsing the options from the question.
    4. Evaluating each option against the constraints and theoretical probability.
    5. Comparing the LLM's choice with the derived correct answer.
    """

    # The LLM's final choice to be verified
    llm_choice = 'A'

    # Define the options based on the question.
    # Note: Option B is malformed in the prompt. We assume it refers to the same
    # path as A but with a different probability, as this is the most logical interpretation.
    options = {
        'A': {'path': [[3, 0, 0], [2, 1, 0], [1, 0, 0]], 'prob': 1/3},
        'B': {'path': [[3, 0, 0], [2, 1, 0], [1, 0, 0]], 'prob': 2/3},
        'C': {'path': [[3, 0, 0], [2, 1, -1], [1, 0, 0]], 'prob': 1/4},
        'D': {'path': [[3, 0, 0], [2, 1, 1], [1, 0, 0]], 'prob': 1/4}
    }

    def is_dipole_transition_allowed(state1, state2):
        """Checks if a transition between two states obeys electric dipole selection rules."""
        n1, l1, m1 = state1
        n2, l2, m2 = state2
        
        delta_n = n2 - n1
        delta_l = l2 - l1
        delta_m = m2 - m1
        
        # For decay, n must decrease
        if delta_n >= 0:
            return False, f"Transition {state1} -> {state2} is not a decay (Δn={delta_n} is not negative)."
        # Δl must be ±1
        if abs(delta_l) != 1:
            return False, f"Transition {state1} -> {state2} violates Δl=±1 rule (Δl={delta_l})."
        # Δm must be 0, ±1
        if abs(delta_m) > 1:
            return False, f"Transition {state1} -> {state2} violates Δm=0,±1 rule (Δm={delta_m})."
        
        return True, ""

    # The initial state |3,0,0> is spherically symmetric (l=0).
    # It can decay to the n=2, l=1 manifold. The possible m' states are -1, 0, 1.
    # There are 3 equally probable intermediate states.
    # The probability of the first step to any single intermediate state is 1/3.
    # The second step |2,1,m'> -> |1,0,0> is the only possible decay, so its probability is 1.
    # Therefore, the total probability for any single valid cascade path is 1/3.
    theoretical_prob = 1/3

    correct_options = []
    error_messages = {}

    for key, option in options.items():
        path = option['path']
        prob = option['prob']
        
        # Check first transition
        is_allowed1, reason1 = is_dipole_transition_allowed(path[0], path[1])
        if not is_allowed1:
            error_messages[key] = f"Option {key} has an invalid path. {reason1}"
            continue
            
        # Check second transition
        is_allowed2, reason2 = is_dipole_transition_allowed(path[1], path[2])
        if not is_allowed2:
            error_messages[key] = f"Option {key} has an invalid path. {reason2}"
            continue
            
        # Check probability
        if not np.isclose(prob, theoretical_prob):
            error_messages[key] = f"Option {key} has a valid path but incorrect probability. Expected {theoretical_prob:.4f}, but got {prob:.4f}."
            continue
            
        # If all checks pass, this is a correct option
        correct_options.append(key)

    if not correct_options:
        return "Error: The code found no correct options among the choices. The question or options may be flawed."
    
    if len(correct_options) > 1:
        return f"Error: The code found multiple correct options: {correct_options}. The question is ambiguous."

    # Check if the LLM's answer matches the one found by the code
    correct_key = correct_options[0]
    if llm_choice == correct_key:
        return "Correct"
    else:
        reason = error_messages.get(llm_choice, f"The correct option is {correct_key}, not {llm_choice}.")
        return f"Incorrect. The provided answer '{llm_choice}' is wrong. Reason: {reason}"

# Execute the check
result = check_correctness()
print(result)