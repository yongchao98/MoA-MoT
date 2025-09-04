import re
import math

def parse_state(state_str):
    """Parses a string like '|n,l,m>' into a tuple of integers (n, l, m)."""
    if not state_str or not state_str.startswith('|'):
        return None
    try:
        # Remove kets, bras, and spaces, then split by comma
        numbers = state_str.strip('|> \rangle').split(',')
        return tuple(map(int, numbers))
    except (ValueError, IndexError):
        return None

def is_transition_allowed(state1, state2):
    """Checks if an electric dipole transition between two states is allowed."""
    if not state1 or not state2:
        return False
    n1, l1, m1 = state1
    n2, l2, m2 = state2
    
    # For decay, n must decrease
    if n2 >= n1:
        return False
        
    # Check selection rules for l and m
    delta_l = abs(l2 - l1)
    delta_m = abs(m2 - m1)
    
    if delta_l == 1 and delta_m <= 1:
        return True
        
    return False

def check_answer():
    """
    Checks the correctness of the provided answer by applying physics principles.
    """
    # The options as described in the prompt
    options = {
        'A': {'path_str': ['|3,0,0>', '|2,1,1>', '|1,0,0>'], 'prob': 1/4},
        'B': {'path_str': ['|3,0,0>', '|2,1,0>', '|1,0,0>'], 'prob': 1/3},
        'C': {'path_str': ['>', '|2,1,0>', '|1,0,0>'], 'prob': 2/3},
        'D': {'path_str': ['|3,0,0>', '|2,1,-1>', '|1,0,0>'], 'prob': 1/4}
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer = 'B'

    # --- Theoretical Calculation ---
    initial_state = (3, 0, 0)
    final_state = (1, 0, 0)
    
    # Find all possible intermediate states for the first transition
    possible_intermediate_states = []
    # Assuming intermediate n=2 as it's the only one that allows a cascade to n=1
    n_intermediate = 2
    # From l=0, must go to l=1
    l_intermediate = 1
    # From m=0, can go to m=-1, 0, 1
    for m_intermediate in [-1, 0, 1]:
        state = (n_intermediate, l_intermediate, m_intermediate)
        if is_transition_allowed(initial_state, state):
            possible_intermediate_states.append(state)
            
    # The number of equally likely paths determines the probability
    num_paths = len(possible_intermediate_states)
    if num_paths == 0:
        return "Error: No valid intermediate states found."
        
    # The probability for any single path is 1 / (number of paths)
    # because the initial state is spherically symmetric (l=0).
    theoretical_prob = 1.0 / num_paths

    # --- Verification ---
    correct_options = []
    for option_letter, data in options.items():
        path_is_valid = True
        
        # Parse path from strings to tuples
        path = [parse_state(s) for s in data['path_str']]
        
        # Check if path is well-formed and matches start/end states
        if None in path or path[0] != initial_state or path[-1] != final_state:
            path_is_valid = False
        else:
            # Check each transition in the path
            if not is_transition_allowed(path[0], path[1]) or not is_transition_allowed(path[1], path[2]):
                path_is_valid = False

        # Check if the probability matches the theoretical one
        prob_is_correct = math.isclose(data['prob'], theoretical_prob)
        
        if path_is_valid and prob_is_correct:
            correct_options.append(option_letter)

    # --- Final Judgement ---
    if llm_answer in correct_options:
        if len(correct_options) == 1:
            return "Correct"
        else:
            return f"Incorrect. The provided answer '{llm_answer}' is one of multiple correct options found: {correct_options}. The question may be ambiguous."
    else:
        # Provide a reason for the failure
        if llm_answer not in options:
            return f"Incorrect. The answer '{llm_answer}' is not a valid option."
        
        data = options[llm_answer]
        path = [parse_state(s) for s in data['path_str']]
        path_is_valid = True
        if None in path or not is_transition_allowed(path[0], path[1]) or not is_transition_allowed(path[1], path[2]):
            path_is_valid = False
        
        if not path_is_valid:
            return f"Incorrect. The path in option {llm_answer} is not physically allowed by the E1 selection rules."
        
        prob_is_correct = math.isclose(data['prob'], theoretical_prob)
        if not prob_is_correct:
            return f"Incorrect. The path in option {llm_answer} is valid, but its probability is wrong. The option states {data['prob']:.2f}, but the correct probability is 1/{num_paths} (~{theoretical_prob:.2f})."
        
        return f"Incorrect. An unknown error occurred while checking answer {llm_answer}."

# Run the check and print the result
result = check_answer()
print(result)