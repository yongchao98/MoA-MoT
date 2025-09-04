import re
from fractions import Fraction

def check_correctness():
    """
    Checks the correctness of the LLM's answer based on the physics of atomic transitions.
    """

    # --- Define the problem constraints and physics principles ---

    # The initial and final states of the overall process
    question_initial_state = {'n': 3, 'l': 0, 'm': 0}
    question_final_state = {'n': 1, 'l': 0, 'm': 0}

    # Physics Principle 1: Dipole Selection Rules
    def is_dipole_transition_allowed(initial_state, final_state):
        """Checks if a transition is allowed by E1 selection rules for decay."""
        if final_state['n'] >= initial_state['n']:
            return False, f"Energy must decrease (n_final={final_state['n']} >= n_initial={initial_state['n']})."
        
        delta_l = final_state['l'] - initial_state['l']
        if abs(delta_l) != 1:
            return False, f"Δl must be ±1, but it was {delta_l}."
            
        delta_m = final_state['m'] - initial_state['m']
        if abs(delta_m) > 1:
            return False, f"Δm must be 0 or ±1, but it was {delta_m}."
            
        return True, "Allowed"

    # Physics Principle 2: Branching Ratios for s-state decay
    # The initial state |3,0,0> is an s-state (l=0), which is spherically symmetric.
    # It can decay to the n=2, l=1 manifold.
    # The possible intermediate m' states are -1, 0, 1.
    # Due to symmetry, the probability of decaying to each of these 3 states is equal.
    num_intermediate_paths = 3
    prob_step1 = Fraction(1, num_intermediate_paths)
    
    # The intermediate state |2,1,m'> can only decay to |1,0,0> via E1 transition.
    prob_step2 = 1
    
    correct_probability = prob_step1 * prob_step2

    # --- Parse and evaluate the LLM's chosen answer ---

    llm_answer_choice = 'A'
    
    options = {
        'A': "|3,0,0>→|2,1,0>→|1,0,0> and \\frac{1}{3}",
        'B': "|3,0,0>→|2,1,-1>→|1,0,0> and \\frac{1}{4}",
        'C': ">→|2,1,0>→|1,0,0> and \\frac{2}{3}",
        'D': "|3,0,0>→|2,1,1>→|1,0,0> and \\frac{1}{4}"
    }

    def parse_state(s):
        match = re.match(r'\|(\d+),(\d+),(-?\d+)\rangle', s)
        if not match: return None
        return {'n': int(match.group(1)), 'l': int(match.group(2)), 'm': int(match.group(3))}

    def parse_option(option_str):
        pattern = r"(\|.+?⟩)\s*→\s*(\|.+?⟩)\s*→\s*(\|.+?⟩)\s*and\s*\\frac\{(\d+)\}\{(\d+)\}"
        match = re.search(pattern, option_str)
        if not match: return {'valid_format': False}
        
        s1_str, s2_str, s3_str, num, den = match.groups()
        return {
            'valid_format': True,
            'initial': parse_state(s1_str),
            'intermediate': parse_state(s2_str),
            'final': parse_state(s3_str),
            'probability': Fraction(int(num), int(den))
        }

    chosen_option_str = options.get(llm_answer_choice)
    if not chosen_option_str:
        return f"Invalid answer choice '{llm_answer_choice}'."

    parsed_option = parse_option(chosen_option_str)

    # 1. Check if the option format is valid
    if not parsed_option['valid_format']:
        return f"The chosen option '{llm_answer_choice}' is malformed: '{chosen_option_str}'"

    # 2. Check if the states in the option match the question
    if parsed_option['initial'] != question_initial_state:
        return f"Constraint violated: The initial state in option {llm_answer_choice} {parsed_option['initial']} does not match the question's initial state {question_initial_state}."
    if parsed_option['final'] != question_final_state:
        return f"Constraint violated: The final state in option {llm_answer_choice} {parsed_option['final']} does not match the question's final state {question_final_state}."

    # 3. Check if the transition path is allowed by selection rules
    allowed1, reason1 = is_dipole_transition_allowed(parsed_option['initial'], parsed_option['intermediate'])
    if not allowed1:
        return f"Constraint violated: The first transition in option {llm_answer_choice} is not allowed. Reason: {reason1}"
    
    allowed2, reason2 = is_dipole_transition_allowed(parsed_option['intermediate'], parsed_option['final'])
    if not allowed2:
        return f"Constraint violated: The second transition in option {llm_answer_choice} is not allowed. Reason: {reason2}"

    # 4. Check if the probability is correct
    if parsed_option['probability'] != correct_probability:
        return f"Constraint violated: The probability in option {llm_answer_choice} is {parsed_option['probability']}, but the correct probability for any valid path is {correct_probability}."

    return "Correct"

# Run the check
result = check_correctness()
print(result)