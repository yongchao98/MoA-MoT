import math
import re

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the quantum mechanics question.

    The logic is as follows:
    1.  Define the physical constraints (electric dipole selection rules).
    2.  Define a function to check if a transition between two states is allowed.
    3.  Define a function to calculate the theoretical probability of a specific decay path.
        - This is based on the spherical symmetry of the initial |n,0,0> state, which leads to
          equally probable decay channels to the intermediate |n',1,m'> states.
    4.  Parse each multiple-choice option to extract its decay path and given probability.
    5.  For each option, verify if its path is valid and its probability is correct.
    6.  Identify the single correct option based on this physical analysis.
    7.  Compare the identified correct option with the LLM's provided answer.
    8.  Return "Correct" if they match, or a detailed reason if they don't.
    """

    # The final answer provided by the LLM analysis.
    llm_answer = "A"

    # The options as presented in the question.
    # Note: The LLM answers sometimes re-label the options. We will check against the original question's options.
    question_options_text = {
        "A": r"|3,0,0\rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle and \frac{1}{3}",
        "B": r"|3,0,0\rangle\rightarrow|2,1,-1\rangle\rightarrow|1,0,0\rangle and \frac{1}{4}",
        "C": r"\rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle  and \frac{2}{3}",
        "D": r"|3,0,0\rangle\rightarrow|2,1,1\rangle\rightarrow|1,0,0\rangle and \frac{1}{4}"
    }

    def parse_state(state_str):
        """Parses a state string like '|n,l,m>' into a tuple (n, l, m)."""
        try:
            return tuple(map(int, state_str.strip('|⟩').split(',')))
        except (ValueError, AttributeError):
            return None

    def parse_option(option_text):
        """Parses the full text of an option to extract path and probability."""
        states_str = re.findall(r'\|[0-9,-]+\rangle', option_text)
        
        # Handle malformed option C
        if len(states_str) == 2 and "rightarrow|2,1,0" in option_text:
             path = [None, parse_state(states_str[0]), parse_state(states_str[1])]
        elif len(states_str) == 3:
            path = [parse_state(s) for s in states_str]
        else: # Malformed path
            path = [None, None, None]

        prob_match = re.search(r'\\frac\{(\d+)\}\{(\d+)\}', option_text)
        if prob_match:
            prob = int(prob_match.group(1)) / int(prob_match.group(2))
        else:
            prob = None
            
        return {"path": path, "prob": prob}

    def is_transition_allowed(initial_state, final_state):
        """Checks if a transition is allowed by E1 selection rules."""
        if initial_state is None or final_state is None:
            return False
        n_i, l_i, m_i = initial_state
        n_f, l_f, m_f = final_state

        # Rule 1: Principal quantum number must decrease for emission
        delta_n_ok = n_f < n_i
        # Rule 2: Orbital angular momentum must change by +/- 1
        delta_l_ok = abs(l_f - l_i) == 1
        # Rule 3: Magnetic quantum number must change by 0 or +/- 1
        delta_m_ok = abs(m_f - m_i) <= 1

        return delta_n_ok and delta_l_ok and delta_m_ok

    def get_correct_probability(initial_state):
        """Calculates the theoretical probability for a specific decay path from a symmetric state."""
        # This logic applies because the initial state |3,0,0> is spherically symmetric (l=0).
        if initial_state[1] != 0:
            return None # More complex calculation needed for non-symmetric states.

        # The intermediate state must have n=2 and l=1.
        # The possible m values are -1, 0, 1.
        num_possible_channels = 3
        
        # The probability of the first step is 1 / (number of channels).
        prob_step1 = 1 / num_possible_channels
        
        # The probability of the second step (|2,1,m'> -> |1,0,0>) is 1, as it's the only
        # allowed dipole decay to a lower energy level.
        prob_step2 = 1
        
        return prob_step1 * prob_step2

    # --- Main Verification Logic ---
    
    true_correct_options = []
    for key, text in question_options_text.items():
        data = parse_option(text)
        path = data["path"]
        given_prob = data["prob"]

        # Check 1: Path validity
        if not (path and len(path) == 3 and all(isinstance(s, tuple) or s is None for s in path)):
            continue # Skip malformed options
        
        path_is_valid = is_transition_allowed(path[0], path[1]) and is_transition_allowed(path[1], path[2])
        
        if not path_is_valid:
            continue

        # Check 2: Probability correctness
        correct_prob = get_correct_probability(path[0])
        
        if given_prob is not None and correct_prob is not None and math.isclose(given_prob, correct_prob):
            true_correct_options.append(key)

    if len(true_correct_options) != 1:
        return f"Analysis Error: Found {len(true_correct_options)} correct options ({true_correct_options}) based on physics. Expected 1."

    true_correct_answer = true_correct_options[0]

    if llm_answer == true_correct_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is '{llm_answer}', but the physically correct answer is '{true_correct_answer}'.\n"
                f"Reasoning: The selection rules allow three decay paths via intermediate states |2,1,-1>, |2,1,0>, and |2,1,1>. "
                f"Because the initial state |3,0,0> is spherically symmetric, the probability of decaying into any one of these three channels is equal, i.e., 1/3. "
                f"Option {true_correct_answer} is the only one that lists a valid path with the correct probability of 1/3.")

# Execute the check and print the result
print(check_correctness())