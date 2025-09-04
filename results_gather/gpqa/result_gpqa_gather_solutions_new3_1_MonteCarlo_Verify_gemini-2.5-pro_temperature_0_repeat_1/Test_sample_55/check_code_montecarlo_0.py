import re
import math

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer to the physics question.

    The function verifies two things:
    1.  The transition route proposed in the chosen option is physically allowed
        according to the electric dipole selection rules.
    2.  The probability (branching ratio) given in the chosen option is correct.
    """

    # --- Helper functions for parsing ---

    def parse_state(state_str):
        """Parses a state string like '|n,l,m>' into a dictionary."""
        state_str = state_str.strip().replace('>', '').replace('<', '')
        if not state_str.startswith('|'):
            return None
        try:
            parts = state_str[1:].split(',')
            return {'n': int(parts[0]), 'l': int(parts[1]), 'm': int(parts[2])}
        except (ValueError, IndexError):
            return None

    def parse_route(route_str, default_initial_state=None):
        """Parses a route string into a list of state dictionaries."""
        route_str = route_str.replace(r'\rightarrow', '->')
        parts = route_str.split('->')
        states = [parse_state(p) for p in parts]
        
        # Handle malformed first state by using the default from the question
        if len(states) > 0 and states[0] is None and default_initial_state is not None:
            states[0] = default_initial_state
            
        if None in states:
            return None  # Indicates a parsing failure
        return states

    def parse_prob(prob_str):
        """Parses a probability string like '\\frac{1}{3}' into a float."""
        match = re.search(r'\\frac\{(\d+)\}\{(\d+)\}', prob_str)
        if match:
            return float(match.group(1)) / float(match.group(2))
        try:
            if '/' in prob_str:
                num, den = prob_str.split('/')
                return float(num) / float(den)
            return float(prob_str)
        except (ValueError, TypeError):
            return None

    # --- Physics Logic ---

    def check_dipole_transition(state1, state2):
        """Checks if a transition between two states is an allowed dipole transition."""
        if not all(k in state1 and k in state2 for k in ['n', 'l', 'm']):
             return False, "State information is incomplete."
        # Rule: n must decrease for decay
        if not state2['n'] < state1['n']:
            return False, f"Principal quantum number n must decrease (from {state1['n']} to {state2['n']})."
        
        # Rule: delta l = +/- 1
        delta_l = abs(state2['l'] - state1['l'])
        if delta_l != 1:
            return False, f"Change in orbital quantum number l must be +/- 1, but it was {state2['l'] - state1['l']}."

        # Rule: delta m = 0, +/- 1
        delta_m = abs(state2['m'] - state1['m'])
        if delta_m > 1:
            return False, f"Change in magnetic quantum number m must be 0 or +/- 1, but it was {state2['m'] - state1['m']}."
            
        return True, ""

    # --- Main Checking Logic ---

    # 1. Define the correct physics
    # For a transition from a spherically symmetric l=0 state, the branching ratio to each of the three
    # possible m' states (for l'=1) is 1/3. The subsequent decay from |2,1,m'> to |1,0,0> has a
    # probability of 1. Therefore, the total probability for any valid cascade is 1/3.
    correct_probability = 1/3.0
    default_initial_state = {'n': 3, 'l': 0, 'm': 0}

    # 2. Extract the chosen option from the LLM's answer
    chosen_option_match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not chosen_option_match:
        return "Incorrect: The final answer format is invalid. It should be <<<X>>> where X is A, B, C, or D."
    chosen_option_letter = chosen_option_match.group(1)

    # 3. Parse all options from the question text
    options = {}
    option_lines = re.findall(r'^[A-D]\)\s*(.*)', question, re.MULTILINE)
    
    if len(option_lines) != 4:
        return "Error: Could not parse the four options (A, B, C, D) from the question text."

    for i, line in enumerate(option_lines):
        letter = chr(ord('A') + i)
        parts = line.split(' and ')
        if len(parts) != 2:
            return f"Error: Could not parse option {letter}. Expected 'route and probability' format."
        
        route = parse_route(parts[0].strip(), default_initial_state=default_initial_state)
        prob = parse_prob(parts[1].strip())

        if route is None: return f"Error: Could not parse the route for option {letter}: '{parts[0].strip()}'"
        if prob is None: return f"Error: Could not parse the probability for option {letter}: '{parts[1].strip()}'"
            
        options[letter] = {'route': route, 'prob': prob}

    # 4. Check the chosen option
    chosen_option = options.get(chosen_option_letter)
    if not chosen_option:
        return f"Incorrect: The chosen option '{chosen_option_letter}' does not exist in the question."

    route = chosen_option['route']
    given_prob = chosen_option['prob']

    # Check 4a: Route validity
    if len(route) != 3:
        return f"Incorrect: The route for option {chosen_option_letter} does not have two steps (three states)."
    
    initial_state, intermediate_state, final_state = route[0], route[1], route[2]

    is_valid_1, reason_1 = check_dipole_transition(initial_state, intermediate_state)
    if not is_valid_1:
        return f"Incorrect: The chosen option {chosen_option_letter} has an invalid transition path. The first step from |{initial_state['n']},{initial_state['l']},{initial_state['m']}> to |{intermediate_state['n']},{intermediate_state['l']},{intermediate_state['m']}> is not allowed. Reason: {reason_1}"

    is_valid_2, reason_2 = check_dipole_transition(intermediate_state, final_state)
    if not is_valid_2:
        return f"Incorrect: The chosen option {chosen_option_letter} has an invalid transition path. The second step from |{intermediate_state['n']},{intermediate_state['l']},{intermediate_state['m']}> to |{final_state['n']},{final_state['l']},{final_state['m']}> is not allowed. Reason: {reason_2}"

    # Check 4b: Probability correctness
    if not math.isclose(given_prob, correct_probability, rel_tol=1e-5):
        return f"Incorrect: The transition route in option {chosen_option_letter} is valid, but the probability is wrong. The correct probability for this cascade is 1/3 (approx 0.333), but the option states {given_prob}."

    # 5. If all checks pass
    return "Correct"

# --- Execution ---
# We extract the question and the final answer from the provided text.
question_text = """
The |3,0,0\rangle state in the standard notation |n,l,m\rangle of the H -atom in the non-relativistic theory decays to the state |1,0,0\rangle via two dipole transition. The transition route and the corresponding probability are (use latex),

A) \rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle  and \frac{2}{3}
B) |3,0,0\rangle\rightarrow|2,1,-1\rangle\rightarrow|1,0,0\rangle and \frac{1}{4}
C) |3,0,0\rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle and \frac{1}{3}
D) |3,0,0\rangle\rightarrow|2,1,1\rangle\rightarrow|1,0,0\rangle
and \frac{1}{4}
"""

llm_output = """
Based on a step-by-step analysis of the problem and a synthesis of the provided answers, here is a detailed breakdown.

### Step-by-step Analysis

1.  **Identify the Physical Process:** The problem describes a two-step radiative decay of a hydrogen atom from the |3,0,0⟩ state to the |1,0,0⟩ state. Each step is an electric dipole (E1) transition.

2.  **Recall the Selection Rules:** For an electric dipole transition between an initial state |n,l,m⟩ and a final state |n',l',m'⟩, the following selection rules must be met:
    *   The change in the orbital angular momentum quantum number must be Δl = ±1.
    *   The change in the magnetic quantum number must be Δm = 0, ±1.
    *   For spontaneous emission (decay), the energy must decrease, so n' < n.

3.  **Analyze the First Transition:**
    *   **Initial State:** |3,0,0⟩ (n=3, l=0, m=0).
    *   **Intermediate State:** |n',l',m'⟩.
    *   Applying the Δl rule: Since l=0, the intermediate state must have l' = 1.
    *   Applying the Δm rule: Since m=0, the intermediate state can have m' = -1, 0, or +1.
    *   Applying the Δn rule: n' must be less than 3. The options suggest n'=2.
    *   **Conclusion:** The possible intermediate states are |2,1,-1⟩, |2,1,0⟩, and |2,1,1⟩.

4.  **Analyze the Second Transition:**
    *   **Intermediate State:** |2,1,m'⟩.
    *   **Final State:** |1,0,0⟩.
    *   Applying the Δl rule: The change is Δl = 0 - 1 = -1. This is allowed.
    *   Applying the Δm rule: The change is Δm = 0 - m'. For m' ∈ {-1, 0, 1}, the change Δm will be +1, 0, or -1, respectively. All are allowed.
    *   **Conclusion:** All three potential decay routes are physically possible according to the selection rules.

5.  **Calculate the Transition Probability:** The probability of a specific cascade route is the product of the branching ratios of each step.
    *   **Probability of the second step:** Once the atom is in any of the |2,1,m'⟩ states, the only possible E1 decay to a lower energy level is to the |1,0,0⟩ ground state. Therefore, the probability of this step is 1.
    *   **Probability of the first step:** This is the branching ratio for the decay from |3,0,0⟩ to a specific |2,1,m'⟩ state. The initial state |3,0,0⟩ is spherically symmetric (l=0). In the absence of external fields, there is no preferred direction in space. This means the decay must populate each of the three possible final magnetic sublevels (m' = -1, 0, 1) with equal probability. Since there are three equally likely outcomes, the probability for any single one is 1/3.
    *   **Total Probability:** The total probability for any of the three specific routes is P(route) = P(step 1) × P(step 2) = (1/3) × 1 = 1/3.

### Evaluation of Options

We are looking for an option that presents a valid route with a probability of 1/3.

*   A) |3,0,0⟩→|2,1,0⟩→|1,0,0⟩ and 2/3. The route is valid, but the probability is incorrect.
*   B) |3,0,0⟩→|2,1,-1⟩→|1,0,0⟩ and 1/4. The route is valid, but the probability is incorrect.
*   C) |3,0,0⟩→|2,1,0⟩→|1,0,0⟩ and 1/3. The route is valid, and the probability is correct.
*   D) |3,0,0⟩→|2,1,1⟩→|1,0,0⟩ and 1/4. The route is valid, but the probability is incorrect.

The only option that correctly identifies a valid transition route and its corresponding probability is C.

<<<C>>>
"""

result = check_correctness(question_text, llm_output)
print(result)