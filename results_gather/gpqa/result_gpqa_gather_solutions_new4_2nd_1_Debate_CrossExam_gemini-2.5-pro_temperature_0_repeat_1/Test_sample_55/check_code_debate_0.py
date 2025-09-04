import re
from fractions import Fraction

def parse_state(state_str):
    """
    Parses a quantum state string like "|n,l,m>" into a tuple of integers (n, l, m).
    Returns None if the string is malformed.
    """
    if not isinstance(state_str, str):
        return None
    # Remove bra-ket notation and whitespace
    cleaned_str = state_str.strip().replace('|', '').replace('>', '')
    if not cleaned_str:
        return None
    try:
        parts = [int(x) for x in cleaned_str.split(',')]
        if len(parts) == 3:
            return tuple(parts)
    except (ValueError, IndexError):
        return None
    return None

def parse_path(path_str):
    """
    Parses a transition path string like "state1 -> state2 -> state3"
    into a list of state tuples. Returns None if any part is malformed.
    """
    if not isinstance(path_str, str):
        return None
    parts = path_str.split('->')
    states = [parse_state(p) for p in parts]
    # Check if any state failed to parse
    if any(s is None for s in states):
        return None
    return states

def is_valid_dipole_transition(initial_state, final_state):
    """
    Checks if a transition between two states is allowed by the electric
    dipole (E1) selection rules for a decay process.
    """
    n_i, l_i, m_i = initial_state
    n_f, l_f, m_f = final_state
    
    # Rule 1: Principal quantum number must decrease for decay
    rule_n = (n_f < n_i)
    
    # Rule 2: Orbital angular momentum quantum number must change by +/- 1
    rule_l = (abs(l_f - l_i) == 1)
    
    # Rule 3: Magnetic quantum number must change by 0 or +/- 1
    rule_m = (abs(m_f - m_i) in [0, 1])
    
    return rule_n and rule_l and rule_m

def check_correctness(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given quantum mechanics question.

    Args:
        question_text (str): The text of the multiple-choice question.
        llm_answer_text (str): The text of the LLM's answer, including the final choice.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """
    # --- Step 1: Parse the question and the LLM's answer ---

    # Extract the final choice (e.g., 'A') from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<X>>> in the provided response."
    llm_choice = match.group(1)

    # Parse the options from the question text
    # This is a simplified parser for the specific format in the question
    options_raw = re.findall(r'\n([A-D])\)\s*(.*)', question_text)
    if not options_raw:
        return "Error: Could not parse the options from the question text."

    options = {}
    for key, content in options_raw:
        parts = content.split(' and ')
        if len(parts) == 2:
            path_str = parts[0].replace('\\rangle', '>') # Normalize latex
            prob_str = parts[1]
            options[key] = {'path_str': path_str, 'prob_str': prob_str}
        else:
            options[key] = {'path_str': content, 'prob_str': ''} # Handle malformed options

    # --- Step 2: Define the correct physics principles ---

    # For a decay from a spherically symmetric state (|n,0,0>), the probability
    # is distributed equally among the possible degenerate intermediate states.
    # First transition: |3,0,0> -> |2,1,m'>. There are 3 possible m' values (-1, 0, 1).
    # So, the probability for any one path is 1/3.
    # Second transition: |2,1,m'> -> |1,0,0>. This is the only possible decay, so P=1.
    # Total probability for any valid route is (1/3) * 1 = 1/3.
    correct_probability = Fraction(1, 3)

    # --- Step 3: Evaluate all options to find the correct one(s) ---
    
    true_correct_options = []
    for key, data in options.items():
        path = parse_path(data['path_str'])
        
        # Check 1: Is the path syntactically valid?
        if path is None:
            continue # Skip malformed options like B

        # Check 2: Does the path follow selection rules?
        path_is_valid = True
        if len(path) != 3: # Expecting a 2-step transition (3 states)
            path_is_valid = False
        else:
            for i in range(len(path) - 1):
                if not is_valid_dipole_transition(path[i], path[i+1]):
                    path_is_valid = False
                    break
        
        if not path_is_valid:
            continue

        # Check 3: Is the probability correct?
        try:
            claimed_probability = Fraction(data['prob_str'])
            if claimed_probability == correct_probability:
                true_correct_options.append(key)
        except (ValueError, ZeroDivisionError):
            continue

    # --- Step 4: Compare LLM's choice with the ground truth ---

    if not true_correct_options:
        return "Error: Analysis shows that none of the provided options are correct."

    if len(true_correct_options) > 1:
        # This case is for questions with multiple correct answers
        if llm_choice in true_correct_options:
            return "Correct"
        else:
            return f"Incorrect. The LLM chose {llm_choice}, but the correct answers are {true_correct_options}."

    # Standard case: one correct answer
    correct_answer_key = true_correct_options[0]
    if llm_choice == correct_answer_key:
        return "Correct"
    else:
        # Provide a detailed reason for the error
        chosen_option_data = options.get(llm_choice)
        if not chosen_option_data:
             return f"Incorrect. The LLM chose an invalid option key '{llm_choice}'. The correct answer is {correct_answer_key}."

        chosen_path = parse_path(chosen_option_data['path_str'])
        if chosen_path is None:
            return f"Incorrect. The LLM chose option {llm_choice}, which has a malformed path. The correct answer is {correct_answer_key}."

        path_is_valid = True
        if len(chosen_path) != 3:
            path_is_valid = False
        else:
            for i in range(len(chosen_path) - 1):
                if not is_valid_dipole_transition(chosen_path[i], chosen_path[i+1]):
                    path_is_valid = False
                    break
        
        if not path_is_valid:
            return f"Incorrect. The LLM chose option {llm_choice}, whose path violates the E1 selection rules. The correct answer is {correct_answer_key}."
        
        try:
            claimed_prob = Fraction(chosen_option_data['prob_str'])
            if claimed_prob != correct_probability:
                return f"Incorrect. The LLM chose option {llm_choice}. While the path is valid, its stated probability of {claimed_prob} is wrong. The correct probability is {correct_probability}. The correct option is {correct_answer_key}."
        except (ValueError, ZeroDivisionError):
             return f"Incorrect. The LLM chose option {llm_choice}, which has an unparsable probability. The correct answer is {correct_answer_key}."

        # This case should not be reached but is a fallback
        return f"Incorrect. The LLM chose {llm_choice}, but the correct answer is {correct_answer_key}."


# --- Provided inputs for the check ---
question = """
The |3,0,0\\rangle state in the standard notation |n,l,m\\rangle of the H -atom in the non-relativistic theory decays to the state |1,0,0\\rangle via two dipole transition. The transition route and the corresponding probability are (use latex),

A) |3,0,0\\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle and \\frac{1}{3}
B) \\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle  and \\frac{2}{3}
C) |3,0,0\\rangle\\rightarrow|2,1,-1\\rangle\\rightarrow|1,0,0\\rangle and \\frac{1}{4}
D) |3,0,0\\rangle\\rightarrow|2,1,1\\rangle\\rightarrow|1,0,0\\rangle
and \\frac{1}{4}
"""

llm_answer = """
Here is a step-by-step analysis of the question, incorporating the careful points derived from the candidate answers.

### Step 1: Understand the Physical Process and Governing Rules
The problem describes a two-step decay of a hydrogen atom from the initial state $|3,0,0\\rangle$ to the final state $|1,0,0\\rangle$. Each step is an electric dipole (E1) transition. These transitions are governed by a set of selection rules for the quantum numbers $|n, l, m\\rangle$:
*   **Rule for l (Orbital Angular Momentum):** The change in the orbital angular momentum quantum number must be $\\Delta l = \\pm 1$.
*   **Rule for m (Magnetic Quantum Number):** The change in the magnetic quantum number must be $\\Delta m = 0, \\pm 1$.
*   **Rule for n (Principal Quantum Number):** For a decay process (spontaneous emission), the energy must decrease, so the principal quantum number must decrease ($n_{final} < n_{initial}$).

### Step 2: Analyze the First Transition Step
The atom starts in the initial state $|i\\rangle = |3,0,0\\rangle$ and decays to an intermediate state $|k\\rangle = |n',l',m'\\rangle$.
*   **Initial State:** $n=3, l=0, m=0$.
*   **Applying the $\\Delta l$ rule:** Since the initial state has $l=0$, the intermediate state must have $l' = 1$.
*   **Applying the $\\Delta m$ rule:** Since the initial state has $m=0$, the intermediate state can have $m' = 0, +1, \\text{ or } -1$.
*   **Applying the $\\Delta n$ rule:** The principal quantum number must decrease. The options suggest an intermediate state with $n'=2$, which is a valid decay path ($2 < 3$).
*   **Conclusion:** The possible intermediate states are $|2,1,-1\\rangle$, $|2,1,0\\rangle$, and $|2,1,1\\rangle$.

### Step 3: Analyze the Second Transition Step
The atom then decays from an intermediate state $|k\\rangle = |2,1,m'\\rangle$ to the final state $|f\\rangle = |1,0,0\\rangle$.
*   **Applying the $\\Delta l$ rule:** The transition is from $l'=1$ to $l_f=0$. The change is $\\Delta l = 0 - 1 = -1$, which is allowed.
*   **Applying the $\\Delta m$ rule:** The transition is from $m'$ (which can be -1, 0, or 1) to $m_f=0$. In all three cases, the change $\\Delta m = 0 - m'$ is either $0$ or $\\pm 1$, which is allowed.
*   **Conclusion:** All three potential pathways, via any of the $|2,1,m'\\rangle$ states, are physically valid according to the selection rules.

### Step 4: Calculate the Probability (Branching Ratio)
This is the most critical part of the analysis. The probability of a specific cascade route is the product of the probabilities of each step.

*   **Careful Point 1: Spherical Symmetry of the Initial State.** The initial state $|3,0,0\\rangle$ is an *s*-orbital, which is spherically symmetric. This means it has no preferred orientation in space. In the absence of external fields, the decay must be isotropic.
*   **Careful Point 2: Equal Probability for Degenerate States.** Because the decay is isotropic, the atom is equally likely to transition into any of the three possible intermediate states: $|2,1,-1\\rangle$, $|2,1,0\\rangle$, and $|2,1,1\\rangle$. These states are degenerate in energy. Since there are three equally likely outcomes for the first decay step, the probability (or branching ratio) for any *single specific path* is **1/3**.
*   **Careful Point 3: Certainty of the Second Step.** Once the atom is in any of the $|2,1,m'\\rangle$ states, the only possible E1 decay to a lower energy level is to the $|1,0,0\\rangle$ ground state. Therefore, the probability of this second step occurring is **1**.

The total probability for any of the three specific routes is the product of the probabilities of the two steps:
$P(\\text{route}) = P(\\text{step 1}) \\times P(\\text{step 2}) = \\frac{1}{3} \\times 1 = \\frac{1}{3}$.

### Step 5: Evaluate the Given Options
We have determined that any valid route has a probability of 1/3. Let's check the original options:

*   **A) $|3,0,0\\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle$ and $\\frac{1}{3}$**: The transition route is valid, and the calculated probability is correct.
*   **B) $\\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle$ and $\\frac{2}{3}$**: The route is malformed, and the probability is incorrect.
*   **C) $|3,0,0\\rangle\\rightarrow|2,1,-1\\rangle\\rightarrow|1,0,0\\rangle$ and $\\frac{1}{4}$**: The transition route is valid, but the probability is incorrect.
*   **D) $|3,0,0\\rangle\\rightarrow|2,1,1\\rangle\\rightarrow|1,0,0\\rangle$ and $\\frac{1}{4}$**: The transition route is valid, but the probability is incorrect.

Based on the analysis, only option A correctly identifies a valid transition route and its corresponding probability.

<<<A>>>
"""

# Run the check
result = check_correctness(question, llm_answer)
print(result)