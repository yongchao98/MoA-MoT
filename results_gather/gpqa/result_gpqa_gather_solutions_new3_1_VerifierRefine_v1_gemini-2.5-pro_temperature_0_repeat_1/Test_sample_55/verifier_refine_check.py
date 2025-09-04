import re
import numpy as np

def parse_state(state_str):
    """Parses a quantum state string like '|n,l,m>' into a tuple (n, l, m)."""
    match = re.search(r'\|\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(-?\d+)\s*\>', state_str)
    if not match:
        return None
    n, l, m = map(int, match.groups())
    return (n, l, m)

def check_dipole_transition(state1, state2):
    """Checks if a transition between two states is allowed by E1 selection rules."""
    n1, l1, m1 = state1
    n2, l2, m2 = state2
    
    # Rule: For decay, n must decrease.
    if n2 >= n1:
        return False, f"Transition fails: n must decrease (from {n1} to {n2})."
        
    # Rule: Δl = ±1
    if abs(l2 - l1) != 1:
        return False, f"Transition fails: Δl must be ±1 (was {l2 - l1})."
        
    # Rule: Δm = 0, ±1
    if abs(m2 - m1) > 1:
        return False, f"Transition fails: Δm must be 0 or ±1 (was {m2 - m1})."
        
    return True, "Transition is allowed."

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer by verifying the transition route
    and the associated probability.
    """
    # Extract the chosen option letter from the final answer
    final_answer_match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not final_answer_match:
        return "Incorrect: The final answer format is invalid. Expected '<<<X>>>'."
    
    chosen_option_letter = final_answer_match.group(1)
    
    # Define the options as provided in the question
    # Note: Option D has a typo in the prompt, which is handled here.
    options = {
        'A': {'route_str': '|3,0,0>->|2,1,1>->|1,0,0>', 'prob_str': '1/4'},
        'B': {'route_str': '|3,0,0>->|2,1,-1>->|1,0,0>', 'prob_str': '1/4'},
        'C': {'route_str': '|3,0,0>->|2,1,0>->|1,0,0>', 'prob_str': '1/3'},
        'D': {'route_str': '|3,0,0>->|2,1,0>->|1,0,0>', 'prob_str': '2/3'}
    }
    
    chosen_option = options.get(chosen_option_letter)
    if not chosen_option:
        return f"Incorrect: The chosen option '{chosen_option_letter}' is not valid."

    # --- Part 1: Verify the transition route ---
    route_str = chosen_option['route_str']
    states_str = route_str.split('->')
    
    try:
        initial_state = parse_state(states_str[0])
        intermediate_state = parse_state(states_str[1])
        final_state = parse_state(states_str[2])
        if not all([initial_state, intermediate_state, final_state]):
            raise ValueError("Could not parse all states.")
    except (ValueError, IndexError):
        return f"Incorrect: Could not parse the quantum states from the route '{route_str}' in option {chosen_option_letter}."

    # Check first transition
    is_allowed1, reason1 = check_dipole_transition(initial_state, intermediate_state)
    if not is_allowed1:
        return f"Incorrect: The route in option {chosen_option_letter} is invalid. {reason1}"

    # Check second transition
    is_allowed2, reason2 = check_dipole_transition(intermediate_state, final_state)
    if not is_allowed2:
        return f"Incorrect: The route in option {chosen_option_letter} is invalid. {reason2}"

    # --- Part 2: Verify the probability ---
    # The correct probability is 1/3, based on the symmetry of the initial l=0 state.
    correct_prob = 1/3
    
    # Get probability from the chosen option
    try:
        num, den = map(int, chosen_option['prob_str'].split('/'))
        option_prob = num / den
    except (ValueError, ZeroDivisionError):
        return f"Incorrect: Could not parse the probability '{chosen_option['prob_str']}' in option {chosen_option_letter}."

    # Compare probabilities
    if not np.isclose(option_prob, correct_prob):
        return (f"Incorrect: The probability for option {chosen_option_letter} is wrong. "
                f"The option gives {chosen_option['prob_str']}, but the correct probability is 1/3. "
                f"This is because the spherically symmetric initial state |3,0,0> decays with equal probability "
                f"to each of the three possible intermediate m' states.")

    # If all checks pass, the answer is correct.
    return "Correct"

# The user's provided response to check
llm_response = """
Here is a step-by-step analysis of the problem to determine the final answer.

### Step-by-step Analysis

1.  **Identify the Process:** The problem describes a two-step decay (a cascade) of a hydrogen atom from the initial state $|3,0,0\rangle$ to the final state $|1,0,0\rangle$. Each step is an electric dipole transition.

2.  **Recall Electric Dipole Selection Rules:** For a transition from a state $|n,l,m\rangle$ to $|n',l',m'\rangle$ to be allowed via an electric dipole interaction, the following selection rules must be met:
    *   The change in the orbital angular momentum quantum number must be $\Delta l = l' - l = \pm 1$.
    *   The change in the magnetic quantum number must be $\Delta m = m' - m = 0, \pm 1$.
    *   For spontaneous emission (decay), the energy must decrease, which means the principal quantum number must decrease, $n' < n$.

3.  **Analyze the First Transition:** The atom decays from the initial state $|3,0,0\rangle$ to an intermediate state $|n',l',m'\rangle$.
    *   **Initial State:** $n=3, l=0, m=0$.
    *   **Applying $\Delta l = \pm 1$:** Since the initial state has $l=0$, the intermediate state must have $l'=1$.
    *   **Applying $\Delta m = 0, \pm 1$:** Since the initial state has $m=0$, the intermediate state can have $m' = -1, 0, \text{ or } 1$.
    *   **Applying $n' < n$:** The principal quantum number must be less than 3. The options suggest the intermediate state has $n'=2$.
    *   **Conclusion:** The possible intermediate states are $|2,1,-1\rangle$, $|2,1,0\rangle$, and $|2,1,1\rangle$.

4.  **Analyze the Second Transition:** The atom decays from the intermediate state $|2,1,m'\rangle$ to the final state $|1,0,0\rangle$.
    *   **Intermediate State:** $n'=2, l'=1, m' \in \{-1, 0, 1\}$.
    *   **Final State:** $n_f=1, l_f=0, m_f=0$.
    *   **Applying $\Delta l = \pm 1$:** The change is $\Delta l = l_f - l' = 0 - 1 = -1$. This is an allowed transition.
    *   **Applying $\Delta m = 0, \pm 1$:** The change is $\Delta m = m_f - m' = 0 - m'$. This is allowed for all possible intermediate values of $m'$ (giving $\Delta m = +1, 0, -1$).
    *   **Conclusion:** All three possible decay routes are physically valid according to the selection rules.

5.  **Calculate the Transition Probability (Branching Ratio):** The probability of a specific cascade route is the product of the probabilities (branching ratios) of each step.
    *   **Probability of the first step ($|3,0,0\rangle \rightarrow |2,1,m'\rangle$):**
        The initial state $|3,0,0\rangle$ has $l=0$, meaning it is spherically symmetric. In the absence of external fields, there is no preferred direction in space. The three possible final states ($|2,1,-1\rangle$, $|2,1,0\rangle$, $|2,1,1\rangle$) are degenerate in energy. Due to the symmetry of the initial state, the probability of decaying into any of these three substates is equal. Since there are three equally likely outcomes, the probability for any specific one is $1/3$.
    *   **Probability of the second step ($|2,1,m'\rangle \rightarrow |1,0,0\rangle$):**
        Once the atom is in a $|2,1,m'\rangle$ state, it can only decay to a lower energy level. The only available lower level is $n=1$, which only has an $l=0$ subshell. The transition to $|1,0,0\rangle$ is the only possible dipole decay channel. Therefore, the probability of this second step occurring is 1.
    *   **Total Probability:** The total probability for any of the three specific routes is the product of the probabilities of the two steps: $P_{total} = P_{step1} \times P_{step2} = \frac{1}{3} \times 1 = \frac{1}{3}$.

6.  **Evaluate the Options:** We need to find the option that lists a valid route with the correct probability of $1/3$.
    *   A) $|3,0,0\rangle\rightarrow|2,1,1\rangle\rightarrow|1,0,0\rangle$ and $\frac{1}{4}$. The route is valid, but the probability is incorrect.
    *   B) $|3,0,0\rangle\rightarrow|2,1,-1\rangle\rightarrow|1,0,0\rangle$ and $\frac{1}{4}$. The route is valid, but the probability is incorrect.
    *   C) $|3,0,0\rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle$ and $\frac{1}{3}$. The route is valid, and the probability is correct.
    *   D) $\rangle\rightarrow|2,1,0\rangle\rightarrow|1,0,0\rangle$ and $\frac{2}{3}$. The route is valid (assuming a typo), but the probability is incorrect.

The only option that correctly identifies a valid transition route and its corresponding probability is C.

<<<C>>>
"""

# Run the check
result = check_correctness("", llm_response)
print(result)