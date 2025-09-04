import math

class QuantumState:
    """A simple class to represent a quantum state |n, l, m>."""
    def __init__(self, n, l, m):
        # Basic validation of quantum numbers
        if not (isinstance(n, int) and n > 0 and
                isinstance(l, int) and l >= 0 and l < n and
                isinstance(m, int) and abs(m) <= l):
            raise ValueError(f"Invalid quantum numbers for state |{n},{l},{m}>")
        self.n = n
        self.l = l
        self.m = m

    def __repr__(self):
        return f"|{self.n},{self.l},{self.m}>"

def is_dipole_transition_allowed(initial: QuantumState, final: QuantumState):
    """
    Checks if an electric dipole (E1) transition is allowed between two states.
    Returns a tuple (bool, str) indicating if allowed and the reason.
    """
    # Constraint: For spontaneous emission (decay), energy must decrease (n must decrease).
    if final.n >= initial.n:
        return False, f"Energy did not decrease (n_final={final.n} >= n_initial={initial.n})"

    # Constraint: The change in orbital angular momentum quantum number must be +/- 1.
    delta_l = final.l - initial.l
    if abs(delta_l) != 1:
        return False, f"Selection rule for l is not satisfied (Δl = {delta_l}, but must be ±1)"

    # Constraint: The change in magnetic quantum number must be 0 or +/- 1.
    delta_m = final.m - initial.m
    if abs(delta_m) > 1:
        return False, f"Selection rule for m is not satisfied (Δm = {delta_m}, but must be 0 or ±1)"

    return True, "Allowed"

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying the physical constraints.
    """
    # --- Problem Definition ---
    # The question asks for the correct decay route and probability for the H-atom
    # transitioning from |3,0,0> to |1,0,0> via two dipole transitions.
    initial_state = QuantumState(3, 0, 0)
    final_state = QuantumState(1, 0, 0)

    # --- Options from the Question ---
    # Note: Option B in the prompt has a typo, but the intended path is clear.
    options = {
        'A': {'path': [QuantumState(3, 0, 0), QuantumState(2, 1, -1), QuantumState(1, 0, 0)], 'prob': 1/4},
        'B': {'path': [QuantumState(3, 0, 0), QuantumState(2, 1, 0), QuantumState(1, 0, 0)], 'prob': 2/3},
        'C': {'path': [QuantumState(3, 0, 0), QuantumState(2, 1, 0), QuantumState(1, 0, 0)], 'prob': 1/3},
        'D': {'path': [QuantumState(3, 0, 0), QuantumState(2, 1, 1), QuantumState(1, 0, 0)], 'prob': 1/4}
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_key = 'C'
    
    chosen_option = options.get(llm_answer_key)
    if not chosen_option:
        return f"Invalid answer key '{llm_answer_key}'. Must be one of {list(options.keys())}."

    path = chosen_option['path']
    given_prob = chosen_option['prob']
    
    # --- Verification Steps ---

    # 1. Verify the first transition in the path
    intermediate_state = path[1]
    is_allowed, reason = is_dipole_transition_allowed(initial_state, intermediate_state)
    if not is_allowed:
        return f"The first transition in the path {initial_state} -> {intermediate_state} is not allowed. Reason: {reason}"

    # 2. Verify the second transition in the path
    is_allowed, reason = is_dipole_transition_allowed(intermediate_state, final_state)
    if not is_allowed:
        return f"The second transition in the path {intermediate_state} -> {final_state} is not allowed. Reason: {reason}"

    # 3. Verify the probability
    # The probability of the cascade is determined by the branching ratio of the first step.
    # The initial state |3,0,0> is spherically symmetric (l=0). Therefore, it decays with
    # equal probability to each of the possible m' sublevels of the intermediate l'=1 state.
    
    # Find all possible intermediate states for the first decay step.
    # From selection rules, the intermediate state must have l'=1. The question implies n'=2.
    n_intermediate, l_intermediate = 2, 1
    possible_m_values = range(-l_intermediate, l_intermediate + 1)
    
    # The number of equally likely decay channels for the first step.
    num_channels = len(list(possible_m_values))
    
    # The probability of the second step is 1, as |2,1,m'> can only decay to |1,0,0> via E1.
    # So, the total probability of the cascade is the probability of the first step.
    correct_prob = 1.0 / num_channels
    
    # Compare the calculated probability with the one given in the option.
    if not math.isclose(correct_prob, given_prob):
        return (f"The transition path is valid, but the probability is incorrect. "
                f"Calculated probability is 1/{num_channels} ≈ {correct_prob:.4f}, "
                f"but the answer states it is {given_prob:.4f}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_answer()
print(result)