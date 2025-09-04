import math

class QuantumState:
    """A simple class to represent a quantum state |n,l,m>."""
    def __init__(self, n, l, m):
        if not (isinstance(n, int) and n > 0 and
                isinstance(l, int) and 0 <= l < n and
                isinstance(m, int) and -l <= m <= l):
            raise ValueError(f"Invalid quantum numbers for state |{n},{l},{m}>")
        self.n = n
        self.l = l
        self.m = m

    def __repr__(self):
        return f"|{self.n},{self.l},{self.m}>"

def check_dipole_selection_rules(initial_state, final_state):
    """Checks if a transition is allowed by electric dipole selection rules."""
    delta_l = final_state.l - initial_state.l
    delta_m = final_state.m - initial_state.m

    if abs(delta_l) != 1:
        return False, f"Invalid transition from {initial_state} to {final_state}: delta_l = {delta_l}, but must be +/- 1."
    if abs(delta_m) > 1:
        return False, f"Invalid transition from {initial_state} to {final_state}: delta_m = {delta_m}, but must be 0 or +/- 1."
    return True, ""

def calculate_cascade_probability(initial, intermediate, final):
    """
    Calculates the probability of a two-step cascade decay, applying physics principles.
    """
    # 1. Check if the path is allowed by selection rules
    is_allowed1, reason1 = check_dipole_selection_rules(initial, intermediate)
    if not is_allowed1:
        return 0, f"First step ({initial} -> {intermediate}) is forbidden: {reason1}"
    
    is_allowed2, reason2 = check_dipole_selection_rules(intermediate, final)
    if not is_allowed2:
        return 0, f"Second step ({intermediate} -> {final}) is forbidden: {reason2}"

    # 2. Calculate probability of the first step
    # For an initial l=0 state, the probability is distributed equally among all
    # allowed final m' states for the destination l'.
    prob_step1 = 0
    if initial.l == 0:
        # The destination l' is intermediate.l. The number of m' states is 2*l'+1.
        num_channels = 2 * intermediate.l + 1
        prob_step1 = 1.0 / num_channels
    else:
        # A more complex calculation involving Clebsch-Gordan coefficients would be needed,
        # but is not required for this specific problem.
        return None, "Probability calculation for initial l!=0 is not implemented."

    # 3. Probability of the second step is 1, as it's the only allowed decay to a lower energy level.
    prob_step2 = 1.0
    
    total_prob = prob_step1 * prob_step2
    return total_prob, "Calculation successful"

def check_llm_answer():
    """
    Checks the correctness of the LLM's answer by applying physical constraints.
    """
    # Define the problem's states and options
    initial_state = QuantumState(3, 0, 0)
    final_state = QuantumState(1, 0, 0)

    options = {
        'A': {'intermediate': QuantumState(2, 1, 1), 'prob': 1/4},
        'B': {'intermediate': QuantumState(2, 1, -1), 'prob': 1/4},
        'C': {'intermediate': QuantumState(2, 1, 0), 'prob': 2/3},
        'D': {'intermediate': QuantumState(2, 1, 0), 'prob': 1/3},
    }
    
    llm_provided_answer = 'D'
    
    # Find the theoretically correct option
    correct_option_key = None
    for key, data in options.items():
        calculated_prob, reason = calculate_cascade_probability(
            initial_state, data['intermediate'], final_state
        )
        
        if reason != "Calculation successful":
            # This path is invalid, so it can't be the correct option
            continue

        # Check if the calculated probability matches the option's probability
        if math.isclose(calculated_prob, data['prob']):
            if correct_option_key is not None:
                # This would mean multiple options are correct, indicating a flawed question
                return "Error: Multiple options were found to be correct."
            correct_option_key = key

    if correct_option_key is None:
        return "Error: No option was found to be correct based on physics principles."

    # Compare the LLM's answer with the derived correct answer
    if llm_provided_answer == correct_option_key:
        return "Correct"
    else:
        correct_prob = options[correct_option_key]['prob']
        llm_prob = options[llm_provided_answer]['prob']
        return (f"Incorrect. The provided answer is '{llm_provided_answer}', but the correct answer is '{correct_option_key}'.\n"
                f"Reason: The probability for the transition path in option '{llm_provided_answer}' is {llm_prob:.3f}, "
                f"but the calculated probability based on the Wigner-Eckart theorem is {correct_prob:.3f}.")

# Run the check
result = check_llm_answer()
print(result)