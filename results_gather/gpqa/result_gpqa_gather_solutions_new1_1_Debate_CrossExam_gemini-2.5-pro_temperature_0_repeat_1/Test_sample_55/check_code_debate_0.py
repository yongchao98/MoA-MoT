import re
import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer to the H-atom decay problem.

    The function verifies two main constraints:
    1.  The validity of the transition path based on electric dipole selection rules.
    2.  The correctness of the transition probability based on the symmetry of the initial state.
    """

    # The final answer provided by the LLM to be checked.
    llm_answer = "C"

    # The options as described in the LLM's reasoning section.
    # Note: The option lettering in the various LLM answers is inconsistent.
    # We use the lettering from the final provided answer's reasoning.
    options = {
        "A": {"route_str": "|3,0,0>→|2,1,0>→|1,0,0>", "prob_str": "2/3"},
        "B": {"route_str": "|3,0,0>→|2,1,1>→|1,0,0>", "prob_str": "1/4"},
        "C": {"route_str": "|3,0,0>→|2,1,0>→|1,0,0>", "prob_str": "1/3"},
        "D": {"route_str": "|3,0,0>→|2,1,-1>→|1,0,0>", "prob_str": "1/4"}
    }

    def parse_state(state_str):
        """Parses a state string like '|n,l,m>' into a tuple (n, l, m)."""
        try:
            numbers = re.findall(r'-?\d+', state_str)
            return tuple(map(int, numbers))
        except (TypeError, ValueError):
            return None

    def check_selection_rules(state1, state2):
        """Checks if a dipole transition from state1 to state2 is allowed."""
        if not state1 or not state2:
            return False, "Invalid state format"
        n1, l1, m1 = state1
        n2, l2, m2 = state2

        # For spontaneous emission, energy must decrease (n must decrease)
        if n2 >= n1:
            return False, f"n does not decrease ({n1} -> {n2})"

        # Selection rule for l: Δl = ±1
        if abs(l2 - l1) != 1:
            return False, f"Δl is not ±1 (l changed from {l1} to {l2})"

        # Selection rule for m: Δm = 0, ±1
        if abs(m2 - m1) > 1:
            return False, f"Δm is not 0 or ±1 (m changed from {m1} to {m2})"

        return True, "Allowed"

    # Theoretical probability for any single path from a spherically symmetric state
    # to one of three degenerate sublevels.
    theoretical_prob = 1/3

    correct_options = []

    for option_key, value in options.items():
        route_str = value["route_str"]
        prob_str = value["prob_str"]

        states_str = route_str.split('→')
        if len(states_str) != 3:
            continue

        initial_state = parse_state(states_str[0])
        intermediate_state = parse_state(states_str[1])
        final_state = parse_state(states_str[2])

        try:
            given_prob = eval(prob_str)
        except Exception:
            continue

        # Constraint 1: Check if the path is valid
        allowed1, _ = check_selection_rules(initial_state, intermediate_state)
        allowed2, _ = check_selection_rules(intermediate_state, final_state)
        is_valid_path = allowed1 and allowed2

        # Constraint 2: Check if the probability is correct
        is_correct_prob = math.isclose(given_prob, theoretical_prob)

        if is_valid_path and is_correct_prob:
            correct_options.append(option_key)

    # Final Evaluation
    if llm_answer not in options:
        return f"The provided answer '{llm_answer}' is not one of the options A, B, C, D."

    if llm_answer in correct_options:
        if len(correct_options) == 1:
            return "Correct"
        else:
            # This case should not happen for this specific problem
            return f"The answer '{llm_answer}' is correct, but the question is ambiguous as options {correct_options} are also correct."
    else:
        # The LLM's answer is wrong. Provide a reason.
        if not correct_options:
             return "Incorrect. None of the options satisfy all constraints. The provided answer is wrong."
        
        correct_answer = correct_options[0]
        chosen_option_data = options[llm_answer]
        
        # Re-check the chosen option to generate a specific error message
        states_str = chosen_option_data["route_str"].split('→')
        initial = parse_state(states_str[0])
        intermediate = parse_state(states_str[1])
        final = parse_state(states_str[2])
        prob = eval(chosen_option_data["prob_str"])

        allowed1, reason1 = check_selection_rules(initial, intermediate)
        allowed2, reason2 = check_selection_rules(intermediate, final)

        if not allowed1:
            return f"Incorrect. The answer is {llm_answer}, but its path is invalid. The first transition is not allowed because {reason1}. The correct answer is {correct_answer}."
        if not allowed2:
            return f"Incorrect. The answer is {llm_answer}, but its path is invalid. The second transition is not allowed because {reason2}. The correct answer is {correct_answer}."
        if not math.isclose(prob, theoretical_prob):
            return f"Incorrect. The answer is {llm_answer}. The transition path is physically allowed, but the probability is wrong. The given probability is {chosen_option_data['prob_str']}, but the theoretical probability should be 1/3. The correct answer is {correct_answer}."
        
        return f"Incorrect. The correct answer is {correct_answer}, but the provided answer was {llm_answer}."

# Run the check and print the result
result = check_correctness()
print(result)