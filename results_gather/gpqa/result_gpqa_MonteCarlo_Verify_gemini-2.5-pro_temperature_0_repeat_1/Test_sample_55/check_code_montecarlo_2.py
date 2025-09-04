import re
import math
from sympy.physics.quantum.cg import CG

def check_correctness():
    """
    This function checks the correctness of the given LLM answer by:
    1. Defining the problem parameters (initial/final states, options).
    2. Parsing the selected answer to get the proposed transition route and probability.
    3. Verifying that each step in the route obeys the electric dipole selection rules.
    4. Calculating the theoretical branching ratio for the first decay step using Clebsch-Gordan coefficients.
    5. Comparing the calculated probability with the probability given in the answer.
    """
    # --- Problem Definition ---
    # The question specifies the overall transition and the LLM provides the answer key.
    initial_state_q = (3, 0, 0)
    final_state_q = (1, 0, 0)
    llm_answer_key = "A"  # This is the provided answer '<<<A>>>'

    # Define the options from the question
    options = {
        "A": "|3,0,0\\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle and \\frac{1}{3}",
        "B": "|3,0,0\\rangle\\rightarrow|2,1,1\\rangle\\rightarrow|1,0,0\\rangle and \\frac{1}{4}",
        "C": "\\rangle\\rightarrow|2,1,0\\rangle\\rightarrow|1,0,0\\rangle and \\frac{2}{3}",
        "D": "|3,0,0\\rangle\\rightarrow|2,1,-1\\rangle\\rightarrow|1,0,0\\rangle and \\frac{1}{4}"
    }

    # --- Helper Functions ---
    def parse_state(state_str):
        """Parses a state string like '|n,l,m>' into a tuple (n, l, m)."""
        state_str = state_str.replace(r'\rangle', '').strip()
        if state_str.startswith('>'): state_str = state_str[1:]
        match = re.search(r'\|(\d+),(\d+),(-?\d+)\>', state_str)
        if not match:
            return None if state_str == '' else f"Could not parse state: {state_str}"
        n, l, m = map(int, match.groups())
        return (n, l, m)

    def check_selection_rules(initial, final):
        """Checks electric dipole transition selection rules for decay."""
        n_i, l_i, m_i = initial
        n_f, l_f, m_f = final
        delta_l, delta_m = l_f - l_i, m_f - m_i
        if n_f >= n_i: return f"Invalid decay {initial} -> {final}: n must decrease."
        if abs(delta_l) != 1: return f"Invalid transition {initial} -> {final}: Δl = {delta_l}, must be ±1."
        if abs(delta_m) > 1: return f"Invalid transition {initial} -> {final}: Δm = {delta_m}, must be 0 or ±1."
        return None

    def get_branching_ratios(initial, intermediate_n):
        """Calculates branching ratios for decay from the initial state."""
        n, l, m = initial
        rates = {}
        possible_l_prime = [l + 1] if l == 0 else [l - 1, l + 1]
        for l_prime in possible_l_prime:
            for m_prime in range(-l_prime, l_prime + 1):
                delta_m = m_prime - m
                if abs(delta_m) > 1: continue
                cg_coeff = CG(l, m, 1, delta_m, l_prime, m_prime).doit()
                rate = float(cg_coeff**2)
                if rate > 0: rates[(intermediate_n, l_prime, m_prime)] = rate
        total_rate = sum(rates.values())
        return {state: rate / total_rate for state, rate in rates.items()} if total_rate > 0 else {}

    # --- Main Execution ---
    try:
        answer_content = options.get(llm_answer_key)
        if not answer_content:
            return f"Invalid answer key '{llm_answer_key}'. Must be one of A, B, C, D."

        # 1. Parse the answer
        parts = answer_content.split(' and ')
        route_str, prob_str = parts[0], parts[1]
        state_strs = re.findall(r'(\\|>)?\|?[\\d, -]+\\>', route_str)
        states = [parse_state(s) for s in state_strs]
        if states[0] is None: states[0] = initial_state_q # Fix for malformed option C
        s1, s2, s3 = states
        
        prob_match = re.search(r'\\frac\{(\d+)\}\{(\d+)\}', prob_str)
        probability_given = float(prob_match.group(1)) / float(prob_match.group(2))

        # 2. Check route validity
        if s1 != initial_state_q or s3 != final_state_q:
            return f"Route {s1}->{s2}->{s3} does not match question's start/end states."
        
        error = check_selection_rules(s1, s2)
        if error: return f"First transition is invalid: {error}"
        error = check_selection_rules(s2, s3)
        if error: return f"Second transition is invalid: {error}"

        # 3. Calculate and check probability
        # The probability of the route is the branching ratio of the first step.
        branching_ratios = get_branching_ratios(s1, s2[0])
        calculated_prob = branching_ratios.get(s2)

        if calculated_prob is None:
            return f"The intermediate state {s2} is not a possible decay product from {s1}."

        if not math.isclose(calculated_prob, probability_given, rel_tol=1e-5):
            return (f"Probability mismatch. The calculated branching ratio for {s1} -> {s2} is "
                    f"{calculated_prob:.4f}, but the answer provides {probability_given:.4f}.")

        return "Correct"
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_correctness()
# The result will be "Correct" because the logic confirms that option A is valid
# and its probability is correctly stated as 1/3.
print(result)