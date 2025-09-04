import math

def check_angular_momentum_coupling_probability(initial_state, measured_outcome):
    """
    Checks the correctness of a probability calculation in angular momentum coupling.

    This function verifies the answer based on the fundamental selection rule for
    Clebsch-Gordan coefficients: m_total = m1 + m2. If this rule is violated,
    the probability of the measurement is zero.

    Args:
        initial_state (dict): A dictionary with the quantum numbers of the initial
                              coupled state, e.g., {'l1': 1, 'l2': 1, 'l': 2, 'm': -1}.
        measured_outcome (dict): A dictionary with the quantum numbers of the
                                 measured uncoupled state, e.g., {'m1': -1, 'm2': -1}.

    Returns:
        tuple: A tuple containing:
               - bool: True if the calculated probability matches the expected outcome, False otherwise.
               - str: A message explaining the reasoning.
    """
    l1 = initial_state.get('l1')
    l2 = initial_state.get('l2')
    l = initial_state.get('l')
    m_total = initial_state.get('m')

    m1_measured = measured_outcome.get('m1')
    m2_measured = measured_outcome.get('m2')

    # --- Constraint Checks ---
    # 1. Check if the measured quantum numbers are valid for the given l1 and l2.
    if abs(m1_measured) > l1:
        return False, f"Measured m1={m1_measured} is not a valid state for a particle with l1={l1}."
    if abs(m2_measured) > l2:
        return False, f"Measured m2={m2_measured} is not a valid state for a particle with l2={l2}."

    # 2. Check if the total angular momentum quantum numbers are valid.
    if l not in range(abs(l1 - l2), l1 + l2 + 1):
         return False, f"The total angular momentum l={l} is not possible for coupling l1={l1} and l2={l2}."
    if abs(m_total) > l:
        return False, f"The initial state m={m_total} is not a valid state for a system with l={l}."

    # 3. The core selection rule: m_total must equal the sum of m1 and m2.
    m_sum_measured = m1_measured + m2_measured
    
    # The probability is the square of the Clebsch-Gordan coefficient.
    # The coefficient is zero if the selection rule is violated.
    if m_total != m_sum_measured:
        calculated_probability = 0
        reason = (f"The selection rule m_total = m1 + m2 is violated. "
                  f"The initial state has m_total = {m_total}, but the measured state "
                  f"has a sum m1 + m2 = {m1_measured} + {m2_measured} = {m_sum_measured}. "
                  f"Therefore, the probability must be 0.")
        return calculated_probability, reason
    else:
        # If the rule holds, the probability is non-zero, but calculating the exact
        # value requires a Clebsch-Gordan table. For this problem, this branch is not needed.
        # We can't determine the exact non-zero probability without more information.
        # However, the question is designed to be solved by the selection rule.
        return "Non-zero, requires Clebsch-Gordan coefficients", "The selection rule m_total = m1 + m2 is satisfied. The probability is non-zero."


# --- Problem Setup ---
# From the question: "two electrons are in p orbital" -> l1=1, l2=1
# Initial state: |l1, l2, l, m> = |1,1, 2, -1>
initial_state_params = {'l1': 1, 'l2': 1, 'l': 2, 'm': -1}

# Measurement: "eigenvalues of both L1z and L2z as -ħ" -> m1=-1, m2=-1
measured_outcome_params = {'m1': -1, 'm2': -1}

# The final answer provided by the LLM analysis
llm_answer_choice = 'A'
answer_options = {'A': 0, 'B': 1/2, 'C': 1, 'D': 2/3}
llm_answer_value = answer_options.get(llm_answer_choice)

# --- Verification ---
calculated_prob, reason_message = check_angular_momentum_coupling_probability(
    initial_state_params,
    measured_outcome_params
)

print(f"Analyzing the problem...")
print(f"Initial state: m_total = {initial_state_params['m']}")
print(f"Measured outcome: m1 = {measured_outcome_params['m1']}, m2 = {measured_outcome_params['m2']}")
print("-" * 20)
print(f"Reasoning: {reason_message}")
print("-" * 20)
print(f"Calculated Probability: {calculated_prob}")
print(f"LLM's Answer ('{llm_answer_choice}'): {llm_answer_value}")
print("-" * 20)

# --- Final Verdict ---
if llm_answer_value is None:
    print(f"Error: The provided answer choice '{llm_answer_choice}' is not a valid option.")
elif isinstance(calculated_prob, str):
     print(f"Cannot definitively check correctness. The calculated probability is {calculated_prob}.")
# Using math.isclose for robust floating-point comparison
elif math.isclose(llm_answer_value, calculated_prob):
    print("Correct")
else:
    print(f"Incorrect. The provided answer corresponds to a probability of {llm_answer_value}, but the calculated probability is {calculated_prob}.")
