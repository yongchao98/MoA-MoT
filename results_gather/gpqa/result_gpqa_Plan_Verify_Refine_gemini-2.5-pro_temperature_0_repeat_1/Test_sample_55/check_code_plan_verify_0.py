import sympy
from sympy.physics.wigner import clebsch_gordan
from fractions import Fraction

def check_correctness_of_llm_answer():
    """
    Checks the correctness of the LLM's answer for the H-atom decay problem.

    The logic involves:
    1. Defining a function to check E1 selection rules.
    2. Identifying all possible intermediate states for the two-step decay.
    3. Calculating the branching ratios for the first decay step using Clebsch-Gordan coefficients.
       The probability of the entire cascade is equal to this branching ratio, as the second
       decay step has a probability of 1 (it's the only allowed E1 decay channel).
    4. Comparing the calculated path probabilities with the options provided in the question.
    5. Verifying that the LLM's chosen answer (B) is consistent with the calculations.
    """

    # --- Problem Definition ---
    initial_state = (3, 0, 0)
    final_state = (1, 0, 0)
    llm_answer_key = "B"

    # Parse the options from the question. Note the typo in D is interpreted as in the LLM answer.
    options = {
        "A": {"path": [(3, 0, 0), (2, 1, 1), (1, 0, 0)], "prob": Fraction(1, 4)},
        "B": {"path": [(3, 0, 0), (2, 1, 0), (1, 0, 0)], "prob": Fraction(1, 3)},
        "C": {"path": [(3, 0, 0), (2, 1, -1), (1, 0, 0)], "prob": Fraction(1, 4)},
        "D": {"path": [(3, 0, 0), (2, 1, 0), (1, 0, 0)], "prob": Fraction(2, 3)},
    }

    # --- Helper Function for Selection Rules ---
    def check_e1_transition(initial, final):
        n_i, l_i, m_i = initial
        n_f, l_f, m_f = final
        if not (n_f < n_i and abs(l_f - l_i) == 1 and abs(m_f - m_i) <= 1):
            return False
        return True

    # --- Step 1: Identify all valid intermediate states ---
    possible_intermediates = []
    # Iterate through all possible quantum numbers for an intermediate state
    for n_prime in range(1, initial_state[0]):
        for l_prime in range(n_prime):
            for m_prime in range(-l_prime, l_prime + 1):
                intermediate = (n_prime, l_prime, m_prime)
                # Check if both steps of the decay are allowed
                if check_e1_transition(initial_state, intermediate) and \
                   check_e1_transition(intermediate, final_state):
                    possible_intermediates.append(intermediate)
    
    # Expected intermediates are the 2p states
    expected_intermediates = [(2, 1, -1), (2, 1, 0), (2, 1, 1)]
    if sorted(possible_intermediates) != sorted(expected_intermediates):
        return f"Logic Error: The code failed to identify the correct set of intermediate states. Found {possible_intermediates}."

    # --- Step 2: Calculate branching ratios for the first decay step ---
    # The rate is proportional to the square of the Clebsch-Gordan coefficient.
    # The radial integral part is common to all transitions and cancels in the ratio.
    # Transition: |l=0, m=0> -> |l'=1, m'>
    rates = {}
    total_rate_proportional = 0
    for inter_state in possible_intermediates:
        l_i, m_i = initial_state[1], initial_state[2]
        l_f, m_f = inter_state[1], inter_state[2]
        
        # The dipole operator is a rank 1 tensor. The relevant CG coefficient is
        # <l_i, m_i; 1, q | l_f, m_f>, where q = m_f - m_i.
        q = m_f - m_i
        cg_coeff = clebsch_gordan(l_i, 1, l_f, m_i, q, m_f)
        
        # Rate is proportional to the coefficient squared
        rate = float(cg_coeff**2)
        rates[inter_state] = rate
        total_rate_proportional += rate

    if total_rate_proportional == 0:
        return "Calculation Error: Total rate is zero, cannot compute branching ratios."

    # --- Step 3: Verify the LLM's chosen answer ---
    selected_option = options.get(llm_answer_key)
    if not selected_option:
        return f"The LLM's answer '{llm_answer_key}' is not a valid option key."

    path_to_check = selected_option["path"]
    intermediate_to_check = path_to_check[1]
    
    # Check if the path is valid
    if not (check_e1_transition(path_to_check[0], path_to_check[1]) and \
            check_e1_transition(path_to_check[1], path_to_check[2])):
        return f"The path in the selected option {llm_answer_key} is not a valid E1 decay path."

    # Check the probability
    calculated_prob = Fraction(rates[intermediate_to_check] / total_rate_proportional).limit_denominator()
    claimed_prob = selected_option["prob"]

    if calculated_prob == claimed_prob:
        return "Correct"
    else:
        return (f"Incorrect probability for option {llm_answer_key}. "
                f"The calculated probability for the path {path_to_check} is {calculated_prob}, "
                f"but the option claims {claimed_prob}.")

# Execute the check and print the result
result = check_correctness_of_llm_answer()
print(result)