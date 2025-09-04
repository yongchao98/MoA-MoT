import sympy
from sympy.physics.quantum.cg import CG
from sympy import S

def check_hydrogen_decay_answer():
    """
    This function checks the correctness of the provided answer for the H-atom decay problem.
    It verifies both the selection rules and the calculated transition probability.
    """
    
    # --- Problem Definition ---
    # States are represented as dictionaries: {'n': val, 'l': val, 'm': val}
    initial_state = {'n': 3, 'l': 0, 'm': 0}
    final_state = {'n': 1, 'l': 0, 'm': 0}

    # The options given in the question. Probabilities are converted to SymPy fractions for exact comparison.
    options = {
        'A': {'intermediate': {'n': 2, 'l': 1, 'm': 1}, 'prob': S('1/4')},
        'B': {'intermediate': {'n': 2, 'l': 1, 'm': -1}, 'prob': S('1/4')},
        'C': {'intermediate': {'n': 2, 'l': 1, 'm': 0}, 'prob': S('2/3')},
        'D': {'intermediate': {'n': 2, 'l': 1, 'm': 0}, 'prob': S('1/3')}
    }
    
    # The answer provided by the other LLM that we need to check.
    llm_answer_choice = 'D'
    
    # --- Verification Step 1: Selection Rules ---
    # For an electric dipole transition, the selection rules are:
    # delta_l = +/- 1
    # delta_m = 0, +/- 1
    # For spontaneous emission, energy must decrease (delta_n > 0 for n_initial - n_final).
    
    llm_option_data = options[llm_answer_choice]
    
    # Transition 1: initial -> intermediate
    s1_i, s1_f = initial_state, llm_option_data['intermediate']
    delta_l1 = abs(s1_f['l'] - s1_i['l'])
    delta_m1 = abs(s1_f['m'] - s1_i['m'])
    delta_n1 = s1_i['n'] - s1_f['n']
    
    if not (delta_l1 == 1 and delta_m1 <= 1 and delta_n1 > 0):
        return f"Incorrect. The first transition in the proposed path {llm_answer_choice} violates the electric dipole selection rules."

    # Transition 2: intermediate -> final
    s2_i, s2_f = llm_option_data['intermediate'], final_state
    delta_l2 = abs(s2_f['l'] - s2_i['l'])
    delta_m2 = abs(s2_f['m'] - s2_i['m'])
    delta_n2 = s2_i['n'] - s2_f['n']

    if not (delta_l2 == 1 and delta_m2 <= 1 and delta_n2 > 0):
        return f"Incorrect. The second transition in the proposed path {llm_answer_choice} violates the electric dipole selection rules."

    # --- Verification Step 2: Transition Probability ---
    # The total probability of the cascade is P(step 1) * P(step 2).
    # P(step 2) is 1, as the decay from |2,1,m'> to |1,0,0> is the only allowed spontaneous emission.
    # Therefore, the total probability is the branching ratio of the first step.
    
    # We calculate the branching ratio for the decay |3,0,0> -> |2,1,m'>.
    # Due to the spherical symmetry of the initial |l=0> state, the decay rates to each of the
    # |l=1, m'> sub-levels are equal. The probability for each is 1/3.
    # We can formally verify this using Clebsch-Gordan coefficients.
    
    l_initial = S(0)
    m_initial = S(0)
    l_intermediate = S(1)
    
    # The possible final m states for the first transition
    possible_m_intermediate = [-1, 0, 1]
    relative_rates = {}
    
    for m_prime in possible_m_intermediate:
        # The transition rate is proportional to the square of the Clebsch-Gordan coefficient
        # <l_i, m_i; 1, q | l_f, m_f>, where q = m_f - m_i.
        # Here, this is <0, 0; 1, m' | 1, m'>.
        cg_coeff = CG(l_initial, m_initial, 1, m_prime, l_intermediate, m_prime)
        relative_rates[m_prime] = cg_coeff.doit()**2
        
    total_relative_rate = sum(relative_rates.values())
    
    # Calculate the branching ratio for the specific path chosen in the answer
    m_of_interest = llm_option_data['intermediate']['m']
    calculated_probability = relative_rates[m_of_interest] / total_relative_rate
    
    # --- Final Conclusion ---
    stated_probability = llm_option_data['prob']
    
    if calculated_probability == stated_probability:
        return "Correct"
    else:
        # Find the actual correct option for a more informative error message.
        correct_option = "None"
        for name, opt in options.items():
            m_inter = opt['intermediate']['m']
            prob_inter = opt['prob']
            if (relative_rates[m_inter] / total_relative_rate) == prob_inter:
                correct_option = name
                break

        return (f"Incorrect. The answer claims the probability for option {llm_answer_choice} is {stated_probability}, "
                f"but the calculated physical probability is {calculated_probability}.\n"
                f"The branching ratio for the decay |3,0,0> -> |2,1,m'> is 1/3 for each m' in [-1, 0, 1].\n"
                f"The correct option is {correct_option}.")

# Run the check
result = check_hydrogen_decay_answer()
print(result)