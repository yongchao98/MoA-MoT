import sympy
from sympy.physics.quantum.cg import CG
import numpy as np

def check_correctness():
    """
    Checks the correctness of the provided answer for the H-atom decay problem.
    It calculates the theoretical probabilities for the decay paths and compares them
    with the options given in the question.
    """
    # Define the initial state and the possible intermediate states
    initial_state = {'l': 0, 'm': 0}
    intermediate_l = 1
    possible_intermediate_m = [-1, 0, 1]

    # Options from the multiple-choice question
    options = {
        'A': {'intermediate_m': 0, 'prob': 1/3},
        'B': {'intermediate_m': 1, 'prob': 1/4},
        'C': {'intermediate_m': 0, 'prob': 2/3},
        'D': {'intermediate_m': -1, 'prob': 1/4},
    }
    llm_answer = 'A'

    # Calculate the relative rates using the Wigner-Eckart theorem.
    # The rate is proportional to the square of the Clebsch-Gordan coefficient.
    # The radial integral is the same for all m' and cancels in the branching ratio.
    relative_rates = {}
    for m_k in possible_intermediate_m:
        # For the transition |l_i, m_i> -> |l_k, m_k>, the photon has q = m_k - m_i
        q = m_k - initial_state['m']
        # The relevant CG coefficient is <l_i, m_i; 1, q | l_k, m_k>
        cg_coeff = CG(j1=initial_state['l'], m1=initial_state['m'], 
                      j2=1, m2=q, 
                      j3=intermediate_l, m3=m_k).doit()
        relative_rates[m_k] = float(cg_coeff**2)

    total_relative_rate = sum(relative_rates.values())

    if total_relative_rate == 0:
        return "Calculation error: Total rate is zero, which means the transition is forbidden. This contradicts the problem statement."

    # Calculate the branching ratio (probability) for each path
    calculated_probabilities = {m: rate / total_relative_rate for m, rate in relative_rates.items()}

    # Check which options are correct based on our calculation
    correct_options = []
    for key, value in options.items():
        m_intermediate = value['intermediate_m']
        prob_given = value['prob']
        
        if m_intermediate in calculated_probabilities:
            prob_calculated = calculated_probabilities[m_intermediate]
            if np.isclose(prob_given, prob_calculated):
                correct_options.append(key)

    # Final verification of the LLM's answer
    if llm_answer in correct_options:
        if len(correct_options) == 1:
            return "Correct"
        else:
            # This case would mean the question has multiple correct answers
            return f"The LLM's answer {llm_answer} is correct, but incomplete. Options {correct_options} are all correct."
    else:
        # The LLM's answer is incorrect. Provide the reason.
        llm_m = options[llm_answer]['intermediate_m']
        llm_prob = options[llm_answer]['prob']
        calc_prob = calculated_probabilities[llm_m]
        return (f"The answer {llm_answer} is incorrect. "
                f"For the transition path via |2,1,{llm_m}>, the calculated probability is {calc_prob:.4f} (which is 1/3), "
                f"but the answer states the probability is {llm_prob:.4f}.")

# Execute the check
result = check_correctness()
print(result)