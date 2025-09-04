import sympy
from sympy.physics.quantum.cg import CG
import numpy as np

def check_correctness_of_answer():
    """
    This function checks the correctness of the given answer to the quantum mechanics question.

    The question asks for the correct decay route and probability for a two-dipole transition
    from the |3,0,0> state to the |1,0,0> state in a hydrogen atom.

    The provided answer is A: |3,0,0> -> |2,1,0> -> |1,0,0> with probability 1/3.

    The check involves three steps:
    1.  Verify that the proposed transition path is allowed by the electric dipole selection rules.
    2.  Calculate the branching ratios for all possible first-step decay channels from |3,0,0>.
        The probability of a specific route is the branching ratio of the first transition,
        as the second transition (|2p> -> |1s>) is the only allowed decay channel from the
        intermediate state.
    3.  Compare the calculated probability for the path in option A with the probability
        given in option A.
    """

    # --- Step 0: Define problem parameters and the answer to check ---
    initial_state = (3, 0, 0) # (n, l, m)
    final_state = (1, 0, 0)
    
    # The answer provided by the other LLM
    llm_answer_key = 'A'
    
    # Define the options from the question
    options = {
        'A': {'intermediate': (2, 1, 0), 'prob': 1/3},
        'B': {'intermediate': (2, 1, 1), 'prob': 1/4},
        'C': {'intermediate': (2, 1, 0), 'prob': 2/3},
        'D': {'intermediate': (2, 1, -1), 'prob': 1/4},
    }
    
    answer_to_check = options[llm_answer_key]
    intermediate_state_to_check = answer_to_check['intermediate']
    prob_to_check = answer_to_check['prob']

    # --- Step 1: Verify selection rules for the proposed path ---
    def is_dipole_transition_allowed(state1, state2):
        """Checks if a transition between two states is allowed by dipole selection rules."""
        n1, l1, m1 = state1
        n2, l2, m2 = state2
        delta_l = l2 - l1
        delta_m = m2 - m1
        # For decay, n must decrease.
        if n2 >= n1:
            return False
        # Dipole selection rules
        if abs(delta_l) == 1 and abs(delta_m) <= 1:
            return True
        return False

    # Check the two transitions in the proposed path
    path_is_valid = (is_dipole_transition_allowed(initial_state, intermediate_state_to_check) and
                     is_dipole_transition_allowed(intermediate_state_to_check, final_state))

    if not path_is_valid:
        return f"Incorrect. The transition path in option {llm_answer_key}, {initial_state} -> {intermediate_state_to_check} -> {final_state}, is not allowed by the electric dipole selection rules."

    # --- Step 2: Calculate branching ratios ---
    # The probability of a specific route is the branching ratio of the first transition.
    # The rate is proportional to the square of the dipole matrix element.
    # By the Wigner-Eckart theorem, the angular part's dependence is given by the
    # square of a Clebsch-Gordan coefficient. The radial part is the same for all
    # m' sublevels of the intermediate state, so it cancels in the ratio.

    # Find all possible intermediate states for the first transition from |3,0,0>
    # Rule: n'=2, l'=1, m' in {-1, 0, 1}
    possible_intermediates = [(2, 1, -1), (2, 1, 0), (2, 1, 1)]
    
    # Calculate relative rates for the first transition to each possible intermediate state
    rates = {}
    l_i, m_i = initial_state[1], initial_state[2] # l=0, m=0

    for inter_state in possible_intermediates:
        l_f, m_f = inter_state[1], inter_state[2]
        q = m_f - m_i
        
        # CG(j1, m1, j2, m2, j3, m3) corresponds to <j1, m1; j2, m2 | j3, m3>
        # Here, j1=l_i, m1=m_i, j2=1 (dipole photon), m2=q, j3=l_f, m3=m_f
        cg_coeff = CG(j1=l_i, m1=m_i, j2=1, m2=q, j3=l_f, m3=m_f).doit()
        
        # Rate is proportional to the square of the CG coefficient
        rate = float(cg_coeff**2)
        rates[inter_state] = rate

    total_rate = sum(rates.values())

    if total_rate == 0:
        return "Incorrect. The calculation shows a total transition rate of zero from the initial state, which is physically incorrect for this scenario."

    # Calculate the branching ratio for the intermediate state in the answer
    calculated_prob = rates[intermediate_state_to_check] / total_rate

    # --- Step 3: Compare calculated probability with the answer's probability ---
    if not np.isclose(calculated_prob, prob_to_check, atol=1e-9):
        return (f"Incorrect. The probability for the path in option {llm_answer_key} is stated as {prob_to_check:.3f}, "
                f"but the calculated branching ratio is {calculated_prob:.3f}. "
                f"The relative rates for transitions from |3,0,0> to |2,1,-1>, |2,1,0>, and |2,1,1> are all equal (proportional to {rates[(2,1,-1)]}, {rates[(2,1,0)]}, and {rates[(2,1,1)]} respectively), "
                f"leading to a branching ratio of 1/3 for each path. Therefore, option A is correct, but the provided answer logic might be flawed if it selected a different option.")

    # If we reach here, the path is valid and the probability is correct.
    return "Correct"

# print(check_correctness_of_answer())