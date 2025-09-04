import numpy as np
from sympy.physics.quantum.cg import CG

def check_angular_momentum_coupling():
    """
    Checks the correctness of the answer to the quantum mechanics problem.

    The problem asks for the joint probability of measuring m1=-1 and m2=-1
    for a system in the coupled state |l=2, m=-1>, which is formed by
    coupling two p-orbitals (l1=1, l2=1).

    The probability is the square of the Clebsch-Gordan coefficient:
    P = |<l1, m1; l2, m2 | l, m>|^2
    """
    
    # --- Problem Parameters ---
    # Initial coupled state: |l1, l2, l, m> = |1, 1, 2, -1>
    l1 = 1
    l2 = 1
    l = 2
    m = -1

    # Desired measurement outcome (uncoupled state)
    # Eigenvalue of L_1z is -ħ -> m1 = -1
    # Eigenvalue of L_2z is -ħ -> m2 = -1
    m1_measured = -1
    m2_measured = -1

    # The final answer provided by the LLM to be checked
    # Options: A) 1/2, B) 1, C) 2/3, D) 0
    llm_answer_choice = 'D'
    options = {'A': 1/2, 'B': 1, 'C': 2/3, 'D': 0}
    llm_answer_value = options[llm_answer_choice]

    # --- Verification using the Selection Rule ---
    # A fundamental selection rule for Clebsch-Gordan coefficients is that
    # they are zero unless m = m1 + m2.
    
    m_sum_measured = m1_measured + m2_measured
    
    if m != m_sum_measured:
        # The rule is violated, so the probability must be 0.
        calculated_probability = 0
    else:
        # If the rule were satisfied, we would calculate the probability using
        # the square of the Clebsch-Gordan coefficient.
        # This branch is not strictly necessary for this problem but shows the general method.
        cg_coefficient = CG(j1=l1, m1=m1_measured, j2=l2, m2=m2_measured, j=l, m=m).doit()
        calculated_probability = float(cg_coefficient)**2

    # --- Final Check ---
    # Compare the calculated probability with the LLM's answer.
    if np.isclose(calculated_probability, llm_answer_value):
        return "Correct"
    else:
        reason = (
            f"The provided answer is {llm_answer_choice}, which corresponds to a probability of {llm_answer_value}.\n"
            f"However, the correct probability is {calculated_probability}.\n"
            f"Reason: The calculation is based on the conservation of the z-component of angular momentum, which requires m = m1 + m2.\n"
            f"In this case, the initial state has m = {m}.\n"
            f"The desired measurement outcome corresponds to m1 = {m1_measured} and m2 = {m2_measured}, giving a sum m1 + m2 = {m_sum_measured}.\n"
            f"Since m ({m}) is not equal to m1 + m2 ({m_sum_measured}), the selection rule is violated, and the probability of this outcome is exactly 0."
        )
        return reason

# Execute the check
result = check_angular_momentum_coupling()
print(result)