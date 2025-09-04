import sympy
from sympy.physics.quantum.cg import CG
import math

def check_correctness():
    """
    Checks the correctness of the LLM's answer by calculating the exact probability
    using Clebsch-Gordan coefficients.
    """
    # --- Parameters from the question ---
    # Two electrons in p-orbitals -> l1=1, l2=1
    l1 = 1
    l2 = 1
    # The system is in the coupled state |l, m> = |2, -1>
    l_total = 2
    m_total = -1

    # --- Measurement to be checked ---
    # We want the probability of measuring L_1z eigenvalue as -ħ and L_2z as -ħ.
    # This corresponds to measuring m1 = -1 and m2 = -1.
    m1_measured = -1
    m2_measured = -1

    # --- LLM's Answer ---
    # The provided answer is 'A', which corresponds to a probability of 0.
    llm_answer_prob = 0.0

    # --- Verification using the selection rule ---
    # The Clebsch-Gordan coefficient <l1, m1; l2, m2 | l, m> is non-zero
    # only if m = m1 + m2. The probability is the square of this coefficient.
    # Let's check this fundamental condition first, as the LLM did.
    m_sum_of_measured_states = m1_measured + m2_measured

    if m_sum_of_measured_states != m_total:
        # The selection rule is violated. The probability must be 0.
        # This matches the reasoning and conclusion of the LLM.
        # We can confirm this with a direct calculation.
        calculated_prob = 0.0
    else:
        # The selection rule is satisfied. We must calculate the non-zero probability
        # using the Clebsch-Gordan coefficient.
        try:
            cg_coefficient = CG(l1, m1_measured, l2, m2_measured, l_total, m_total).doit()
            # The probability is the square of the coefficient's value.
            calculated_prob = float(cg_coefficient**2)
        except Exception as e:
            return f"An error occurred during the sympy calculation: {e}"

    # --- Final Check ---
    # Compare the calculated probability with the LLM's answer.
    if math.isclose(calculated_prob, llm_answer_prob, abs_tol=1e-9):
        return "Correct"
    else:
        return (f"Incorrect. The question asks for the probability of measuring m1={m1_measured} and m2={m2_measured} "
                f"from the state |l={l_total}, m={m_total}>. The selection rule m = m1 + m2 dictates the outcome. "
                f"Here, m1+m2 = {m_sum_of_measured_states} and m = {m_total}. Since these are not equal, the probability must be 0. "
                f"The calculated probability is {calculated_prob}, but the LLM's answer corresponds to {llm_answer_prob}.")

# Run the check
result = check_correctness()
print(result)