import math

def check_annihilation_velocity():
    """
    Checks the correctness of the answer for the proton-antiproton annihilation problem.

    The function recalculates the velocity of particle A based on the principles
    of energy conservation in special relativity and compares it to the provided answer.
    """

    # --- Define problem constants and options ---
    # Rest mass energy of a proton in MeV (using CODATA 2018 value)
    m_p_c2 = 938.27208816
    # Rest mass energy of particle A in MeV
    m_A_c2 = 300.0
    
    # Options provided in the question
    options = {
        "A": 0.91,
        "B": 0.77,
        "C": 0.86,
        "D": 0.96
    }
    
    # The final answer from the LLM to be checked
    llm_answer_letter = "B"

    # --- Step 1: Calculate the initial energy (E_initial) ---
    # The initial system is a proton and an antiproton, assumed to be at rest ("slowly moving").
    # An antiproton has the same mass as a proton.
    E_initial = 2 * m_p_c2

    # --- Step 2: Set up the final energy equation (E_final) ---
    # The final system consists of 4 particles of type A (2 A+ and 2 A-).
    # By symmetry and conservation of momentum, they all have the same speed v,
    # and thus the same Lorentz factor, gamma.
    # E_final = 4 * gamma * m_A_c2

    # --- Step 3: Apply conservation of energy to find the Lorentz factor (gamma) ---
    # E_initial = E_final  =>  2 * m_p_c2 = 4 * gamma * m_A_c2
    try:
        gamma = E_initial / (4 * m_A_c2)
    except ZeroDivisionError:
        return "Incorrect: The mass of particle A cannot be zero."

    # --- Step 4: Calculate the velocity (beta = v/c) from gamma ---
    # The Lorentz factor is defined as gamma = 1 / sqrt(1 - beta^2).
    # Rearranging for beta gives: beta = sqrt(1 - 1/gamma^2).
    if gamma < 1:
        return f"Incorrect: Calculation resulted in a Lorentz factor (gamma) of {gamma:.4f}, which is less than 1 and physically impossible."
    
    beta_squared = 1 - (1 / (gamma**2))
    calculated_beta = math.sqrt(beta_squared)

    # --- Step 5: Check the correctness of the LLM's answer ---
    # Get the value corresponding to the LLM's chosen letter
    llm_answer_value = options.get(llm_answer_letter)
    if llm_answer_value is None:
        return f"Incorrect: The provided answer letter '{llm_answer_letter}' is not a valid option."

    # Check if the calculated value (beta) is close to the LLM's answer value.
    # A tolerance of 0.005 is appropriate for values given to two decimal places.
    if math.isclose(calculated_beta, llm_answer_value, abs_tol=0.005):
        return "Correct"
    else:
        # Find the closest correct option
        closest_option_letter = min(options, key=lambda k: abs(options[k] - calculated_beta))
        closest_option_value = options[closest_option_letter]
        
        return (f"Incorrect: The calculated velocity is v = {calculated_beta:.4f}c. "
                f"This value rounds to {round(calculated_beta, 2)}c, which corresponds to option {closest_option_letter} ({closest_option_value}c). "
                f"The provided answer was {llm_answer_letter} ({llm_answer_value}c).")

# Execute the check and print the result
result = check_annihilation_velocity()
if result == "Correct":
    print("Correct")
else:
    print(result)