import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided LLM answer for the physics problem.
    It recalculates the velocity of particle A and compares it to the given answer.
    """

    # --- Problem Constants and Given Data ---
    # Rest mass energy of a proton (CODATA 2018 value) in MeV
    m_p_c2 = 938.272
    # Rest mass energy of particle A in MeV
    m_A_c2 = 300
    # Number of A particles produced
    num_A_particles = 4

    # --- Physics Calculation ---
    # The problem is solved using the conservation of energy.
    # The initial energy is the rest mass energy of the proton-antiproton pair at rest.
    # E_initial = 2 * m_p_c2
    # The final energy is the total energy of the four 'A' particles.
    # E_final = num_A_particles * gamma * m_A_c2
    # By setting E_initial = E_final, we can solve for the Lorentz factor, gamma.

    try:
        # Calculate the Lorentz factor (gamma)
        gamma = (2 * m_p_c2) / (num_A_particles * m_A_c2)
    except ZeroDivisionError:
        return "Calculation Error: Division by zero. The mass of particle A or the number of particles cannot be zero."

    # The Lorentz factor must be >= 1 for a physical velocity.
    if gamma < 1:
        return f"Calculation Error: The calculated Lorentz factor gamma is {gamma:.4f}, which is less than 1 and physically impossible."

    # Calculate the velocity (beta = v/c) from the Lorentz factor.
    # gamma = 1 / sqrt(1 - beta^2)  =>  beta = sqrt(1 - 1/gamma^2)
    calculated_beta = math.sqrt(1 - 1 / gamma**2)

    # --- Verification of the LLM's Answer ---
    # The LLM's final answer is <<<B>>>.
    # The options provided in the LLM's response are:
    # A) 0.96c, B) 0.77c, C) 0.91c, D) 0.86c
    llm_answer_letter = 'B'
    
    options = {
        'A': 0.96,
        'B': 0.77,
        'C': 0.91,
        'D': 0.86
    }

    # Get the numerical value corresponding to the LLM's answer letter.
    llm_answer_value = options.get(llm_answer_letter)
    if llm_answer_value is None:
        return f"Invalid Answer Format: The answer letter '{llm_answer_letter}' is not one of the valid options (A, B, C, D)."

    # Compare the calculated result with the LLM's answer.
    # We check if the calculated value rounds to the answer's value.
    # A tolerance of 0.005 is appropriate for checking rounding to two decimal places.
    tolerance = 0.005
    if abs(calculated_beta - llm_answer_value) < tolerance:
        return "Correct"
    else:
        # If incorrect, provide a detailed reason.
        reason = (
            f"Incorrect. The calculation shows the velocity should be approximately {calculated_beta:.4f}c. "
            f"This value rounds to {round(calculated_beta, 2)}c. "
            f"The provided answer is '{llm_answer_letter}', which corresponds to a velocity of {llm_answer_value}c. "
            f"The calculated result does not match the answer."
        )
        return reason

# print(check_correctness_of_answer())