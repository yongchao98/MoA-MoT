def check_answer_correctness():
    """
    Checks the correctness of the final answer for the nuclear decay problem.

    The function breaks the problem down into two verifiable constraints:
    1. The nature of the energy spectrum (continuous vs. discrete).
    2. The change in the endpoint energy (increases vs. decreases).

    It then evaluates the provided final answer against these physical constraints.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'B'

    # --- Constraint 1: Nature of the Energy Spectrum ---
    # A continuous spectrum results from a decay with 3 or more final-state particles
    # sharing the energy. A discrete spectrum results from a 2-body decay.
    # Original decay (2A -> 2B + 2E + 2V): 4 light final-state particles (2E, 2V).
    # Variant decay (2A -> 2B + 2E + M): 3 light final-state particles (2E, M).
    # Since the number of light particles is > 2 in the variant decay, the spectrum remains continuous.
    correct_spectrum_nature = "continuous"

    # --- Constraint 2: Change in the Endpoint Energy ---
    # The endpoint is the maximum available kinetic energy, determined by Q = (m_initial - m_final)c^2.
    # Let Q_nuc = (m(2A) - m(2B))c^2.
    # Endpoint_original = Q_nuc - (2*m_E + 2*m_V)c^2
    # Endpoint_variant = Q_nuc - (2*m_E + m_M)c^2
    # Given: m_M = 0 (massless) and m_V > 0 ("light particle", not massless).
    # Therefore, Endpoint_variant = Q_nuc - 2*m_E*c^2
    # And Endpoint_original = Endpoint_variant - 2*m_V*c^2
    # Since 2*m_V*c^2 > 0, Endpoint_variant > Endpoint_original. The endpoint increases.
    correct_endpoint_change = "increases"

    # Define the properties of each multiple-choice option
    options = {
        'A': {'spectrum': 'discrete', 'endpoint': 'decreases'},
        'B': {'spectrum': 'continuous', 'endpoint': 'increases'},
        'C': {'spectrum': 'continuous', 'endpoint': 'decreases'},
        'D': {'spectrum': 'discrete', 'endpoint': 'increases'}
    }

    # Check if the provided answer is a valid option
    if llm_final_answer not in options:
        return f"Invalid Answer: The final answer '{llm_final_answer}' is not one of the possible options (A, B, C, D)."

    # Get the properties of the provided answer
    answer_properties = options[llm_final_answer]

    # Verify against Constraint 1
    if answer_properties['spectrum'] != correct_spectrum_nature:
        return (f"Incorrect. The answer claims the spectrum becomes '{answer_properties['spectrum']}', "
                f"but the correct physical principle dictates it should remain '{correct_spectrum_nature}'. "
                f"This is because the variant decay (2A -> 2B + 2E + M) still has three light particles in the final state sharing the energy, "
                f"which results in a continuous spectrum.")

    # Verify against Constraint 2
    if answer_properties['endpoint'] != correct_endpoint_change:
        return (f"Incorrect. The answer claims the endpoint '{answer_properties['endpoint']}', "
                f"but the correct physical principle dictates it '{correct_endpoint_change}'. "
                f"This is because replacing two massive particles (2V) with one massless particle (M) reduces the total rest mass of the final products. "
                f"By E=mc^2, this mass difference is converted into additional kinetic energy, increasing the maximum possible energy (the endpoint) for the E particles.")

    # If both constraints are satisfied
    return "Correct"

# Run the check and print the result
result = check_answer_correctness()
print(result)