def check_decay_spectrum_answer(final_answer: str):
    """
    Checks the correctness of the answer to the nuclear decay spectrum question.

    The function verifies two key physical principles:
    1. Spectrum Type (Continuous vs. Discrete): A decay with more than two
       particles in the final state results in a continuous energy spectrum for
       any subset of those particles.
    2. Endpoint Change: The endpoint is determined by the Q-value of the reaction,
       which depends on the total rest mass of the final products. Less final
       rest mass means a higher Q-value and a higher endpoint.

    Args:
        final_answer: The letter ('A', 'B', 'C', or 'D') provided as the final answer.

    Returns:
        A string indicating "Correct" or the reason for the incorrectness.
    """

    # --- Step 1: Determine the correct nature of the spectrum ---
    # The variant decay is 2A -> 2B + 2E + M.
    # The final state has 3 light particles (2E, 1M) sharing the energy.
    # Rule: If the number of light final-state particles > 2, the spectrum is continuous.
    num_light_particles_variant = 3
    if num_light_particles_variant > 2:
        correct_spectrum_type = "continuous"
    else:
        correct_spectrum_type = "discrete"

    # --- Step 2: Determine the correct change in the endpoint ---
    # The endpoint is determined by the Q-value: Q = (m_initial - m_final)c^2
    # We only need to compare the rest mass of the particles that change.
    # Original decay has two 'V' particles, which are "light" (implying m_V > 0).
    # Variant decay has one 'M' particle, which is "massless" (m_M = 0).
    # Let's use placeholder values to represent this relationship.
    mass_V = 0.1  # Represents a small, non-zero mass
    mass_M = 0.0   # Represents a zero mass

    # The change in Q-value is proportional to the negative change in final rest mass.
    # Final rest mass of changing particles (Original): 2 * mass_V
    # Final rest mass of changing particles (Variant): mass_M
    if (2 * mass_V) > mass_M:
        # The final rest mass has DECREASED in the variant decay.
        # Therefore, the energy released (Q-value) has INCREASED.
        correct_endpoint_change = "increases"
    elif (2 * mass_V) < mass_M:
        correct_endpoint_change = "decreases"
    else:
        correct_endpoint_change = "remains the same"

    # --- Step 3: Map options to their claims and check the provided answer ---
    options = {
        'A': {'spectrum': 'continuous', 'endpoint': 'increases'},
        'B': {'spectrum': 'discrete', 'endpoint': 'increases'},
        'C': {'spectrum': 'continuous', 'endpoint': 'decreases'},
        'D': {'spectrum': 'discrete', 'endpoint': 'decreases'}
    }

    if final_answer not in options:
        return f"Invalid answer option '{final_answer}'. Please provide A, B, C, or D."

    provided_claims = options[final_answer]

    # Check against the derived correct properties
    if provided_claims['spectrum'] != correct_spectrum_type:
        return (f"Incorrect. The answer claims the spectrum becomes '{provided_claims['spectrum']}', "
                f"but it should be '{correct_spectrum_type}'. This is because the variant decay "
                f"still has {num_light_particles_variant} (>2) light particles in the final state, "
                f"which share the decay energy continuously.")

    if provided_claims['endpoint'] != correct_endpoint_change:
        return (f"Incorrect. The answer claims the endpoint '{provided_claims['endpoint']}', "
                f"but it actually '{correct_endpoint_change}'. This is because replacing two massive V particles "
                f"(mass = 2*m_V) with one massless M particle (mass = 0) reduces the total final rest mass, "
                f"thus increasing the total energy (Q-value) released.")

    return "Correct"

# The final answer provided by the LLM was 'A'.
llm_final_answer = 'A'

# Run the check
result = check_decay_spectrum_answer(llm_final_answer)
print(result)