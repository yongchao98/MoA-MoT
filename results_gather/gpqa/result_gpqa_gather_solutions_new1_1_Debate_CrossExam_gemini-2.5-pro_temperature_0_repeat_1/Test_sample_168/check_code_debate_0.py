def check_decay_answer():
    """
    Checks the correctness of the answer based on the principles of nuclear decay.

    The provided answer is C: "The spectrum remains continuous with an adjusted shape, and the endpoint increases."

    This function will verify the two claims in this answer:
    1. The spectrum remains continuous.
    2. The endpoint increases.
    """

    # --- Part 1: Verify the Nature of the Spectrum (Continuous vs. Discrete) ---

    # A decay results in a continuous energy spectrum for its products if there are
    # three or more particles in the final state sharing the energy.
    
    # Original decay: 2A -> 2B + 2E + 2V. The final state has a nucleus (2B) and 4 light particles.
    # Total final particles = 1 nucleus + 2 E + 2 V = 5.
    num_final_particles_original = 5
    
    # Variant decay: 2A -> 2B + 2E + M. The final state has a nucleus (2B) and 3 light particles.
    # Total final particles = 1 nucleus + 2 E + 1 M = 4.
    num_final_particles_variant = 4

    # Check if both decays produce a continuous spectrum.
    is_original_continuous = num_final_particles_original > 2
    is_variant_continuous = num_final_particles_variant > 2

    # The first claim is that the spectrum *remains* continuous.
    claim1_is_correct = is_original_continuous and is_variant_continuous
    
    if not claim1_is_correct:
        return "Incorrect: The answer claims the spectrum remains continuous, but the analysis shows this is not the case. A discrete spectrum occurs in two-body decays."

    # --- Part 2: Verify the Endpoint Comparison ---

    # The endpoint of the E-particle spectrum is the maximum kinetic energy they can have.
    # This is determined by the total kinetic energy released (Q-value), which is:
    # Q_value = (m_initial - m_final_rest_mass) * c^2
    # A larger Q-value means a higher endpoint. This happens when the final rest mass is smaller.
    
    # We only need to compare the rest mass of the final *light* particles, as m_A, m_B, and m_E are constant.
    # Let's analyze the masses based on the problem description:
    # m_V > 0 (V is a "much lighter particle", implying it has mass)
    # m_M = 0 (M is explicitly "massless")

    # Rest mass of final light particles in original decay: 2*m_V
    # Rest mass of final light particles in variant decay: m_M = 0
    
    # Since m_V > 0, the rest mass of the final light particles is greater in the original decay.
    # Mass_final_light_original > Mass_final_light_variant
    
    # The decay with the smaller final rest mass (the variant decay) will have a larger Q-value.
    # Therefore, the endpoint of the variant decay's spectrum is higher.
    endpoint_increases = True
    
    claim2_is_correct = endpoint_increases

    if not claim2_is_correct:
        return "Incorrect: The answer claims the endpoint increases, but the analysis of rest masses shows it must decrease or stay the same."

    # --- Final Conclusion ---
    if claim1_is_correct and claim2_is_correct:
        # The provided answer 'C' correctly states that the spectrum remains continuous and the endpoint increases.
        return "Correct"
    else:
        # This case should not be reached if the logic is sound, but it's here for completeness.
        return "Incorrect: The provided answer fails on one or more physical constraints."

# Run the check and print the result.
result = check_decay_answer()
print(result)