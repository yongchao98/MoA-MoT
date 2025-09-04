def check_decay_spectrum_answer():
    """
    Checks the correctness of the answer to the nuclear decay spectrum question.

    The question compares two decays:
    1. Original: 2A -> 2B + 2E + 2V
    2. Variant:  2A -> 2B + 2E + M

    Constraints from the question:
    - V is a "much lighter particle", implying it has a non-zero rest mass (m_V > 0).
    - M is an "exotic, massless particle", meaning its rest mass is zero (m_M = 0).
    - The original spectrum is continuous.

    The provided answer to check is D:
    "The spectrum remains continuous with an adjusted shape, and the endpoint increases."
    """

    # Let's break down the claims of answer D
    claim_spectrum_is_continuous = True
    claim_endpoint_increases = True

    # --- Verification Step 1: Nature of the Spectrum (Continuous vs. Discrete) ---
    # A continuous energy spectrum for a subset of particles occurs when the total decay
    # energy is shared among 3 or more final-state particles. A discrete spectrum is
    # characteristic of a 2-body final state where kinematics are fixed.

    # In the variant decay (2A -> 2B + 2E + M), we analyze the number of final-state
    # particles that share the bulk of the kinetic energy. The heavy B nucleons have
    # negligible recoil. The energy is shared among the light particles: two E's and one M.
    num_light_final_particles = 2 + 1

    # Check if the spectrum should be continuous
    is_continuous = num_light_final_particles > 2

    if is_continuous != claim_spectrum_is_continuous:
        return (f"Incorrect: The answer claims the spectrum remains continuous. "
                f"However, the analysis of the number of final state particles is wrong. "
                f"The variant decay has {num_light_final_particles} light final-state particles. "
                f"A value > 2 implies a continuous spectrum, a value of 2 implies a discrete one. "
                f"The code's conclusion is that the spectrum is continuous, which matches the claim, but this check is in place for rigor.")

    # --- Verification Step 2: Change in the Endpoint ---
    # The endpoint of the spectrum is the maximum available kinetic energy, which is
    # determined by the Q-value of the reaction: Q = (m_initial - m_final) * c^2.
    # We only need to compare the change in the final rest mass.

    # Let's represent the masses symbolically. We only care about their properties (zero or non-zero).
    m_V_is_positive = True  # From "much lighter particle"
    m_M_is_zero = True      # From "massless"

    # The change in Q-value depends on the difference in the rest mass of the
    # particles that are different between the two decays: 2V vs. M.
    # Change in final rest mass energy = (Mass_final_variant - Mass_final_original) * c^2
    # Change = (m_M - 2*m_V) * c^2

    # Since m_M = 0 and m_V > 0, the change is (0 - 2*m_V) * c^2, which is negative.
    # This means the final rest mass of the variant decay is LESS than the original.
    final_mass_decreased = True

    # According to E=mc^2, if the final rest mass decreases, more mass is converted
    # into energy, so the total kinetic energy available (the Q-value/endpoint) must increase.
    endpoint_increases = final_mass_decreased

    if endpoint_increases != claim_endpoint_increases:
        return (f"Incorrect: The answer claims the endpoint increases. "
                f"This is based on the principle that a lower final rest mass yields a higher Q-value. "
                f"The code's analysis shows the final rest mass decreases (by 2*m_V), so the endpoint should increase. "
                f"The claim and the analysis do not match.")

    # --- Final Conclusion ---
    # Both claims made by answer D are verified by the physical principles outlined.
    return "Correct"

# Run the check and print the result
result = check_decay_spectrum_answer()
print(result)