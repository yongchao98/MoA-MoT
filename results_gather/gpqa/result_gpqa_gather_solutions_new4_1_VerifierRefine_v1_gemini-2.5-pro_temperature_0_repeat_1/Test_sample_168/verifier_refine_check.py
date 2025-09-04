import re

def check_correctness_of_nuclear_decay_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to the nuclear decay problem.

    The function verifies two key physical principles:
    1. The nature of the energy spectrum (continuous vs. discrete).
    2. The change in the energy endpoint (increase vs. decrease).
    """

    # --- Step 1: Define the correct physical conclusions based on the question ---

    # Principle 1: Nature of the spectrum.
    # A decay with >2 final state particles has a continuous spectrum.
    # Original decay (2A -> 2B + 2E + 2V): 6 particles -> continuous.
    # Variant decay (2A -> 2B + 2E + M): 5 particles -> continuous.
    # Therefore, the spectrum must remain continuous.
    correct_spectrum_nature = "continuous"

    # Principle 2: Endpoint of the spectrum.
    # Endpoint is determined by the Q-value: Q = (m_initial - m_final)c^2.
    # Q is maximized when m_final is minimized.
    # We compare the mass of the light final particles: (2*m_V) vs. (m_M).
    # Given: m_M = 0 and m_V > 0.
    # The total final rest mass is lower in the variant decay.
    # This means Q_variant > Q_original.
    # Therefore, the endpoint must increase.
    correct_endpoint_change = "increases"

    # --- Step 2: Parse the LLM's provided answer ---

    # Find the final answer choice, e.g., <<<B>>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."

    llm_choice = match.group(1)

    # Map the option choices to their physical claims
    options = {
        'A': ("discrete", "decreases"),
        'B': ("continuous", "increases"),
        'C': ("discrete", "increases"),
        'D': ("continuous", "decreases")
    }

    if llm_choice not in options:
        return f"Error: Invalid option '{llm_choice}' found."

    llm_spectrum_nature, llm_endpoint_change = options[llm_choice]

    # --- Step 3: Compare the LLM's claims to the correct conclusions ---

    errors = []

    # Check the first condition: the nature of the spectrum.
    if llm_spectrum_nature != correct_spectrum_nature:
        error_msg = (
            f"The nature of the spectrum is incorrect. "
            f"The answer claims the spectrum becomes '{llm_spectrum_nature}', "
            f"but it should remain '{correct_spectrum_nature}'. "
            f"This is because the variant decay (2A -> 2B + 2E + M) still results in a multi-body final state (5 particles), "
            f"so the energy is shared continuously."
        )
        errors.append(error_msg)

    # Check the second condition: the change in the endpoint.
    if llm_endpoint_change != correct_endpoint_change:
        error_msg = (
            f"The change in the endpoint is incorrect. "
            f"The answer claims the endpoint '{llm_endpoint_change}', "
            f"but it should '{correct_endpoint_change}'. "
            f"This is because the Q-value of the reaction increases. "
            f"The final mass is lower in the variant decay because two massive particles (2V) are replaced by one massless particle (M). "
            f"A lower final mass means a higher Q-value and thus a higher energy endpoint."
        )
        errors.append(error_msg)

    if not errors:
        return "Correct"
    else:
        # Join all found errors into a single message.
        return "\n".join(errors)

# The provided answer from the LLM to be checked.
llm_answer_to_check = """
Here is a step-by-step analysis of the problem to determine the final answer.

### Step 1: Analyze the Nature of the Energy Spectrum (Continuous vs. Discrete)

The type of energy spectrum of a decay product is determined by the number of particles in the final state.
*   **Two-body decay:** If a decay results in only two final particles, conservation of energy and momentum fully constrain the kinematics. Each particle is emitted with a specific, fixed energy, resulting in a **discrete** energy spectrum (a sharp peak).
*   **Three-or-more-body decay:** If a decay results in three or more final particles, the available energy can be shared among them in a continuous range of ways. This results in a **continuous** energy spectrum for each particle.

Let's apply this to the two processes in the question:
1.  **Original Decay (`2A -> 2B + 2E + 2V`):** The final state consists of six particles (two B, two E, two V). Since this is a multi-body decay, the energy is shared, and the spectrum of the E particles is continuous, as stated in the problem.
2.  **Variant Decay (`2A -> 2B + 2E + M`):** The final state now consists of five particles (two B, two E, one M). This is still a multi-body decay. The energy released is shared between the E particles and the M particle. Therefore, the energy spectrum of the E particles will **remain continuous**.

This conclusion eliminates options A and C, which incorrectly state that the spectrum becomes discrete.

### Step 2: Analyze the Endpoint of the Spectrum

The endpoint of the energy spectrum represents the maximum possible kinetic energy that the particles of interest (the two E particles combined) can have. This maximum energy is determined by the total energy released in the decay, known as the Q-value. The Q-value is the difference between the total rest mass of the initial particles and the total rest mass of the final particles, converted to energy via `E = mc²`.

Let's calculate the Q-value for each decay, which corresponds to the total kinetic energy released.
*   **Q-value for the Original Decay (`Q_original`):**
    The total mass converted to kinetic energy is the initial mass minus the final mass.
    `Q_original = (Mass_initial - Mass_final_1) * c²`
    `Q_original = (2m_A - (2m_B + 2m_E + 2m_V)) * c²`

*   **Q-value for the Variant Decay (`Q_variant`):**
    `Q_variant = (Mass_initial - Mass_final_2) * c²`
    `Q_variant = (2m_A - (2m_B + 2m_E + m_M)) * c²`

Now, we use the information given about the particles:
*   The particle M is **massless**, so `m_M = 0`.
*   The particle V is a "much lighter particle," which implies it has a real, non-zero mass (`m_V > 0`), just like a neutrino.

Let's compare the two Q-values:
`Q_variant = (2m_A - 2m_B - 2m_E) * c²`
`Q_original = (2m_A - 2m_B - 2m_E) * c² - 2m_V * c²`

The difference is `Q_variant - Q_original = 2m_V * c²`. Since `m_V > 0`, this difference is a positive value. This means `Q_variant > Q_original`.

The total energy released in the variant decay is greater than in the original decay. This extra energy comes from the fact that no energy is needed to create the rest mass of the V particles. This larger amount of available kinetic energy means the maximum possible energy for the E particles is higher. Therefore, the endpoint of the spectrum **increases**.

### Step 3: Final Conclusion

1.  The spectrum **remains continuous** because the variant decay is still a multi-body process.
2.  The shape of the spectrum is **adjusted** because the phase space of a 5-body final state is different from that of a 6-body final state.
3.  The endpoint of the spectrum **increases** because replacing two massive V particles with one massless M particle reduces the total rest mass of the products, thus increasing the total energy released as kinetic energy.

This combination of findings matches option B.

<<<B>>>
"""

# Run the check
result = check_correctness_of_nuclear_decay_answer(llm_answer_to_check)
print(result)