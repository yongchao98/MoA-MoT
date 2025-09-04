import re

def check_correctness(question: str, final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the nuclear decay problem.

    The function performs a step-by-step logical check based on physics principles:
    1.  Determines if the energy spectrum should be continuous or discrete based on the number of final-state particles.
    2.  Determines if the energy endpoint should increase or decrease by analyzing the Q-value of the decay.
    3.  Compares the derived correct option with the provided final answer.

    Args:
        question: The original question string, which defines the options A, B, C, D.
        final_answer_text: The text containing the final answer in the format <<<X>>>.

    Returns:
        "Correct" if the answer is correct.
        A string explaining the reason for the error if the answer is incorrect.
    """
    # --- Step 0: Parse the provided answer and define options ---
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<X>>> in the provided text."
    
    provided_answer = match.group(1)

    # Define the options based on the question text
    options = {
        'A': {'spectrum': 'continuous', 'endpoint': 'decreases'},
        'B': {'spectrum': 'continuous', 'endpoint': 'increases'},
        'C': {'spectrum': 'discrete', 'endpoint': 'increases'},
        'D': {'spectrum': 'discrete', 'endpoint': 'decreases'}
    }

    # --- Step 1: Analyze the nature of the spectrum (Continuous vs. Discrete) ---
    # The variant decay is 2A -> 2B + 2E + M.
    # The final state has 5 particles (2B, 2E, 1M).
    # A decay with 3 or more final particles results in a continuous spectrum.
    num_final_particles_variant = 5
    if num_final_particles_variant > 2:
        correct_spectrum_type = "continuous"
    else:
        correct_spectrum_type = "discrete"

    # --- Step 2: Analyze the endpoint of the spectrum (Increase vs. Decrease) ---
    # The endpoint is determined by the Q-value: Q = (m_initial - m_final) * c^2.
    # A larger Q-value means a higher endpoint.
    # Q_orig corresponds to final mass: m(2B) + 2m(E) + 2m(V)
    # Q_var corresponds to final mass: m(2B) + 2m(E) + m(M)
    # We are given m(V) > 0 and m(M) = 0.
    # Therefore, the final mass in the variant decay is smaller than in the original decay.
    # A smaller final mass means a larger Q-value.
    # Q_var > Q_orig, so the endpoint increases.
    correct_endpoint_change = "increases"

    # --- Step 3: Determine the correct option based on the analysis ---
    derived_correct_option = None
    for option_key, properties in options.items():
        if properties['spectrum'] == correct_spectrum_type and properties['endpoint'] == correct_endpoint_change:
            derived_correct_option = option_key
            break
    
    if derived_correct_option is None:
        return "Internal logic error: Could not determine a correct option based on the physics."

    # --- Step 4: Compare derived correct option with the provided final answer ---
    if provided_answer == derived_correct_option:
        return "Correct"
    else:
        reason = (
            f"The final answer is {provided_answer}, but the correct answer should be {derived_correct_option}. "
            f"Reasoning: 1) The spectrum remains '{correct_spectrum_type}' because the variant decay is still a multi-body process (>2 particles in the final state). "
            f"2) The endpoint '{correct_endpoint_change}' because replacing two massive particles (2V) with one massless particle (M) reduces the total final rest mass, which increases the energy released (Q-value)."
        )
        return reason

# Extract the question and the final answer from the provided text
question_text = """
A) The spectrum remains continuous with an adjusted shape, and the endpoint decreases.
B) The spectrum remains continuous with an adjusted shape, and the endpoint increases.
C) The spectrum becomes discrete, and the endpoint increases.
D) The spectrum becomes discrete, and the endpoint decreases.
"""

final_answer_text_from_prompt = """
Here is a step-by-step analysis of the question, listing all careful points, followed by the final answer.

### Step-by-step Analysis

The problem asks us to compare the total energy spectrum of the 'E' particles in a new decay mode to an original one. To do this, we need to analyze two key features of the spectrum: its nature (continuous or discrete) and its endpoint (the maximum energy).

**1. Analyze the Nature of the Spectrum (Continuous vs. Discrete)**

*   **Careful Point:** The type of energy spectrum of a decay product is determined by the number of particles in the final state.
    *   A **two-body decay** (e.g., A → B + C) results in a **discrete** energy spectrum. Conservation of energy and momentum fully constrain the kinematics, so each product particle is emitted with a single, fixed energy.
    *   A **three-or-more-body decay** (e.g., A → B + C + D) results in a **continuous** energy spectrum. The total available energy is shared among the products in a continuous range of combinations, so any single product's energy is not fixed.

*   **Applying the Principle:**
    *   **Original Decay (`2A -> 2B + 2E + 2V`):** The final state has six particles. Ignoring the heavy, slow-recoiling B nucleons, we have four light particles (2E and 2V) sharing the decay energy. This is a multi-body decay, which is why the problem correctly states the spectrum is continuous.
    *   **Variant Decay (`2A -> 2B + 2E + M`):** The final state has five particles. Ignoring the B nucleons, we have three light particles (2E and M) sharing the decay energy. Since this is still a multi-body decay (more than two particles), the energy will continue to be shared in a variable way.

*   **Conclusion 1:** The energy spectrum of the E particles will **remain continuous**. This immediately eliminates options C and D, which state the spectrum becomes discrete.

**2. Analyze the Endpoint of the Spectrum (Increase vs. Decrease)**

*   **Careful Point:** The endpoint of the spectrum is the maximum possible kinetic energy for the particles of interest (the two E particles combined). This maximum energy is determined by the total energy released in the decay, known as the Q-value. The Q-value is the energy equivalent of the mass difference between the initial and final states: `Q = (m_initial - m_final) * c²`. A larger Q-value means more energy is available for the kinetic energy of the products, and thus a higher endpoint.

*   **Applying the Principle:**
    *   Let `m_A`, `m_B`, `m_E`, `m_V`, and `m_M` be the masses of the respective particles.
    *   **Q-value for Original Decay (`Q_orig`):**
        `Q_orig = (2m_A - (2m_B + 2m_E + 2m_V)) * c²`
    *   **Q-value for Variant Decay (`Q_var`):**
        `Q_var = (2m_A - (2m_B + 2m_E + m_M)) * c²`

*   **Comparing the Q-values:**
    *   We are given that particle M is **massless**, so `m_M = 0`.
    *   We are told V is a "much lighter particle," which in the context of nuclear decays implies it has a real, albeit small, positive mass. So, `m_V > 0`.
    *   The total rest mass of the final products is smaller in the variant decay because the mass of two V particles (`2m_V`, which is positive) is replaced by the mass of one M particle (`m_M`, which is zero).
    *   Since `m_final_variant < m_final_original`, more of the initial mass is converted into kinetic energy in the variant decay.
    *   Mathematically, `Q_var - Q_orig = (2m_V - m_M)c² = (2m_V - 0)c² = 2m_V c²`. Since `m_V > 0`, this difference is positive, meaning `Q_var > Q_orig`.

*   **Conclusion 2:** The total energy released in the variant decay is greater than in the original decay. Therefore, the maximum possible energy for the E particles, which is the endpoint of their spectrum, **increases**.

**Final Synthesis**

1.  The spectrum remains **continuous** because the decay is still a multi-body process.
2.  The shape of the spectrum is **adjusted** because the decay kinematics and available phase space are different for a 3-light-particle final state compared to a 4-light-particle one.
3.  The endpoint **increases** because replacing two particles with mass (2V) with a single massless particle (M) reduces the total rest mass of the products, thereby increasing the total energy released as kinetic energy.

This combination of findings corresponds to option B.

<<<B>>>
"""

# Run the check
result = check_correctness(question_text, final_answer_text_from_prompt)
print(result)