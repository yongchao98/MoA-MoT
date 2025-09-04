import re
import math

def check_astronomy_answer(llm_response_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the astronomy luminosity problem.

    The function performs the following steps:
    1.  Calculates the correct luminosity ratio from the physical principles given in the question.
    2.  Parses the LLM's response to find the chosen option (e.g., 'B') and the list of option values.
    3.  Determines which option is numerically closest to the correctly calculated answer.
    4.  Compares the LLM's chosen option with the correct option.
    5.  Returns "Correct" if they match, or an explanation of the error if they don't.
    """
    # --- Step 1: Calculate the correct answer from first principles ---

    # Given values from the problem statement
    radius_ratio_val = 1.5
    velocity_star2 = 700  # in km/s
    speed_of_light = 299792.458  # in km/s

    # The luminosity ratio L1/L2 is given by the Stefan-Boltzmann law:
    # L1/L2 = (R1/R2)^2 * (T1/T2)^4
    radius_term = radius_ratio_val**2

    # The temperature ratio T1/T2 must be derived from the Doppler effect and Wien's Law.
    # The problem states the *observed* peak wavelengths are the same: lambda_obs_1 = lambda_obs_2.
    # For Star 1 (v=0): lambda_obs_1 = lambda_emitted_1
    # For Star 2 (v=700): lambda_obs_2 = lambda_emitted_2 * (1 + v2/c)
    # So, lambda_emitted_1 = lambda_emitted_2 * (1 + v2/c)
    # From Wien's Law, T is proportional to 1/lambda_emitted, so T1/T2 = lambda_emitted_2 / lambda_emitted_1.
    # Therefore, T1/T2 = 1 / (1 + v2/c).
    temperature_ratio = 1 / (1 + velocity_star2 / speed_of_light)
    temperature_term = temperature_ratio**4

    # The final correct luminosity ratio
    correct_luminosity_ratio = radius_term * temperature_term

    # --- Step 2: Parse the LLM's response ---

    # Extract the final chosen option letter (e.g., 'B' from '<<<B>>>')
    final_answer_match = re.search(r'<<<([A-D])>>>', llm_response_text)
    if not final_answer_match:
        return "Incorrect: The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."
    
    chosen_option_letter = final_answer_match.group(1)

    # Extract the list of options and their values from the text.
    # This regex looks for the pattern "A) ~2.25", "B) ~2.23", etc. on separate lines.
    options_text_match = re.search(
        r'A\)\s*~?([\d\.]+).*?\n.*?'
        r'B\)\s*~?([\d\.]+).*?\n.*?'
        r'C\)\s*~?([\d\.]+).*?\n.*?'
        r'D\)\s*~?([\d\.]+)',
        llm_response_text,
        re.DOTALL | re.IGNORECASE
    )
    
    if not options_text_match:
        return "Incorrect: Could not parse the options A, B, C, D from the provided text to verify the answer."

    try:
        option_values = {
            'A': float(options_text_match.group(1)),
            'B': float(options_text_match.group(2)),
            'C': float(options_text_match.group(3)),
            'D': float(options_text_match.group(4)),
        }
    except (ValueError, IndexError):
        return "Incorrect: Failed to parse numerical values from the options list."

    chosen_option_value = option_values.get(chosen_option_letter)
    if chosen_option_value is None:
        return f"Incorrect: The chosen option '{chosen_option_letter}' is not a valid option (A, B, C, or D)."

    # --- Step 3: Verify the correctness of the chosen option ---

    # Find which of the given options is numerically closest to the calculated correct answer.
    closest_option_letter = min(option_values, key=lambda k: abs(option_values[k] - correct_luminosity_ratio))
    
    if chosen_option_letter == closest_option_letter:
        # The LLM chose the option that is indeed closest to the calculated value.
        # This confirms the answer is correct.
        return "Correct"
    else:
        # The LLM chose the wrong option.
        correct_option_value = option_values[closest_option_letter]
        return (
            f"Incorrect: The correct calculation for the luminosity ratio L1/L2 yields approximately {correct_luminosity_ratio:.4f}. "
            f"This value is derived from the formula L1/L2 = (R1/R2)^2 * (T1/T2)^4, where the temperature ratio T1/T2 = 1 / (1 + v2/c) due to the Doppler effect. "
            f"Among the given options, the closest value is {correct_option_value} (Option {closest_option_letter}). "
            f"The provided answer chose Option {chosen_option_letter} ({chosen_option_value}), which is not the closest and therefore incorrect."
        )

# The final response from the prompt to be checked
llm_response_to_check = """
Here is a step-by-step analysis of the question, highlighting the careful points required to arrive at the correct answer.

### Step-by-Step Analysis

1.  **Identify the Goal and the Primary Formula**: The question asks for the ratio of the luminosity of Star 1 to Star 2 (L₁ / L₂). For stars radiating as black bodies, the Stefan-Boltzmann Law gives the luminosity as L = 4πR²σT⁴. The ratio is therefore:
    L₁ / L₂ = (4πR₁²σT₁⁴) / (4πR₂²σT₂⁴)
    After canceling the constants (4, π, σ), the working formula is:
    **L₁ / L₂ = (R₁/R₂)² * (T₁/T₂)⁴**

2.  **Calculate the Radius Component**: This step is straightforward. The problem states that R₁ = 1.5 * R₂.
    Therefore, the squared ratio of the radii is:
    **(R₁/R₂)² = (1.5)² = 2.25**

3.  **Determine the Temperature Ratio (The Most Critical Step)**: This requires careful interpretation of the problem statement and the application of multiple physical principles.

    *   **Careful Point 1: "Appeared brightest" refers to Observed Wavelength.** The phrase "the wavelengths at which the stars appeared brightest" refers to the peak wavelength as *observed* from Earth (λ_obs), not the intrinsic wavelength *emitted* by the star (λ_emitted). This distinction is crucial because one of the stars is moving relative to the observer.

    *   **Careful Point 2: The Doppler Effect must be applied.** Star 2 has a radial velocity (v₂ = 700 km/s), which shifts its light's wavelength. Star 1 has zero radial velocity.
        *   For Star 1 (v₁ = 0): There is no Doppler shift. Therefore, λ_obs,1 = λ_emitted,1.
        *   For Star 2 (v₂ = 700 km/s): Its light is redshifted (wavelength is stretched). Therefore, λ_obs,2 = λ_emitted,2 * (1 + v₂/c).

    *   **Careful Point 3: Relate Temperature to Emitted Wavelength.** Wien's Displacement Law states that a star's temperature is inversely proportional to its *emitted* peak wavelength (T ∝ 1/λ_emitted). This means the temperature ratio is the inverse of the emitted wavelength ratio: **T₁/T₂ = λ_emitted,2 / λ_emitted,1**.

    *   **Careful Point 4: Combine the information to find the temperature ratio.** Since the problem states the observed wavelengths are the same (λ_obs,1 = λ_obs,2), we can set the expressions from the Doppler effect equal:
        λ_emitted,1 = λ_emitted,2 * (1 + v₂/c)
        Rearranging this gives the ratio we need for the temperature calculation:
        λ_emitted,2 / λ_emitted,1 = 1 / (1 + v₂/c)
        Therefore, the temperature ratio is:
        **T₁/T₂ = 1 / (1 + v₂/c)**

4.  **Identify Extraneous Information**: The problem provides the mass ratio of the stars (M₁ = 1.5 * M₂). For a black body radiation problem based on radius and temperature, the mass is irrelevant. This information is included to test a complete understanding of the governing principles.

5.  **Final Calculation**: Substitute the derived ratios for radius and temperature back into the main luminosity formula:
    L₁ / L₂ = (2.25) * [1 / (1 + v₂/c)]⁴
    L₁ / L₂ = 2.25 / (1 + v₂/c)⁴

    Now, plug in the numerical values (v₂ = 700 km/s and c ≈ 300,000 km/s):
    L₁ / L₂ = 2.25 / (1 + 700/300,000)⁴
    L₁ / L₂ = 2.25 / (1.002333...)⁴
    L₁ / L₂ ≈ 2.25 / 1.009355
    **L₁ / L₂ ≈ 2.229**

6.  **Conclusion**: The calculated luminosity ratio is approximately 2.229. Comparing this to the given options:
    *   A) ~2.25 (This is the incorrect answer obtained by ignoring the Doppler effect)
    *   B) ~2.23
    *   C) ~2.35
    *   D) ~2.32

The calculated value of 2.229 is closest to 2.23.

<<<B>>>
"""

# Run the check
result = check_astronomy_answer(llm_response_to_check)
print(result)