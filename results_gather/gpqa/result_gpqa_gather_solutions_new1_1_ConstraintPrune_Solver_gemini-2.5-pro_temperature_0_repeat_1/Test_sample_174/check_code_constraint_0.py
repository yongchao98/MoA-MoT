import math
import re

def check_physics_answer(final_answer_text):
    """
    Checks the correctness of the answer for the oscillating charge distribution problem.

    The function verifies two physical principles:
    1. The fraction of maximum power at 30 degrees for dipole radiation.
    2. The wavelength dependence of dipole radiation power.
    """
    # --- Step 1: Define the options and the correct physical principles ---

    # The options provided in the question
    options = {
        'A': {'fraction': 1/2, 'lambda_power': -4},
        'B': {'fraction': 1/4, 'lambda_power': -4},
        'C': {'fraction': 1/4, 'lambda_power': -3},
        'D': {'fraction': 3/4, 'lambda_power': -6}
    }

    # Principle 1: Angular dependence of power for an electric dipole
    # Power is proportional to sin^2(theta). Max power is at theta=90 deg.
    # Fraction at 30 deg is sin^2(30) / sin^2(90) = (0.5)^2 / 1^2 = 0.25
    theta_rad = math.radians(30)
    expected_fraction = math.sin(theta_rad)**2

    # Principle 2: Wavelength dependence of power for an electric dipole
    # Power is proportional to omega^4, and omega is proportional to 1/lambda.
    # Therefore, power is proportional to (1/lambda)^4 = lambda^-4.
    expected_lambda_power = -4

    # --- Step 2: Parse the provided final answer ---
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Incorrect: The final answer is not in the required format '<<<X>>>'."

    chosen_option_key = match.group(1)
    if chosen_option_key not in options:
        return f"Incorrect: The chosen option '{chosen_option_key}' is not a valid option (A, B, C, D)."

    chosen_option_values = options[chosen_option_key]
    chosen_fraction = chosen_option_values['fraction']
    chosen_lambda_power = chosen_option_values['lambda_power']

    # --- Step 3: Compare the chosen answer with the correct values ---

    # Check the fraction part
    if not math.isclose(chosen_fraction, expected_fraction):
        return (f"Incorrect: The fraction of power is wrong. "
                f"The correct fraction at 30 degrees should be sin^2(30) = {expected_fraction:.2f}. "
                f"The chosen option '{chosen_option_key}' has a fraction of {chosen_fraction}.")

    # Check the wavelength dependence part
    if chosen_lambda_power != expected_lambda_power:
        return (f"Incorrect: The wavelength dependence is wrong. "
                f"The power should be proportional to lambda^{expected_lambda_power}. "
                f"The chosen option '{chosen_option_key}' has a dependence of lambda^{chosen_lambda_power}.")

    # If both checks pass, the answer is correct.
    return "Correct"

# The final answer provided by the LLM
final_answer = """
Here is a step-by-step analysis of the problem, incorporating the insights from the candidate answers.

### Step-by-Step Analysis

1.  **Identify the Physical Model:** The problem describes an "oscillating charge distribution" with a "spheroid in shape" and a symmetry axis along the z-axis. In electromagnetism, the simplest and most dominant form of radiation from a compact, oscillating source is electric dipole radiation. The symmetry described suggests the effective electric dipole moment oscillates along the z-axis. This is the standard model for such a problem unless specific conditions (like a vanishing dipole moment) force consideration of higher-order multipoles (e.g., magnetic dipole or electric quadrupole).

2.  **Analyze the Angular Dependence of Radiated Power:**
    For an electric dipole oscillating along the z-axis, the time-averaged power radiated per unit solid angle (`dP/dΩ`) is not uniform. It follows a specific angular distribution given by:
    `dP/dΩ ∝ sin²(θ)`
    where `θ` is the angle measured from the axis of oscillation (the z-axis).

    *   The problem states that the maximum power radiated per unit solid angle is `A`. This maximum occurs where `sin²(θ)` is at its maximum value.
    *   The maximum value of `sin²(θ)` is 1, which occurs at `θ = 90°` (in the plane perpendicular to the oscillation axis).
    *   Therefore, we can express the power at any angle `θ` as `f(θ) = A * sin²(θ)`.
    *   We need to find the fraction of `A` that is radiated at `θ = 30°`. We calculate `f(30°)`:
        `f(30°) = A * sin²(30°)`
    *   Using the standard trigonometric value `sin(30°) = 1/2`, we get:
        `f(30°) = A * (1/2)² = A * (1/4)`
    *   The fraction of `A` radiated at this angle is **1/4**.

3.  **Analyze the Wavelength (λ) Dependence of Radiated Power:**
    The power radiated by an oscillating electric dipole is also strongly dependent on the frequency of oscillation. The full formula shows that the power is proportional to the fourth power of the angular frequency (`ω`):
    `dP/dΩ ∝ ω⁴`

    *   The angular frequency `ω` is related to the wavelength `λ` by the formula `ω = 2πc/λ`, where `c` is the speed of light. This means `ω` is inversely proportional to `λ`.
    *   Substituting this relationship into the power formula:
        `dP/dΩ ∝ (1/λ)⁴ = λ⁻⁴`
    *   Therefore, a possible form for the wavelength dependence of the function `f` is **λ⁻⁴**. This is a characteristic feature of electric dipole radiation and is also responsible for phenomena like Rayleigh scattering (which explains why the sky is blue).

4.  **Combine Results and Select the Correct Option:**
    *   The fraction of maximum power at `θ = 30°` is **1/4**.
    *   The wavelength dependence of the power function `f` is **λ⁻⁴**.

    Now, we match this result with the given options:
    A) 1/2, λ⁻⁴
    B) 1/4, λ⁻⁴
    C) 1/4, λ⁻³
    D) 3/4, λ⁻⁶

    The only option that matches both of our findings is **B**.

### Careful Points to Note

*   **Model Validity:** The electric dipole model is the correct leading-order approximation for this system. Higher-order radiation types like electric quadrupole have different dependencies (e.g., power `∝ λ⁻⁶`), but they are typically much weaker and only become dominant if the electric dipole radiation is forbidden by symmetry, which is not the case here.
*   **Squaring the Sine:** A common mistake is to use `sin(30°) = 1/2` instead of `sin²(30°) = 1/4`. This is a crucial step to get the correct fraction.
*   **Frequency vs. Wavelength:** Another critical point is the conversion from frequency dependence (`ω⁴`) to wavelength dependence (`λ⁻⁴`). Forgetting the inverse relationship between `ω` and `λ` leads to an incorrect exponent.
*   **Matching:** As seen in the candidate answers, even with the correct physics, a final error can be made by incorrectly matching the derived result `(1/4, λ⁻⁴)` to the lettered options.

The analysis consistently and robustly points to a single correct answer.

<<<B>>>
"""

# Run the check
result = check_physics_answer(final_answer)
print(result)