import math
import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the physics problem.

    The problem involves electric dipole radiation, where:
    1. Power is proportional to sin^2(theta).
    2. Power is proportional to lambda^-4.
    """

    # --- Step 1: Define the options from the question ---
    options = {
        'A': {'fraction': 3/4, 'exponent': -6},
        'B': {'fraction': 1/4, 'exponent': -3},
        'C': {'fraction': 1/4, 'exponent': -4},
        'D': {'fraction': 1/2, 'exponent': -4}
    }

    # --- Step 2: Calculate the correct physical values ---
    # 1. Calculate the fraction of maximum power at 30 degrees.
    # Power is proportional to sin^2(theta). Max power is at theta=90 deg.
    # Fraction = sin^2(30) / sin^2(90) = (0.5)^2 / 1^2 = 0.25
    theta_rad = math.radians(30)
    expected_fraction = math.sin(theta_rad)**2

    # 2. Determine the wavelength dependence exponent.
    # Power is proportional to omega^4, and omega is proportional to 1/lambda.
    # So, Power is proportional to (1/lambda)^4 = lambda^-4.
    expected_exponent = -4

    # --- Step 3: Parse the LLM's chosen answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer is not in the required format '<<<X>>>'."
    
    chosen_option_key = match.group(1)
    if chosen_option_key not in options:
        return f"Incorrect: The chosen option '{chosen_option_key}' is not a valid option (A, B, C, or D)."

    chosen_option = options[chosen_option_key]
    llm_fraction = chosen_option['fraction']
    llm_exponent = chosen_option['exponent']

    # --- Step 4: Compare the LLM's answer with the correct values ---
    fraction_correct = math.isclose(llm_fraction, expected_fraction)
    exponent_correct = (llm_exponent == expected_exponent)

    if fraction_correct and exponent_correct:
        return "Correct"
    else:
        error_messages = []
        if not fraction_correct:
            error_messages.append(
                f"The fraction of maximum power is incorrect. "
                f"For an electric dipole, the power is proportional to sin^2(theta). "
                f"At theta=30 degrees, the fraction should be sin^2(30) = 0.25 (or 1/4). "
                f"The chosen option '{chosen_option_key}' has a fraction of {llm_fraction}."
            )
        if not exponent_correct:
            error_messages.append(
                f"The wavelength dependence is incorrect. "
                f"For electric dipole radiation, power is proportional to omega^4, which is lambda^-4. "
                f"The expected exponent is -4. "
                f"The chosen option '{chosen_option_key}' has an exponent of {llm_exponent}."
            )
        return "Incorrect: " + " ".join(error_messages)

# Provided question and LLM answer
question = """
Consider an oscillating charge distribution, which is spheroid in shape (with the symmetry axis along z-axis). Because it oscillates, it would radiate at a wavelength \lambda. In the radiation zone, the radiated power per unit solid angle will be a function f of \lambda and the angle \theta (with z-axis), along which the power is measured. If the maximum power thus radiated is A, the fraction of A that will be radiated at an angle \theta = 30^0 and a possible form of f are, respectively,

A) 3/4, \lambda^(-6)
B) 1/4, \lambda^(-3)
C) 1/4, \lambda^(-4)
D) 1/2, \lambda^(-4)
"""

llm_answer = """
Here is a step-by-step analysis of the problem, synthesizing the reasoning from the provided candidate answers and identifying the correct physical principles.

### Step-by-Step Analysis

1.  **Identify the Physical Model:** The problem describes an "oscillating charge distribution" with a "spheroid in shape" and a symmetry axis along the z-axis. For a compact oscillating source (where the size is much smaller than the radiated wavelength λ), the most dominant form of radiation is electric dipole radiation. The symmetry described implies the electric dipole moment oscillates along the z-axis. This is the standard model used in virtually all candidate answers and is the correct starting point.

2.  **Analyze the Angular Dependence of Radiated Power:**
    *   For an electric dipole oscillating along the z-axis, the time-averaged power radiated per unit solid angle (`dP/dΩ`) is not uniform. It follows the relationship:
        `dP/dΩ ∝ sin²(θ)`
        where `θ` is the angle measured from the axis of oscillation (the z-axis).
    *   The maximum power, `A`, is radiated where `sin²(θ)` is maximum. This occurs at `θ = 90°` (the plane perpendicular to the oscillation axis), because `sin²(90°) = 1`.
    *   The question asks for the power at `θ = 30°` as a fraction of the maximum power `A`. We need to calculate the value of `sin²(30°)`.
    *   The sine of 30 degrees is a standard value: `sin(30°) = 1/2`.
    *   Therefore, `sin²(30°) = (1/2)² = 1/4`.
    *   The power radiated at 30° is `1/4` of the maximum power `A`. This result is consistently derived by nearly all the agents.

3.  **Analyze the Wavelength (λ) Dependence of Radiated Power:**
    *   The function `f(λ, θ)` represents the power per unit solid angle. Its dependence on the wavelength comes from the fundamental physics of dipole radiation.
    *   The total power radiated by an oscillating electric dipole is proportional to the fourth power of the angular frequency of oscillation (`ω`). This is a result from the Larmor formula.
        `P_total ∝ ω⁴`
    *   The power per unit solid angle has the same frequency dependence: `f ∝ ω⁴`.
    *   The angular frequency `ω` is inversely proportional to the wavelength `λ` (`ω = 2πc/λ`, where `c` is the speed of light).
    *   Substituting this relationship into the power dependence gives:
        `f ∝ (1/λ)⁴ = λ⁻⁴`
    *   Therefore, a possible form for the function `f` must include a `λ⁻⁴` term. This is also consistently derived by the vast majority of the agents.

### Careful Points and Synthesis

*   **Model Validity:** The electric dipole model is the correct and standard approximation for this type of problem. Higher-order multipoles (like quadrupole) exist but are typically much weaker and have different dependencies (e.g., quadrupole power `∝ ω⁶ ∝ λ⁻⁶`). The question is clearly targeting the dominant, lowest-order radiation.
*   **Angular Calculation:** A common mistake is to forget to square the sine. The power is proportional to `sin²(θ)`, not `sin(θ)`. `sin(30°) = 1/2`, but `sin²(30°) = 1/4`.
*   **Frequency vs. Wavelength:** Another key point is the conversion from frequency dependence (`ω⁴`) to wavelength dependence (`λ⁻⁴`). The inverse relationship is crucial.
*   **Candidate Answer Inconsistency:** The provided candidate answers show a strong consensus on the physics: the fraction is **1/4** and the wavelength dependence is **λ⁻⁴**. However, they disagree on which option letter (A, B, C, or D) this corresponds to. This suggests they were working with different or scrambled versions of the options. The correct approach is to ignore their final letter choice and use their correct reasoning to select from the original options provided in the question.

### Conclusion

Based on the physical principles of electric dipole radiation:
1.  The fraction of maximum power at `θ = 30°` is **1/4**.
2.  The wavelength dependence of the radiated power is **λ⁻⁴**.

Now, we match this result to the options given in the problem description:
A) 3/4, λ⁻⁶
B) 1/4, λ⁻³
C) 1/4, λ⁻⁴
D) 1/2, λ⁻⁴

The only option that matches both derived quantities is **C**.

<<<C>>>
"""

# Run the check
result = check_answer(question, llm_answer)
print(result)