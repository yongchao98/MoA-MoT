import math
import re

def check_correctness_of_answer():
    """
    Checks the correctness of the LLM's answer to the physics question.

    The function verifies the two key physical principles:
    1. The angular dependence of electric dipole radiation power (proportional to sin^2(theta)).
    2. The wavelength dependence of electric dipole radiation power (proportional to lambda^-4).

    It then compares the result with the option chosen by the LLM.
    """

    # 1. Define the problem parameters and calculate the correct physical result.
    angle_deg = 30

    # Principle 1: Angular Dependence
    # The power radiated per unit solid angle for an electric dipole is dP/dΩ ∝ sin²(θ).
    # The maximum power A is at θ = 90°, where sin²(90°) = 1.
    # The fraction of maximum power at a given angle is sin²(θ) / sin²(90°) = sin²(θ).
    angle_rad = math.radians(angle_deg)
    expected_fraction = math.sin(angle_rad)**2

    # Principle 2: Wavelength Dependence
    # The power radiated by an oscillating electric dipole is P ∝ ω⁴.
    # Since ω ∝ 1/λ, the power is P ∝ (1/λ)⁴ = λ⁻⁴.
    expected_lambda_dependence = "lambda**-4"

    # 2. Define the options as provided in the question.
    options = {
        "A": (3/4, "lambda**-6"),
        "B": (1/2, "lambda**-4"),
        "C": (1/4, "lambda**-4"),
        "D": (1/4, "lambda**-3")
    }

    # 3. Parse the LLM's final answer from the provided text.
    llm_answer_text = """
    Here is a step-by-step analysis of the problem, considering the points raised by the various candidate answers.

    ### Step-by-Step Analysis

    1.  **Choose the Correct Physical Model:**
        The problem describes an "oscillating charge distribution" with a spheroid shape and a symmetry axis along the z-axis. This system radiates electromagnetic waves. The standard method to analyze such radiation is the multipole expansion.
        *   **Electric Dipole Model:** This is the lowest-order and most fundamental radiation mechanism. It arises from a time-varying separation of charge, creating an oscillating electric dipole moment. For a charge distribution oscillating along its axis of symmetry, this is the most natural and expected dominant form of radiation. The vast majority of the candidate answers correctly identify this as the appropriate model.
        *   **Electric Quadrupole Model:** This is a higher-order radiation mechanism. It would only be the dominant term if the electric dipole moment were zero due to a specific symmetry in the oscillation (e.g., a shape deformation from prolate to oblate where the center of charge remains fixed). While physically possible, this is a more specific and less common scenario.

        **Conclusion on Model:** In the absence of information forbidding a dipole moment, the electric dipole model is the standard and most appropriate assumption, as it represents the lowest-order, dominant radiation.

    2.  **Analyze the Angular Dependence of Radiated Power:**
        For an electric dipole oscillating along the z-axis, the time-averaged power radiated per unit solid angle (`dP/dΩ`) is given by the well-known formula:
        `dP/dΩ ∝ sin²(θ)`
        where `θ` is the angle measured from the axis of oscillation (the z-axis).

    3.  **Calculate the Fraction of Maximum Power:**
        *   The maximum radiated power per unit solid angle, `A`, occurs when the `sin²(θ)` term is at its maximum.
        *   The maximum value of `sin²(θ)` is 1, which occurs at `θ = 90°` (in the plane perpendicular to the oscillation axis). So, `A` is proportional to 1.
        *   We need to find the power radiated at `θ = 30°`. This is proportional to `sin²(30°)`.
        *   Using the standard trigonometric value, `sin(30°) = 1/2`.
        *   Therefore, the power at 30° is proportional to `sin²(30°) = (1/2)² = 1/4`.
        *   The fraction of the maximum power `A` radiated at `θ = 30°` is the ratio of the power at 30° to the maximum power: `(proportional to 1/4) / (proportional to 1) = 1/4`.

    4.  **Determine the Wavelength (λ) Dependence:**
        *   The power radiated by an oscillating electric dipole is strongly dependent on its frequency. The formula shows that the radiated power is proportional to the fourth power of the angular frequency (`ω⁴`).
        *   The angular frequency `ω` is inversely related to the wavelength `λ` by the formula `ω = 2πc/λ`, where `c` is the speed of light.
        *   Substituting this relationship into the power dependence, we get:
            `Power ∝ ω⁴ ∝ (1/λ)⁴ = λ⁻⁴`
        *   Therefore, a possible form for the function `f(λ, θ)` includes a wavelength dependence of `λ⁻⁴`. This inverse fourth-power law is a famous characteristic of dipole radiation, also responsible for Rayleigh scattering.

    5.  **Combine Results and Select the Final Answer:**
        *   The fraction of maximum power at `θ = 30°` is **1/4**.
        *   The wavelength dependence is **λ⁻⁴**.
        *   We now match this pair of results, `(1/4, λ⁻⁴)`, with the given options:
            A) 3/4, λ⁻⁶
            B) 1/2, λ⁻⁴
            C) 1/4, λ⁻⁴
            D) 1/4, λ⁻³
        *   The derived result perfectly matches option C.

    ### Careful Points to Note

    *   **Model Choice (Dipole vs. Quadrupole):** The most critical step is choosing the correct physical model. The problem is likely designed to test the understanding of which radiation mechanism is dominant. The electric dipole is the lowest-order term and is the standard assumption unless explicitly forbidden. The result for the electric quadrupole model (`3/4`, `λ⁻⁶`) matches option A, making it a plausible but incorrect distractor for this general problem.
    *   **Squaring the Sine Function:** A very common mistake is to use `sin(30°) = 1/2` as the fraction of power. The power is proportional to `sin²(θ)`, so the value must be squared to get `(1/2)² = 1/4`.
    *   **Angle of Maximum Radiation:** The maximum power is radiated perpendicular to the axis of oscillation (`θ = 90°`), not along the axis (`θ = 0°`), where the power is zero.
    *   **Frequency vs. Wavelength Dependence:** The power's dependence on the fourth power of frequency (`ω⁴`) is fundamental. It is crucial to correctly convert this to a dependence on the inverse fourth power of wavelength (`λ⁻⁴`).

    <<<C>>>
    """
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find the final answer in the format <<<A>>>."

    llm_choice = match.group(1)

    # 4. Check if the LLM's choice is physically correct.
    
    # First, determine the correct option key based on our calculation.
    correct_option_key = None
    for key, (fraction, dep) in options.items():
        # Normalize the lambda dependence string for comparison
        normalized_dep = dep.replace("λ", "lambda").replace("^", "**")
        if math.isclose(fraction, expected_fraction) and normalized_dep == expected_lambda_dependence:
            correct_option_key = key
            break

    if correct_option_key is None:
        return "Failure: Could not find a correct option among the choices based on physical principles. The question's options may be flawed."

    # Compare the LLM's choice with the physically correct choice.
    if llm_choice == correct_option_key:
        return "Correct"
    else:
        llm_fraction, llm_dep = options[llm_choice]
        reason = f"Incorrect: The final answer is {llm_choice}, but the correct answer is {correct_option_key}.\n"
        reason += f"The physically correct result is a fraction of 1/4 ({expected_fraction:.2f}) and a wavelength dependence of λ⁻⁴.\n"
        
        is_fraction_wrong = not math.isclose(llm_fraction, expected_fraction)
        is_dep_wrong = llm_dep.replace("λ", "lambda").replace("^", "**") != expected_lambda_dependence

        if is_fraction_wrong and is_dep_wrong:
            reason += f"The chosen option {llm_choice} has an incorrect fraction ({llm_fraction}) and an incorrect wavelength dependence ({llm_dep.replace('**', '^')})."
        elif is_fraction_wrong:
            reason += f"The chosen option {llm_choice} has the correct wavelength dependence but an incorrect fraction ({llm_fraction})."
        elif is_dep_wrong:
            reason += f"The chosen option {llm_choice} has the correct fraction but an incorrect wavelength dependence ({llm_dep.replace('**', '^')})."
            
        return reason

# Execute the check
result = check_correctness_of_answer()
print(result)