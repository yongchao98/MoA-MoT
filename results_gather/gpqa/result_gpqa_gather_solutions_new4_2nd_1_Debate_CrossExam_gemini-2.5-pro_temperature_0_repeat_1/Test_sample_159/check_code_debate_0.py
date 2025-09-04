import numpy as np
from scipy.special import jn_zeros
import re

def check_diffraction_answer():
    """
    Checks the correctness of the answer for the diffraction problem.

    The problem asks for the angular distance between the first two minima for a
    circular aperture of radius 'a'. The formula for the angular position of the
    m-th minimum is θ_m = (z_m * λ) / (2π * a), where z_m is the m-th non-trivial
    zero of the Bessel function J1(x).

    The angular distance is Δθ = θ₂ - θ₁ = (z₂ - z₁) * λ / (2π * a).
    This code calculates the coefficient C = (z₂ - z₁) / (2π) and compares it
    with the value from the selected option.
    """
    
    # 1. Perform the correct physical calculation
    try:
        # Get the first two non-trivial zeros of the Bessel function J1(x)
        # jn_zeros(n, nt) gives the first nt zeros of Jn(x)
        zeros = jn_zeros(1, 2)
        z1 = zeros[0]  # approx 3.8317
        z2 = zeros[1]  # approx 7.0156
    except ImportError:
        # Fallback if scipy is not installed
        z1 = 3.8317
        z2 = 7.0156

    # Calculate the correct coefficient for the angular distance
    correct_coefficient = (z2 - z1) / (2 * np.pi)

    # 2. Parse the provided LLM answer
    llm_answer_text = """
    Here is a step-by-step analysis of the problem, synthesizing the consensus from the provided answers and highlighting the key points required to arrive at the correct solution.

    ### Step-by-step Analysis

    1.  **Analyze the Aperture's Geometry**: The problem describes an aperture shaped like an N-sided regular polygon with a constant apothem 'a', where N is infinitely large. As the number of sides of a regular polygon with a constant apothem approaches infinity, the polygon's shape converges to a circle. The apothem 'a' becomes the radius of this circle. Therefore, the problem is equivalent to analyzing the diffraction pattern from a circular aperture of radius 'a'.

    2.  **Recall the Physics of Diffraction**: The question asks for the far-field (Fraunhofer) diffraction pattern. For a circular aperture, this pattern is known as an Airy pattern. It consists of a bright central disk surrounded by concentric dark and bright rings.

    3.  **Formulate the Condition for Minima**: The angular positions of the intensity minima (the dark rings) are determined by the zeros of the first-order Bessel function of the first kind, J₁(x). The condition for the m-th minimum is given by:
        `k * a * sin(θ_m) = z_m`
        where:
        *   `k = 2π/λ` is the wave number.
        *   `a` is the radius of the circular aperture.
        *   `θ_m` is the angular position of the m-th minimum.
        *   `z_m` is the m-th non-trivial zero of the Bessel function J₁(x).

    4.  **Apply the Small-Angle Approximation**: The problem states that the related angles are small, so we can use the approximation `sin(θ) ≈ θ`. The formula becomes:
        `(2π/λ) * a * θ_m = z_m`
        Solving for the angular position `θ_m`:
        `θ_m = (z_m * λ) / (2π * a)`

    5.  **Use the Zeros of the Bessel Function**: To find the angular distance between the first two minima, we need the values of the first two non-trivial zeros of J₁(x):
        *   First zero, `z₁ ≈ 3.8317`
        *   Second zero, `z₂ ≈ 7.0156`

    6.  **Calculate the Angular Positions of the Minima**:
        *   The angular position of the first minimum (θ₁) is:
            `θ₁ = (z₁ * λ) / (2π * a) ≈ (3.8317 / 2π) * (λ / a) ≈ 0.6098 * (λ / a)`
        *   The angular position of the second minimum (θ₂) is:
            `θ₂ = (z₂ * λ) / (2π * a) ≈ (7.0156 / 2π) * (λ / a) ≈ 1.1166 * (λ / a)`

    7.  **Calculate the Angular Distance**: The question asks for the angular distance between the first two minima, which is `Δθ = θ₂ - θ₁`.
        `Δθ = θ₂ - θ₁ ≈ (1.1166 - 0.6098) * (λ / a)`
        `Δθ ≈ 0.5068 * (λ / a)`

    8.  **Compare with Options**: The calculated value `Δθ ≈ 0.5068 * (λ / a)` is compared with the given options:
        *   A) 0.506 λ / a
        *   B) 1.220 λ /a
        *   C) 0.500 λ / c (assumed typo for 'a')
        *   D) 0.610 λ / a

        The calculated value is best matched by option A.

    ### Careful Points

    *   **Aperture Geometry**: It is crucial to correctly interpret that the N-sided polygon with a constant apothem 'a' becomes a circle of **radius 'a'**.
    *   **The Specific Question**: The question asks for the *angular distance between* the first two minima (`Δθ = θ₂ - θ₁`), not the position of the first minimum (`θ₁ ≈ 0.610 λ / a`, which is option D) or the angular diameter of the central bright disk (`2θ₁ ≈ 1.220 λ / a`, which is option B). These are common distractors.
    *   **Bessel Function Zeros**: The calculation relies on the specific numerical values of the zeros of the J₁(x) Bessel function. The consistent calculation across all provided analyses confirms the result.

    The analysis shows a clear and unanimous consensus on the physical principles, the calculation steps, and the final numerical result.

    <<<A>>>
    """
    
    # Extract the final answer choice from the text
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find the final answer in the format <<<X>>> in the provided text."
    llm_choice = match.group(1)

    # 3. Define the options from the question
    # Note: Option C has a typo 'c' which is assumed to be 'a'. The coefficient is 0.500.
    options = {
        'A': 0.506,
        'B': 1.220,
        'C': 0.500,
        'D': 0.610
    }
    
    llm_value = options.get(llm_choice)
    if llm_value is None:
        return f"Error: Invalid option '{llm_choice}' was selected."

    # 4. Compare the calculated value with the LLM's chosen value
    tolerance = 1e-3
    if np.isclose(correct_coefficient, llm_value, atol=tolerance):
        return "Correct"
    else:
        # If incorrect, provide a detailed reason and check for common mistakes
        pos_first_min = z1 / (2 * np.pi)
        diameter_central_disk = 2 * pos_first_min
        
        reason = f"Incorrect. The answer chose option {llm_choice} with a value of {llm_value}.\n"
        reason += f"The correct coefficient for the angular distance is (z₂ - z₁) / (2π) ≈ {correct_coefficient:.4f}.\n"
        
        if np.isclose(llm_value, pos_first_min, atol=tolerance):
            reason += f"The chosen value corresponds to the angular position of the first minimum (θ₁ ≈ {pos_first_min:.4f}), which is a common distractor (Option D)."
        elif np.isclose(llm_value, diameter_central_disk, atol=tolerance):
            reason += f"The chosen value corresponds to the angular diameter of the central bright disk (2θ₁ ≈ {diameter_central_disk:.4f}), which is another common distractor (Option B)."
        else:
            reason += f"The correct value {correct_coefficient:.4f} is closest to option A (0.506)."
            
        return reason

# Run the check
print(check_diffraction_answer())