import math
import re

def check_exoplanet_temperature_ratio():
    """
    This function checks the correctness of the provided answer for the exoplanet temperature ratio problem.
    """
    
    # --- Define constants and given values from the question ---
    # Mass of Planet 1 in Earth masses
    M_p1 = 7.0
    # Mass of Planet 2 in Earth masses
    M_p2 = 5.0
    # Doppler shift for Planet 1 in Angstroms
    delta_lambda1 = 0.03
    # Doppler shift for Planet 2 in Angstroms
    delta_lambda2 = 0.04

    # --- Perform the calculation based on physics principles ---
    # The derivation shows that T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)
    # All other terms (star's properties, albedo) cancel out in the ratio.
    
    try:
        calculated_ratio = (M_p2 / M_p1) * (delta_lambda1 / delta_lambda2)
    except ZeroDivisionError:
        return "Error: Division by zero in calculation. Check input values."

    # --- Analyze the provided answer ---
    # The provided answer text is given below.
    # We will extract the final letter choice and compare its value to our calculation.
    
    final_answer_text = """
    Here is a step-by-step analysis of the problem, incorporating careful points and addressing potential pitfalls.

    ### **1. Goal Identification**
    The primary objective is to determine the ratio of the equilibrium temperatures between Planet 1 and Planet 2, specifically `T_eq1 / T_eq2`.

    ### **2. Step-by-Step Derivation**

    **Step 2a: Formulate the Equilibrium Temperature Ratio**
    The equilibrium temperature (`T_eq`) of a planet is given by the formula:
    `T_eq = T_star * (1 - A)^(1/4) * sqrt(R_star / (2 * a))`
    where `T_star` is the star's temperature, `R_star` is its radius, `A` is the planet's albedo, and `a` is the semi-major axis of its orbit.

    To find the ratio `T_eq1 / T_eq2`, we divide the formula for Planet 1 by the formula for Planet 2:
    `T_eq1 / T_eq2 = [T_star * (1 - A₁)^(1/4) * sqrt(R_star / (2 * a₁))] / [T_star * (1 - A₂)^(1/4) * sqrt(R_star / (2 * a₂))]`

    The problem states that both planets orbit the same star (`T_star` and `R_star` are constant) and have the same albedo (`A₁ = A₂`). Therefore, these terms, along with the constant `2`, cancel out, simplifying the ratio significantly:
    `T_eq1 / T_eq2 = sqrt(1/a₁) / sqrt(1/a₂) = sqrt(a₂ / a₁)`

    This shows the temperature ratio is the square root of the inverse ratio of the planets' semi-major axes. The problem now becomes finding the ratio `a₂ / a₁`.

    **Step 2b: Relate Semi-Major Axis to Radial Velocity Data**
    The semi-amplitude of the star's radial velocity (`K`) induced by an orbiting planet is given by:
    `K ≈ (M_p / M_star) * sqrt(G * M_star / a)`

    We can rearrange this to solve for the semi-major axis, `a`:
    `K² ≈ (G * M_p²) / (M_star * a)`
    `a ≈ (G * M_p²) / (M_star * K²)`

    **Step 2c: Calculate the Ratio of Semi-Major Axes**
    Using the expression for `a`, we can find the ratio `a₂ / a₁`. The gravitational constant `G` and the star's mass `M_star` are the same for both planets and will cancel out:
    `a₂ / a₁ = [ (G * M_p2²) / (M_star * K₂²) ] / [ (G * M_p1²) / (M_star * K₁²) ]`
    `a₂ / a₁ = (M_p2² / K₂²) * (K₁² / M_p1²) = (M_p2 / M_p1)² * (K₁ / K₂)²`

    **Step 2d: Combine Formulas and Calculate the Final Ratio**
    We can substitute the expression for `a₂ / a₁` back into our temperature ratio formula from Step 2a:
    `T_eq1 / T_eq2 = sqrt(a₂ / a₁) = sqrt( (M_p2 / M_p1)² * (K₁ / K₂)² )`

    The square root and the squares cancel, leaving a simple final expression:
    `T_eq1 / T_eq2 = (M_p2 / M_p1) * (K₁ / K₂)`

    The radial velocity semi-amplitude `K` is directly proportional to the maximum Doppler shift `Δλ`. Therefore, the ratio `K₁ / K₂` is equal to the ratio `Δλ₁ / Δλ₂`.
    `T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)`

    **Step 2e: Substitute Values and Solve**
    Now, we plug in the values given in the problem:
    *   `M_p1 = 7` Earth masses
    *   `M_p2 = 5` Earth masses
    *   `Δλ₁ = 0.03` Å
    *   `Δλ₂ = 0.04` Å

    `T_eq1 / T_eq2 = (5 / 7) * (0.03 / 0.04)`
    `T_eq1 / T_eq2 = (5 / 7) * (3 / 4)`
    `T_eq1 / T_eq2 = 15 / 28`
    `T_eq1 / T_eq2 ≈ 0.5357`

    The calculated ratio is approximately **0.536**.

    ### **3. Careful Points**

    *   **Extraneous Information:** The problem provides a significant amount of information that is not needed to solve for the *ratio* of the temperatures. The star's mass, radius, effective temperature, the planets' radii, and the specific wavelength of the spectral line are all distractors, as these terms cancel out during the derivation.
    *   **Key Assumptions:** The solution relies on the assumption that the orbital inclination `i` is close to 90° for both planets, so that `sin(i) ≈ 1`. This is a well-justified assumption because the problem states the planets were detected by the transit method, which requires the planets' orbits to be viewed nearly edge-on.
    *   **Index Management:** It is crucial to keep the subscripts (1 for Planet1, 2 for Planet2) consistent throughout the calculation. A mix-up, such as inverting one of the ratios, would lead to an incorrect result. The derivation correctly calculates `T_eq1 / T_eq2`.
    *   **Formula Simplification:** The final calculation can be done in two ways: by first calculating the ratio `a₂/a₁` and then taking the square root, or by simplifying the combined formula first. Both methods yield the same result, `15/28`.

    ### **4. Final Answer**

    The calculated ratio of the equilibrium temperatures is approximately 0.536. Comparing this to the given options:
    A) ~ 0.53
    B) ~ 1.05
    C) ~ 0.98
    D) ~ 1.30

    The calculated value is closest to **0.53**.

    <<<A>>>
    """
    
    # Options from the question
    options = {
        "A": 0.53,
        "B": 1.05,
        "C": 0.98,
        "D": 1.30
    }
    
    # Extract the letter from the final answer
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find the final answer in the format <<<X>>> in the provided text."
    
    answer_letter = match.group(1)
    
    # Get the numerical value corresponding to the chosen letter
    answer_value = options.get(answer_letter)
    if answer_value is None:
        return f"Error: The chosen answer letter '{answer_letter}' is not a valid option."

    # --- Verification ---
    # Check if the calculated value is close to the value of the chosen option.
    # A relative tolerance of 2% is used because the options are approximate ("~").
    if math.isclose(calculated_ratio, answer_value, rel_tol=0.02):
        return "Correct"
    else:
        # Find the correct option
        correct_letter = min(options, key=lambda k: abs(options[k] - calculated_ratio))
        return (f"Incorrect. The calculated ratio is {calculated_ratio:.4f}. "
                f"This value is closest to option {correct_letter} ({options[correct_letter]}). "
                f"The provided answer was {answer_letter} ({answer_value}).")

# Run the check and print the result
result = check_exoplanet_temperature_ratio()
print(result)