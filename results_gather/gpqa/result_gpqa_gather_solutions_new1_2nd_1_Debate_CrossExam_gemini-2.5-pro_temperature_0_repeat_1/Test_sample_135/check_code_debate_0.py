import sympy
from sympy import pi, sin, integrate, symbols

def check_divergence_integral():
    """
    This function verifies the calculation of the volume integral of the divergence
    of the field f(r) = (1/r^2) * r_hat over a sphere of radius R.
    It uses the Divergence Theorem to convert the volume integral into a surface integral.
    """
    
    # 1. Define symbolic variables for the integration
    theta, phi = symbols('theta phi')

    # 2. Define the integrand for the surface integral.
    # According to the Divergence Theorem and the problem statement:
    # f . dS = (1/R^2) * (R^2 * sin(theta)) d(theta) d(phi)
    # The R^2 terms cancel, so the integrand is sin(theta).
    integrand = sin(theta)

    # 3. Perform the double integral over the surface of the sphere.
    # The polar angle theta ranges from 0 to pi.
    # The azimuthal angle phi ranges from 0 to 2*pi.
    try:
        calculated_integral = integrate(integrand, (theta, 0, pi), (phi, 0, 2*pi))
    except Exception as e:
        return f"An error occurred during symbolic integration: {e}"

    # 4. The expected analytical result is 4*pi.
    expected_result = 4 * pi

    # 5. Check if the symbolic calculation matches the expected result.
    if not sympy.simplify(calculated_integral - expected_result) == 0:
        return f"The calculated value of the integral is {calculated_integral}, which does not match the expected value of {expected_result}."

    # 6. Analyze the final answer provided by the LLM.
    # The LLM's final answer correctly calculates the value as 4π.
    # It then maps this value to the options it lists:
    # A) 4/3 π R
    # B) 1
    # C) 4 π
    # D) 0
    # The LLM correctly identifies that 4π corresponds to option C and outputs <<<C>>>.
    
    llm_final_choice = "C"
    options = {
        "A": "4/3 * pi * R",
        "B": "1",
        "C": "4 * pi",
        "D": "0"
    }

    # 7. Verify the LLM's choice.
    # Find which option key corresponds to the correct value.
    correct_option_key = None
    for key, value_str in options.items():
        # We are looking for the string representation of 4*pi
        if value_str == "4 * pi":
            correct_option_key = key
            break
            
    if llm_final_choice == correct_option_key:
        return "Correct"
    else:
        return (f"The final answer's reasoning is correct, but the final choice is wrong. "
                f"The calculated value is {expected_result}, which corresponds to option {correct_option_key}. "
                f"However, the final answer chose option {llm_final_choice}.")

# Run the check
result = check_divergence_integral()
print(result)