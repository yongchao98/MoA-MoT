import sympy

def evaluate_1s_integral():
    """
    This function demonstrates the final calculation step for the integral
    <phi_1s| 1/r |phi_1s> for a Slater-type orbital.
    """
    
    # Define zeta as a symbolic variable using sympy for pretty printing
    zeta = sympy.Symbol('zeta', positive=True)

    # The expression for the integral is derived as 4*zeta^3 multiplied by 
    # the integral of r*exp(-2*zeta*r)dr from r=0 to infinity.

    # This is the pre-factor resulting from the normalization constant and angular integration.
    # It is 4 * pi * (zeta**3 / pi) = 4 * zeta**3
    pre_factor = 4 * zeta**3

    # The result of the radial integral part: Integral(r * exp(-2*zeta*r) dr)
    # is 1 / (4*zeta^2).
    radial_integral_result = 1 / (4 * zeta**2)

    # The final result is the product of the pre-factor and the radial integral result.
    final_result = pre_factor * radial_integral_result

    # Print the breakdown of the final calculation, showing each part of the equation.
    print("The final step of the evaluation combines the pre-factor with the result of the radial integral.")
    print("The equation is: [Pre-Factor] * [Radial Integral Result] = Final Value")
    print("\nHere are the components:")
    print(f"Pre-Factor: {pre_factor}")
    print(f"Radial Integral Result: {radial_integral_result}")
    
    print("\nThe final equation with all numbers and symbols is:")
    # The f-string automatically converts the sympy objects to strings for printing.
    print(f"{pre_factor} * ({radial_integral_result}) = {final_result}")

if __name__ == '__main__':
    evaluate_1s_integral()