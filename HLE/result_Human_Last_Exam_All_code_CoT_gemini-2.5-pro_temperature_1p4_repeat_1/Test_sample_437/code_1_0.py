import numpy as np

def evaluate_1s_integral(zeta):
    """
    Calculates the value of the integral <phi_1s| 1/r |phi_1s> for a given zeta.

    The analytical result of the integral is simply zeta. This function demonstrates
    the final step of the analytical calculation: I = (4*zeta**3) * (1/(4*zeta**2)).
    
    Args:
        zeta (float): The orbital exponent of the 1s Slater orbital.
    
    Returns:
        float: The calculated value of the integral.
    """
    # The integral can be broken down into a pre-factor from the normalization
    # and the angular integration, and the value of the radial integral itself.
    # Pre-factor = 4 * zeta^3
    pre_factor = 4 * (zeta**3)
    
    # Value of the radial integral: integral from 0 to infinity of r * exp(-2*zeta*r) dr
    # This evaluates to 1 / (4 * zeta^2)
    radial_integral_value = 1 / (4 * (zeta**2))
    
    # The final result is the product of these two terms.
    final_result = pre_factor * radial_integral_value
    
    print("The final step of the analytical calculation is: pre_factor * radial_integral_value")
    print(f"For zeta = {zeta}:")
    print(f"Pre-factor (4 * zeta**3) = {pre_factor}")
    print(f"Radial integral value (1 / (4 * zeta**2)) = {radial_integral_value}")
    
    # "output each number in the final equation"
    # We print the equation with the evaluated numbers.
    print(f"Final Equation: {pre_factor} * {radial_integral_value} = {final_result}")
    
    return final_result

if __name__ == "__main__":
    # For a hydrogen atom, the orbital exponent zeta is 1.0.
    zeta_H = 1.0
    
    print(f"Evaluating the integral <phi_1s| 1/r |phi_1s> for zeta = {zeta_H}\n")
    
    result = evaluate_1s_integral(zeta_H)
    
    print(f"\nThe final result is: {result}")