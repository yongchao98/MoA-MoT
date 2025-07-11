import math

def evaluate_slater_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater-Type Orbital (STO).
    """
    # The orbital exponent zeta. For a Hydrogen atom, zeta is 1.0.
    zeta = 1.0

    print("Evaluating the integral <phi_1s| 1/r |phi_1s> for a 1s Slater-Type Orbital.")
    print(f"The chosen orbital exponent is zeta = {zeta:.2f}, corresponding to a Hydrogen atom.\n")
    
    print("The normalized 1s STO is given by: phi_1s(r) = (zeta^3 / pi)^(1/2) * exp(-zeta*r)")
    print("The integral evaluates to the analytical result: zeta\n")

    print("We can show this by calculating the full expression:")
    print("Result = (4 * zeta^3) * Integral(r * exp(-2*zeta*r) dr)")
    print("The integral part evaluates to: 1 / (4 * zeta^2)\n")

    print("Let's compute this step-by-step with our value of zeta:")
    
    # Calculate the components of the final equation: (4 * zeta^3) * (1 / (4 * zeta^2))
    
    # First term: 4 * zeta^3
    four_zeta_cubed = 4 * (zeta**3)
    
    # Denominator of the second term: 4 * zeta^2
    four_zeta_squared = 4 * (zeta**2)
    
    # The full second term: 1 / (4 * zeta^2)
    inv_four_zeta_squared = 1 / four_zeta_squared
    
    # The final result
    result = four_zeta_cubed * inv_four_zeta_squared

    print("Final Equation: (4 * zeta^3) * (1 / (4 * zeta^2)) = Result")
    print(f"Substituting zeta = {zeta:.2f}:")
    print(f"({4} * {zeta:.2f}^3) * (1 / ({4} * {zeta:.2f}^2)) = ?")
    print(f"({4} * {zeta**3:.2f}) * (1 / ({4} * {zeta**2:.2f})) = ?")
    print(f"({four_zeta_cubed:.2f}) * (1 / {four_zeta_squared:.2f}) = ?")
    print(f"{four_zeta_cubed:.2f} * {inv_four_zeta_squared:.2f} = {result:.2f}\n")

    print(f"The final evaluated result of the integral is {result:.2f}")

if __name__ == "__main__":
    evaluate_slater_integral()