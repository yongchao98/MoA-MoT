import math

def evaluate_1s_integral():
    """
    Evaluates the integral <phi_1s| 1/r |phi_1s> for a Slater-type orbital.

    This function demonstrates the step-by-step calculation of the integral,
    which is analytically known to be equal to the orbital exponent, zeta.
    """
    # For a hydrogen atom, the 1s Slater orbital that best approximates
    # the exact wavefunction has an exponent zeta = 1.0.
    zeta = 1.0

    print(f"Evaluating the integral <phi|1/r|phi> for a 1s Slater orbital with exponent zeta = {zeta}")
    
    # The integral is derived as I = (Constant Part) * (Angular Part) * (Radial Part)
    # The full formula is:
    # I = (zeta^3 / pi) * (4 * pi) * (1 / (2*zeta)^2)
    
    # Calculate each component of the equation
    constant_part = (zeta**3) / math.pi
    angular_part = 4 * math.pi
    radial_part = 1 / (2 * zeta)**2
    
    # Calculate the final result by multiplying the parts
    final_result = constant_part * angular_part * radial_part
    
    print("\nThe equation for the integral is:")
    print("I = (zeta^3 / pi) * (4 * pi) * (1 / (4 * zeta^2))")
    
    print("\nEvaluating each number in the final equation:")
    print(f"I = ({constant_part:.4f}) * ({angular_part:.4f}) * ({radial_part:.4f})")
    print(f"I = {final_result:.4f}")
    
    print("\nThe analytical result is zeta.")
    print(f"For zeta = {zeta}, the expected result is {zeta}.")

evaluate_1s_integral()