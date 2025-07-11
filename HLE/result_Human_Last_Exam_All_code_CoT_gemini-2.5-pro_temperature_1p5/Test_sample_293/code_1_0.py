import math

def solve_integral():
    """
    Solves the definite integral by simplifying the sum and using a known result.
    """

    # Symbolic representation of the original integral
    integral_string = "∫[0 to ∞] (Σ[n=1 to ∞] log(cos(x/2^n))) dx"

    # The sum simplifies to log(sin(x)/x)
    simplified_integral_string = "∫[0 to ∞] log(sin(x)/x) dx"

    # The known value of this integral is -π/2.
    # We will represent the components of the final answer.
    sign = "-"
    numerator = "π"
    denominator = 2
    value = -math.pi / 2

    print("The integral to solve is:")
    print(integral_string)
    print("\nStep 1: Simplify the summation.")
    print("The summation Σ log(cos(x/2^n)) can be rewritten as log(Π cos(x/2^n)).")
    print("Using the identity Π[n=1 to ∞] cos(x/2^n) = sin(x)/x, the integral becomes:")
    print(simplified_integral_string)
    print("\nStep 2: Evaluate the simplified integral.")
    print("This is a known definite integral in mathematical analysis.")
    print("The value is -π / 2.")
    
    print("\nThe final equation with its components is:")
    # As requested, printing each part of the final equation
    print(f"The integral evaluates to: {sign} {numerator} / {denominator}")

    # For clarity, also printing the numerical value
    print(f"\nNumerical value: {value}")


solve_integral()