import sympy

def solve_sum():
    """
    Calculates the sum by multiplying the values of zeta(6) and zeta(8).
    It then prints the derivation of the final answer.
    """
    # Define pi as a symbolic object
    pi = sympy.pi

    # sympy.zeta() provides the exact values for even integer arguments
    zeta_6 = sympy.zeta(6)
    zeta_8 = sympy.zeta(8)

    # The total sum is the product of zeta(6) and zeta(8)
    total_sum = zeta_6 * zeta_8

    # Extract the denominators for printing the equation
    # For zeta(2k), the result is (rational) * pi^(2k)
    # We find the rational part by dividing by pi^n
    z6_rational_part = zeta_6 / pi**6
    z8_rational_part = zeta_8 / pi**8
    
    # The numerator of the rational part is 1, so we just need the denominator
    z6_denominator = z6_rational_part.q
    z8_denominator = z8_rational_part.q

    # Calculate the final denominator
    final_denominator = z6_denominator * z8_denominator

    # Print the final calculation step-by-step
    print("The sum is the product of the Riemann zeta functions ζ(6) and ζ(8).")
    print(f"Sum = ζ(6) * ζ(8)")
    print(f"    = (π^6 / {z6_denominator}) * (π^8 / {z8_denominator})")
    print(f"    = π^(6 + 8) / ({z6_denominator} * {z8_denominator})")
    print(f"    = π^14 / {final_denominator}")
    print(f"\nThe result is (1/{final_denominator}) * π^14.")

solve_sum()