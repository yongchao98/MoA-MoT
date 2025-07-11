import math

def solve_expected_value():
    """
    Calculates the expected value E[T] based on the derived formula 3 - 2*e^(1/3).
    It also prints the components of the final equation as requested.
    """
    # The derived formula for the expected value E[T] is 3 - 2 * e^(1/3)
    constant_term = 3
    coefficient = -2
    exponent = 1/3

    # Calculate the value of e^(1/3)
    val_e_pow_1_3 = math.exp(exponent)

    # Calculate the final result
    result = constant_term + coefficient * val_e_pow_1_3

    # Print the final equation and its components
    print("The expected value E[T] is calculated by the formula: 3 - 2 * e^(1/3)")
    print("\nBreaking down the calculation:")
    print(f"The formula is: {constant_term} + ({coefficient}) * e^({exponent:.2f})")
    print(f"First, calculate e^({exponent:.2f}): {val_e_pow_1_3:.8f}")
    print(f"Next, multiply by the coefficient {coefficient}: {coefficient} * {val_e_pow_1_3:.8f} = {coefficient * val_e_pow_1_3:.8f}")
    print(f"Finally, add the constant term {constant_term}: {constant_term} + ({coefficient * val_e_pow_1_3:.8f}) = {result:.8f}")
    print(f"\nThus, the final expected value E[T] is approximately: {result:.8f}")

if __name__ == '__main__':
    solve_expected_value()