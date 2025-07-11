import math

def approximate_integral():
    """
    This function develops an analytical formula that approximates the integral
    I(epsilon) for a small epsilon.

    The integral is given by:
    I(epsilon) = integral from 0 to 15.00 of 1 / (epsilon + 9.0*x^5.0 + 5.0*x^6.0 + 9.0*x^8.0) dx

    The approximation follows the plan outlined above.
    """

    # Parameters from the dominant term 9.0 * x^5.0
    c = 9.0  # Coefficient of the lowest power term
    n = 5.0  # The lowest power of x

    # Calculate the exponent B for the epsilon term
    B = (1.0 / n) - 1.0

    # Calculate the constant coefficient A
    # The value of the related definite integral is (pi/n) / sin(pi/n)
    integral_part = (math.pi / n) / math.sin(math.pi / n)
    A = math.pow(c, -1.0 / n) * integral_part

    # Print the resulting analytical formula with numerical values
    print("The analytical approximation for I(epsilon) for small epsilon is of the form:")
    print("I(epsilon) = A * epsilon^B")
    print("\nCalculated values:")
    print(f"A = {A:.8f}")
    print(f"B = {B:.1f}")

    # Display the final formula
    print("\nThe final formula is:")
    # The user asked to output each number in the final equation.
    print(f"I(epsilon) = {A:.8f} * epsilon^({B:.1f})")
    
    return A

if __name__ == '__main__':
    approximate_integral()