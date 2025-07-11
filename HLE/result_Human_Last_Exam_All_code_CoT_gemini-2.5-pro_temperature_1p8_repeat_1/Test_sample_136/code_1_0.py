import math

def calculate_limit_probability():
    """
    This function calculates the limit of the conditional probability
    P[x_0 was not visited | 0 was not visited].
    The limit is given by the expression: exp(-2*pi + pi^2/2).
    """

    # The final expression for the limit
    final_expression = "exp(-2 * pi + (pi^2) / 2)"

    print(f"The limit is given by the expression: {final_expression}")
    print("-" * 20)

    # Define pi
    pi = math.pi

    # Calculate the terms in the exponent
    term1 = -2 * pi
    term2 = pi**2 / 2
    exponent_value = term1 + term2

    # Calculate the final result
    result = math.exp(exponent_value)

    # Print the breakdown of the calculation as requested
    print(f"To evaluate this, we use the value of pi:")
    print(f"pi = {pi}")
    print("\nCalculating the components of the exponent:")
    print(f"First term (-2 * pi) = {term1}")
    print(f"Second term (pi^2 / 2) = {term2}")
    print(f"Value of the exponent = {exponent_value}")
    print("-" * 20)
    print(f"The final numerical result for the limit is: {result}")

if __name__ == "__main__":
    calculate_limit_probability()