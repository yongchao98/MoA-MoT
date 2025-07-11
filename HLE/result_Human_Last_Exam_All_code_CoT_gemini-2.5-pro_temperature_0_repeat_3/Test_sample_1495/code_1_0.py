import math

def solve_probability():
    """
    This function calculates the exact probability that a random point P inside
    triangle ABC falls within the inner triangle XYZ, which is formed by three
    random cevians.

    The problem is a classic one in geometric probability. The final probability
    is the result of a complex integral derived from Routh's theorem on areas.
    The analytical solution for this integral is known to be 5 - (pi^2 / 2).

    This code calculates and prints the components of this expression and the
    final result.
    """

    # Define the constants in the expression
    five = 5.0
    two = 2.0
    pi = math.pi

    # Calculate the terms of the expression
    pi_squared = pi ** 2
    pi_squared_over_two = pi_squared / two
    
    # Calculate the final probability
    probability = five - pi_squared_over_two

    # Print the explanation and the components of the final equation
    print("The exact probability is given by the analytical formula: 5 - (pi^2 / 2)")
    print("\nHere are the components of the final equation:")
    print(f"The number five: {five}")
    print(f"The number pi: {pi}")
    print(f"The number pi squared (pi^2): {pi_squared}")
    print(f"The number two: {two}")
    print(f"The term (pi^2 / 2): {pi_squared_over_two}")

    # Print the final calculation step-by-step
    print("\nThe final calculation is:")
    print(f"{five} - {pi_squared_over_two} = {probability}")

    # Print the final answer
    print(f"\nThe exact probability is approximately: {probability}")

solve_probability()