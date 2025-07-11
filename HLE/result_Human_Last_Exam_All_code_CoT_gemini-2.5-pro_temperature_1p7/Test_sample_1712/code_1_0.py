import math

def solve_series():
    """
    This function calculates the value of the series sum_{n=2 to inf} (-1)^n * S_n / n,
    where S_n is the n-th harmonic number.
    The value is given by the formula: 1 + (ln(2))^2 / 2 - pi^2 / 12.
    """

    # Get the values for the constants in the formula
    ln2 = math.log(2)
    pi = math.pi
    
    # The numbers used in the final equation
    num_1 = 1
    num_ln2 = ln2
    num_2 = 2
    num_pi = pi
    num_12 = 12

    # Calculate the final result
    result = num_1 + (num_ln2 ** 2) / num_2 - (num_pi ** 2) / num_12

    # Print the equation with the numerical values
    print("The final equation is:")
    print(f"S = {num_1} + ({num_ln2})^2 / {num_2} - ({num_pi})^2 / {num_12}")
    
    # Print the calculated result
    print("\nThe value of the sum is:")
    print(result)

solve_series()