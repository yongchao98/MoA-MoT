import math

def calculate_b():
    """
    Calculates the exact value of the constant b in the asymptotic formula
    for the expected cover and return time of a random walk on a random tree.

    The formula is C(n) ~ b * n^(3/2).
    The constant b is known from probability theory to be sqrt(2 * pi).
    """
    # The equation for the constant b is b = sqrt(2 * pi)
    # Here are the numbers in that equation:
    number_2 = 2
    constant_pi = math.pi

    # Calculate b
    b = math.sqrt(number_2 * constant_pi)

    # Output the result, showing the numbers used in the equation
    print(f"The equation for the constant is b = sqrt({number_2} * pi)")
    print(f"The exact value of b is the square root of (2 * {constant_pi}).")
    print(f"The numerical value of b is approximately: {b}")

if __name__ == "__main__":
    calculate_b()