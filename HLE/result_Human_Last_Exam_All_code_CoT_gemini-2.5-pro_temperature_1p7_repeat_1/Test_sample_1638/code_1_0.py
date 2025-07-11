import math

def calculate_probability():
    """
    Calculates the probability that for a point p chosen uniformly in the unit square,
    the floor of the reciprocal of the distance from p to at least one vertex is 1.
    """

    # The probability is given by the analytical formula:
    # P = (11/12) * pi - sqrt(3)/2 - 3/2
    
    # Define the numerical components of the formula
    coeff_pi = 11.0 / 12.0
    val_pi = math.pi
    term_sqrt3 = math.sqrt(3) / 2.0
    term_const = 3.0 / 2.0

    # Calculate the final probability
    probability = coeff_pi * val_pi - term_sqrt3 - term_const

    # Print the equation and the breakdown of its components
    print("The probability P is given by the equation:")
    print("P = (11/12) * pi - sqrt(3)/2 - 3/2")
    print("\nSubstituting the numerical values:")
    print(f"P = {coeff_pi} * {val_pi} - {term_sqrt3} - {term_const}")
    print("\nFinal Result:")
    print(f"P = {probability}")

# Execute the calculation
calculate_probability()