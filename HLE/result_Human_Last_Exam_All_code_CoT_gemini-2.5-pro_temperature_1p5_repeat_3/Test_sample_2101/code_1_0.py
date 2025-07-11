import math

def solve_electron_escape_problem():
    """
    Calculates the difference between the probability that an electron escapes
    through the hypotenuse of an isosceles right triangle and the probability
    that it escapes through either of the two legs.

    The solution is based on a known result from geometrical probability.
    - The probability of escaping through the hypotenuse, P_H, is 2/pi.
    - The probability of escaping through the legs, P_L, is 1 - P_H.
    - The desired difference is P_H - P_L = 2 * P_H - 1.
    - Substituting the value of P_H gives 2 * (2/pi) - 1 = 4/pi - 1.
    """

    # Numbers used in the final equation
    numerator = 4
    subtrahend = 1
    
    # The mathematical constant Pi
    pi_value = math.pi

    # Calculate the final result
    result = numerator / pi_value - subtrahend

    # Output the final equation with the numbers used
    print(f"The final equation is: {numerator} / \u03C0 - {subtrahend}")

    # Output the calculated difference
    print(f"The difference between the probabilities is: {result}")

solve_electron_escape_problem()