import math

def solve_series_convergence():
    """
    This function calculates the probability that the given geometric series converges.
    It iterates through all possible values of X, Y, and Z, determines if the
    convergence condition |r| < 1 is met, and then computes the probability.
    """
    
    # Define the possible integer values for X, Y, and Z
    # X is in [-9, -1] U [1, 9]
    x_values = list(range(-9, 0)) + list(range(1, 10))
    # Y is in [0, 9]
    y_values = list(range(0, 10))
    # Z is in [0, 9]
    z_values = list(range(0, 10))
    
    # Calculate the total number of possible combinations for (X, Y, Z)
    total_outcomes = len(x_values) * len(y_values) * len(z_values)
    
    # Initialize a counter for the combinations that lead to convergence
    favorable_outcomes = 0
    
    # Iterate through every possible combination of X, Y, and Z
    for x in x_values:
        for y in y_values:
            for z in z_values:
                # Let W = (XYZ/100) + (Z/10).
                # Interpreting XYZ as the value 100*X + 10*Y + Z,
                # W = (100*X + 10*Y + Z)/100 + Z/10
                # W = X + (10*Y + 11*Z) / 100
                w = x + (10 * y + 11 * z) / 100.0
                
                # The common ratio of the geometric series is r = 20*W^2 + 24*W
                r = 20 * w**2 + 24 * w
                
                # A geometric series converges if the absolute value of the common ratio is less than 1
                if abs(r) < 1:
                    favorable_outcomes += 1
                    
    # The final equation for the probability is P = favorable_outcomes / total_outcomes
    # The problem asks to output each number in the final equation.
    # We will print the numerator and denominator of the probability fraction.
    
    # To simplify the fraction, we find the greatest common divisor.
    common_divisor = math.gcd(favorable_outcomes, total_outcomes)
    
    numerator = favorable_outcomes // common_divisor
    denominator = total_outcomes // common_divisor
    
    print(f"{numerator} / {denominator}")

solve_series_convergence()