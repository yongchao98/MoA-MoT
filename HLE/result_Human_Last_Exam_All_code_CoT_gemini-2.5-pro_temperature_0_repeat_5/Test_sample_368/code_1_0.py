import math

def solve_and_calculate_probability():
    """
    This function calculates the probability that the given series converges
    by iterating through all possible values of X, Y, and Z.
    """
    # Define the ranges for the digits X, Y, Z
    x_range = list(range(-9, 0)) + list(range(1, 10))
    y_range = range(10)
    z_range = range(10)

    # The condition for convergence is |r| < 1, where r = 20u^2 + 24u.
    # This leads to the inequality: -1 < 20u^2 + 24u < 1.
    # Solving this gives two disjoint intervals for u.
    # The boundaries are derived from the roots of 20u^2 + 24u - 1 = 0 and 20u^2 + 24u + 1 = 0.
    sqrt41 = math.sqrt(41)
    sqrt31 = math.sqrt(31)

    # Interval 1: (-6-sqrt(41))/10 < u < (-6-sqrt(31))/10
    u_lower1 = (-6 - sqrt41) / 10
    u_upper1 = (-6 - sqrt31) / 10
    
    # Interval 2: (-6+sqrt(31))/10 < u < (-6+sqrt(41))/10
    u_lower2 = (-6 + sqrt31) / 10
    u_upper2 = (-6 + sqrt41) / 10

    favorable_outcomes = 0
    total_outcomes = 0

    # Iterate through all possible combinations of X, Y, Z
    for x in x_range:
        for y in y_range:
            for z in z_range:
                total_outcomes += 1
                
                # Calculate u = X + Y/10 + 11Z/100
                u = x + y / 10.0 + 11 * z / 100.0

                # Check if u falls into one of the convergence intervals
                is_in_interval1 = u_lower1 < u < u_upper1
                is_in_interval2 = u_lower2 < u < u_upper2

                if is_in_interval1 or is_in_interval2:
                    favorable_outcomes += 1

    # To simplify the fraction, find the greatest common divisor (GCD)
    def gcd(a, b):
        while b:
            a, b = b, a % b
        return a

    common_divisor = gcd(favorable_outcomes, total_outcomes)
    numerator = favorable_outcomes // common_divisor
    denominator = total_outcomes // common_divisor

    # Print the final equation for the probability
    print(f"Number of favorable outcomes: {favorable_outcomes}")
    print(f"Total number of outcomes: {total_outcomes}")
    print(f"The probability is P = {favorable_outcomes} / {total_outcomes} = {numerator} / {denominator}")

solve_and_calculate_probability()