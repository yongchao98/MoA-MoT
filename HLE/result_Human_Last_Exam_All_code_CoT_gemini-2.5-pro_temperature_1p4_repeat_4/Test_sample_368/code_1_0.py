import math

def solve():
    """
    This function calculates the probability that the given series converges
    by iterating through all possible values of X, Y, and Z.
    """
    
    # These are the bounds for the variable 'v' for the series to converge.
    # The condition for convergence |20v^2 + 24v| < 1 solves to:
    # v in (-1.2403, -1.1568) U (-0.0432, 0.0403)
    v_lower_1 = (-6 - math.sqrt(41)) / 10
    v_upper_1 = (-6 - math.sqrt(31)) / 10
    v_lower_2 = (-6 + math.sqrt(31)) / 10
    v_upper_2 = (-6 + math.sqrt(41)) / 10

    # Define the possible integer values for X, Y, and Z
    x_values = list(range(-9, 0)) + list(range(1, 10))
    y_values = range(10)
    z_values = range(10)

    favorable_outcomes = 0
    total_outcomes = 0

    # Iterate through all combinations of X, Y, Z
    for x in x_values:
        for y in y_values:
            for z in z_values:
                total_outcomes += 1
                
                # We interpret "XYZ" as the numerical value 100*X + 10*Y + Z.
                # The variable 'v' from the problem statement is calculated as:
                v = (100 * x + 10 * y + z) / 100.0 + z / 10.0

                # Check if v falls into one of the two convergence intervals
                is_convergent = (v > v_lower_1 and v < v_upper_1) or \
                                (v > v_lower_2 and v < v_upper_2)

                if is_convergent:
                    favorable_outcomes += 1

    # The probability is the ratio of favorable outcomes to total outcomes.
    # We find the greatest common divisor to simplify the fraction.
    common_divisor = math.gcd(favorable_outcomes, total_outcomes)
    
    # Print the numbers that form the final probability equation
    print(f"Number of favorable outcomes for convergence: {favorable_outcomes}")
    print(f"Total number of possible outcomes: {total_outcomes}")
    print(f"The simplified probability is the fraction: {favorable_outcomes // common_divisor}/{total_outcomes // common_divisor}")

solve()