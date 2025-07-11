import math

def solve():
    """
    Calculates the probability that the series converges by iterating through all
    possible values of X, Y, and Z.
    """
    favorable_cases = 0
    
    # Define the ranges for X, Y, and Z
    x_values = list(range(-9, 0)) + list(range(1, 10))
    y_values = range(10)
    z_values = range(10)
    
    total_cases = len(x_values) * len(y_values) * len(z_values)
    
    # Pre-calculate the boundaries of the convergence intervals for u
    sqrt41 = math.sqrt(41)
    sqrt31 = math.sqrt(31)
    
    # Interval 1: (-1.2403, -1.1568)
    lower_bound_1 = (-6 - sqrt41) / 10
    upper_bound_1 = (-6 - sqrt31) / 10
    
    # Interval 2: (-0.0432, 0.0403)
    lower_bound_2 = (-6 + sqrt31) / 10
    upper_bound_2 = (-6 + sqrt41) / 10
    
    # Iterate through all possible combinations of (X, Y, Z)
    for x in x_values:
        for y in y_values:
            for z in z_values:
                # Calculate u = X + Y/10 + 11Z/100
                u = x + y / 10.0 + (11 * z) / 100.0
                
                # Check if u falls into one of the convergence intervals
                is_in_interval_1 = (lower_bound_1 < u < upper_bound_1)
                is_in_interval_2 = (lower_bound_2 < u < upper_bound_2)
                
                if is_in_interval_1 or is_in_interval_2:
                    favorable_cases += 1
                    
    print(f"The number of combinations (X, Y, Z) for which the series converges is: {favorable_cases}")
    print(f"The total number of possible combinations for (X, Y, Z) is: {total_cases}")
    print(f"The probability of convergence is the ratio of these two numbers.")
    print(f"P = {favorable_cases} / {total_cases}")

solve()