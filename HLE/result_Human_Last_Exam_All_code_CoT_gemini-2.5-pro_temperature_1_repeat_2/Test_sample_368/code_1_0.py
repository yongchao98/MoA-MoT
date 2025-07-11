import math

def solve():
    """
    This function calculates the number of favorable outcomes for the series to converge
    and the total number of possible outcomes.
    """
    
    # Define the set of possible values for X, Y, and Z
    x_values = list(range(-9, 0)) + list(range(1, 10))
    y_values = range(10)
    z_values = range(10)
    
    # Calculate the total number of outcomes
    total_outcomes = len(x_values) * len(y_values) * len(z_values)
    
    # Initialize a counter for favorable outcomes
    favorable_outcomes = 0
    
    # Iterate through all possible combinations of (X, Y, Z)
    for x in x_values:
        for y in y_values:
            for z in z_values:
                # The number XYZ is interpreted as 100*X + 10*Y + Z.
                # A = (100*X + 10*Y + Z)/100 + Z/10
                # A = X + (10*Y + 11*Z)/100
                A = x + (10 * y + 11 * z) / 100.0
                
                # The common ratio of the geometric series
                r = 20 * A**2 + 24 * A
                
                # The series converges if |r| < 1
                if abs(r) < 1:
                    favorable_outcomes += 1
                    
    # The final equation for the probability is P = favorable_outcomes / total_outcomes.
    # The following print statements output the numbers in this equation.
    print("Numerator of the probability fraction (favorable outcomes):")
    print(favorable_outcomes)
    print("Denominator of the probability fraction (total outcomes):")
    print(total_outcomes)

solve()