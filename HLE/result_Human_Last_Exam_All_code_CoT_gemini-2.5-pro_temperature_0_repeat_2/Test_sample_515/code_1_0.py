import math

def solve_random_walk_probability():
    """
    Calculates the approximate probability for the given random walk problem.
    The probability is given by the formula: P = log(2) / (2 * log(3600)).
    """
    
    # Define the components of the equation
    numerator = math.log(2)
    denominator_log_term = math.log(3600)
    
    # Calculate the final probability
    probability = numerator / (2 * denominator_log_term)
    
    # Print the components of the final equation
    print("The final equation for the probability P is: log(2) / (2 * log(3600))")
    print(f"The value of log(2) is: {numerator}")
    print(f"The value of log(3600) is: {denominator_log_term}")
    
    # Print the final calculation
    print(f"P = {numerator} / (2 * {denominator_log_term})")
    
    # Print the final approximate answer
    print(f"The approximate probability is: {probability}")
    
    # Print the answer rounded to two significant digits
    print(f"The answer with two significant digits is: {probability:.2g}")

solve_random_walk_probability()