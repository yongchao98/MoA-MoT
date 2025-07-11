import math

def calculate_conductance_moment(n):
    """
    Calculates the n-th statistical moment of the dimensionless conductance <g^n>
    for a disordered wire at the critical point of the class D symmetry.
    
    The conductance g follows the arcsine distribution, and its moments are given by the formula:
    <g^n> = Gamma(n + 0.5) / (sqrt(pi) * n!)
    
    Args:
        n (int): The order of the moment (must be a non-negative integer).
        
    Returns:
        float: The value of the n-th moment <g^n>.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("The moment order 'n' must be a non-negative integer.")
    
    # Calculate the numerator: Gamma(n + 0.5)
    numerator = math.gamma(n + 0.5)
    
    # Calculate the denominator: sqrt(pi) * n!
    denominator = math.sqrt(math.pi) * math.factorial(n)
    
    return numerator / denominator

# The problem asks for the ratio between the fourth moment and the average value.
# Average value is the 1st moment, <g>.
n1 = 1
# Fourth moment is <g^4>.
n4 = 4

# Calculate the average value <g>
average_g = calculate_conductance_moment(n1)

# Calculate the fourth moment <g^4>
fourth_moment_g = calculate_conductance_moment(n4)

# Calculate the final ratio
ratio = fourth_moment_g / average_g

# The exact fractional values are <g> = 1/2 and <g^4> = 35/128.
# The exact ratio is (35/128) / (1/2) = 35/64.

print(f"The average value of the conductance is <g> = {average_g}")
print(f"The fourth statistical moment of the conductance is <g^4> = {fourth_moment_g}")
print("The required ratio is <g^4> / <g>.")
print("The final calculation is:")
print(f"{fourth_moment_g} / {average_g} = {ratio}")
