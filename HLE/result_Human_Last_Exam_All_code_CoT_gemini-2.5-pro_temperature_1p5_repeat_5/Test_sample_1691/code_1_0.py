import math

def solve_integral_approximation():
    """
    Calculates the coefficients for the analytical approximation of the integral.
    The integral is of the form I(e) = integral[0 to 15] (1 / (e + 9x^5 + 5x^6 + 9x^8)) dx.
    For small epsilon, this is approximated by integral[0 to inf] (1 / (e + a*x^p)) dx,
    where a*x^p is the dominant term for small x.
    In this case, a = 9.0 and p = 5.0.

    The analytical formula for the approximation is A * epsilon^B, where:
    B = (1/p) - 1
    A = (pi / p) / (a^(1/p) * sin(pi/p))
    """
    
    # Parameters from the dominant term (9.0 * x^5.0)
    a = 9.0
    p = 5.0
    
    # Calculate the exponent B
    B = (1/p) - 1
    
    # Calculate the coefficient A
    A = (math.pi / p) / (a**(1/p) * math.sin(math.pi / p))
    
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    # The prompt requests that we output each number in the final equation.
    print(f"I(epsilon) approx = {A} * epsilon^({B})")

solve_integral_approximation()

# The question asks for the answer in a specific format at the end.
# Based on examples, it seems to request a single numerical value.
# Let's provide the calculated coefficient A.
a = 9.0
p = 5.0
A = (math.pi / p) / (a**(1/p) * math.sin(math.pi / p))
# print(f"<<<{A}>>>")