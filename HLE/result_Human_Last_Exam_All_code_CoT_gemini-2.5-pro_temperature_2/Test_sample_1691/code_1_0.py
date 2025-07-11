import math

def get_approximation_formula():
    """
    Calculates the constants for the analytical approximation of the integral I(epsilon).
    The integral is I(epsilon) = integral from 0 to 15 of 1/(epsilon + 9x^5 + 5x^6 + 9x^8) dx.
    For small epsilon, this is approximated by C * epsilon^a.
    """
    
    # Parameters from the leading term f(x) approx a * x^p near x=0
    a = 9.0
    p = 5.0
    
    # The exponent of epsilon in the approximation is -(p-1)/p
    exponent_val = -(p - 1.0) / p
    
    # The coefficient C is given by a^(-1/p) * integral(0,inf) of 1/(1+u^p) du
    # The integral part is (pi/p)/sin(pi/p)
    integral_const = (math.pi / p) / math.sin(math.pi / p)
    coefficient_c = math.pow(a, -1.0 / p) * integral_const
    
    # Print the derived analytical formula with the computed constants
    print("The analytical formula that approximates I(epsilon) for small epsilon is:")
    print(f"I(epsilon) \u2248 {coefficient_c} * \u03B5^({exponent_val})")
    
if __name__ == '__main__':
    get_approximation_formula()
