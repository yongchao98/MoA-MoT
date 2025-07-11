import numpy as np
from scipy.integrate import quad

def demonstrate_moments():
    """
    Numerically computes the moments of the Gaussian function f(x) = exp(-x^2).
    A Schwartz function f(x) must be the zero function if all its moments are zero.
    Since f(x) = exp(-x^2) is not the zero function, we expect to find non-zero moments.
    """
    
    # Define a non-zero Schwartz function (the Gaussian function)
    f = lambda x: np.exp(-x**2)
    
    print("Calculating moments M_k = integral(x^k * f(x) dx) from -inf to inf for f(x) = exp(-x^2)\n")

    max_k = 10
    for k in range(max_k + 1):
        # Define the integrand for the k-th moment
        integrand = lambda x: (x**k) * f(x)
        
        # Perform the numerical integration
        # quad returns a tuple (result, estimated_error)
        moment, error = quad(integrand, -np.inf, np.inf)
        
        # The theoretical results are 0 for odd k, and non-zero for even k.
        # Due to numerical precision, results for odd k will be very close to zero.
        # We can round very small numbers to 0 for clarity.
        if np.abs(moment) < 1e-9:
            moment = 0.0
            
        print(f"M_{k} = integral(x^{k} * exp(-x^2)) dx = {moment:.4f}")

demonstrate_moments()