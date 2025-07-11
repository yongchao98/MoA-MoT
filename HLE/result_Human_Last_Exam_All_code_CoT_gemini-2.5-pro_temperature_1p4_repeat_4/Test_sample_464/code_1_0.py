import numpy as np
from scipy.integrate import quad
import math

def illustrate_moments():
    """
    This function illustrates that a non-zero Schwartz function must have at least one non-zero moment.
    """
    # Define a non-zero Schwartz function f(x) = (1 - 2x^2) * exp(-x^2).
    # This function is chosen because its 0-th moment is zero.
    f = lambda x: (1 - 2*x**2) * np.exp(-x**2)

    print("Analyzing the moments of the non-zero Schwartz function f(x) = (1 - 2x^2) * exp(-x^2)")
    print("The theorem states that if f(x) is not the zero function, at least one moment must be non-zero.")
    print("-" * 70)
    
    # Calculate moments M_k = integral of x^k * f(x) from -inf to inf
    max_k = 10
    has_nonzero_moment = False
    
    for k in range(max_k + 1):
        # Define the integrand for the k-th moment
        integrand = lambda x: x**k * f(x)
        
        # Perform numerical integration
        moment, error = quad(integrand, -np.inf, np.inf)
        
        # We know from theory that M_k = 0 for odd k, and for even k=2n, M_{2n} = -2n * Gamma(n+1/2)
        # For k=0, M_0 = 0.
        # For k=2, M_2 = -sqrt(pi).
        # For k=4, M_4 = -3*sqrt(pi).
        
        # Check if the moment is non-zero (within a small tolerance)
        if abs(moment) > 1e-9:
            has_nonzero_moment = True
        
        print(f"Moment M_{k} = integral(x^{k} * f(x) dx) = {moment:.8f}")

    print("-" * 70)
    if has_nonzero_moment:
        print("As expected, we found non-zero moments (e.g., M_2, M_4, ...).")
        print("This supports the contrapositive: since the function is not zero, its moments are not all zero.")
        print("Therefore, if all moments were zero, the function would have to be zero.")
    else:
        # This case should not be reached for this function
        print("Something went wrong, all calculated moments are zero.")


if __name__ == "__main__":
    illustrate_moments()
