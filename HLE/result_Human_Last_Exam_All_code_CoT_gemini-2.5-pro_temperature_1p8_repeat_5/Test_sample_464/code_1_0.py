import numpy as np
from scipy.integrate import quad
from scipy.special import hermite

def illustrate_moments():
    """
    This function illustrates that we can construct non-zero Schwartz functions
    for which any finite number of moments are zero. This highlights how
    strong the condition of *all* moments being zero is.

    We use functions of the form f_n(x) = H_n(x) * exp(-x^2), where H_n(x)
    is the n-th Hermite polynomial. These are Schwartz functions.

    A known property of these functions is that the k-th moment,
    integral(x^k * f_n(x) dx), is zero for all integers 0 <= k < n.
    """
    
    # Let's choose n=4. We expect the moments for k=0, 1, 2, 3 to be zero,
    # and the moment for k=4 to be non-zero.
    n = 4
    
    # Get the polynomial object for the n-th Hermite polynomial
    Hn = hermite(n)
    
    # Define our Schwartz function f_n(x) = H_n(x) * e^(-x^2)
    def f_n(x):
        return Hn(x) * np.exp(-x**2)

    print(f"Testing the moments of the Schwartz function f(x) = H_{n}(x) * exp(-x^2)")
    print(f"where H_{n}(x) is the Hermite polynomial of degree {n}:")
    print(f"H_{n}(x) = {Hn}\n")
    print("The k-th moment is the integral of x^k * f(x) from -inf to inf.")
    print("We expect moments for k < 4 to be zero.\n")

    # We will compute the moments for k from 0 up to n+1
    max_k = n + 1
    
    for k in range(max_k + 1):
        # Define the integrand for the k-th moment
        integrand = lambda x: (x**k) * f_n(x)
        
        # Perform the numerical integration from -infinity to +infinity
        moment, error = quad(integrand, -np.inf, np.inf)
        
        # The result of the integration might be a very small number close to
        # zero due to numerical precision. We can round very small results.
        if abs(moment) < 1e-10:
            moment = 0.0
        
        print(f"The moment for k = {k}:")
        # For a clearer representation of the equation for each moment:
        # Since printing the polynomial can be long, we describe the function.
        equation = f"∫ x^{k} * (H_{n}(x) * e^(-x^2)) dx"
        print(f"{equation} ≈ {moment:.4f}")

illustrate_moments()
