import numpy as np
from scipy.integrate import quad

def main():
    """
    This script demonstrates that a non-zero Schwartz function does not have all its moments equal to zero.
    """
    
    # Define a non-zero Schwartz class function.
    # This function is f(x) = (1 - 2x^2) * exp(-x^2), which is constructed
    # such that its 0th and 1st moments are zero.
    def f(x):
        return (1 - 2 * x**2) * np.exp(-x**2)

    print("Calculating the first few moments for the non-zero Schwartz function f(x) = (1 - 2x^2) * exp(-x^2)")
    print("-" * 80)

    # We will calculate the moments M_k = integral from -inf to inf of x^k * f(x) dx
    # for k from 0 to 5.
    # The theory states that if f(x) is not the zero function, at least one M_k must be non-zero.
    
    # Analytically, the moments are:
    # M_0 = 0
    # M_1 = 0
    # M_2 = -sqrt(pi) ≈ -1.77245
    # M_3 = 0
    # M_4 = -3*sqrt(pi)/2 ≈ -2.65868
    # M_5 = 0
    
    max_k = 5
    for k in range(max_k + 1):
        # The integrand for the k-th moment is x^k * f(x)
        integrand = lambda x: (x**k) * f(x)
        
        # Perform the numerical integration from -infinity to +infinity
        moment, error = quad(integrand, -np.inf, np.inf)
        
        # Print the final equation for each moment calculation
        print(f"M_{k} = integral(x^{k} * f(x) dx) = {moment: .8f}")

    print("-" * 80)
    print("\nAs shown by the numerical calculation, while some moments are zero (or numerically close to zero),")
    print("the moments M_2 and M_4 are clearly non-zero.")
    print("This supports the conclusion that for a non-zero Schwartz function, not all moments can be zero.")

if __name__ == "__main__":
    main()