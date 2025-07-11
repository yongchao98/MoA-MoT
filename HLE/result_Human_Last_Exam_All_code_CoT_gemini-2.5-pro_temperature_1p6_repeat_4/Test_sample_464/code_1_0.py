import numpy as np
from scipy.integrate import quad
import math

def demonstrate_moments():
    """
    This function demonstrates that a non-zero Schwartz function has non-zero moments.
    We use the Gaussian function f(x) = exp(-x^2) as our example.
    """
    
    # Define the Schwartz function.
    def f(x):
        """The Gaussian function, a classic example of a Schwartz function."""
        return np.exp(-x**2)

    print("Let's test a non-zero Schwartz function, f(x) = exp(-x^2).")
    print("According to the theorem, this function must have at least one non-zero moment.")
    print("We will compute its moments M_k = integral from -inf to inf of x^k * f(x) dx.\n")

    # Calculate and print moments for k from 0 to 5.
    for k in range(6):
        # The integrand for the k-th moment is x^k * f(x).
        integrand = lambda x: (x**k) * f(x)

        # Perform numerical integration from -infinity to +infinity.
        # quad returns a tuple: (result, estimated_error)
        result, _ = quad(integrand, -np.inf, np.inf)

        # The analytical value for the k-th moment of exp(-x^2) is known.
        # It is 0 for odd k, and Gamma((k+1)/2) for even k.
        if k % 2 != 0:
            analytical_result = 0.0
        else:
            analytical_result = math.gamma((k + 1) / 2.0)
        
        print(f"Moment M_{k} = integral(x^{k} * exp(-x^2) dx):")
        print(f"  - Numerical Result: {result: .8f}")
        print(f"  - Analytical Value: {analytical_result: .8f}")

        # Highlight whether the moment is zero or non-zero.
        if abs(result) > 1e-9:
            print("  -> This moment is NON-ZERO.")
        else:
            print("  -> This moment is effectively ZERO.")
        print("-" * 35)

    print("\nConclusion: As you can see, the even moments (M_0, M_2, M_4) are non-zero.")
    print("This supports the theorem: since f(x) = exp(-x^2) is not the zero function,")
    print("its moments cannot all be zero.")

if __name__ == '__main__':
    demonstrate_moments()
