import numpy as np
from scipy.integrate import quad

def borwein_integral_evaluation():
    """
    Calculates and analyzes the Borwein integrals I_n for n=1 to 9.

    This function numerically evaluates the integral:
    I_n = ∫[0 to ∞] Π[k=1 to n] (sin(x/k)/(x/k)) dx

    It prints the value of I_n for each n and its difference from π/2,
    which helps verify the propositions in the Borwein Integral Problem.
    """

    # Define the sinc function sin(x)/x, handling the case x=0.
    # We use np.sinc which is defined as sin(pi*x)/(pi*x), so we adjust the input.
    def sinc(x):
        return np.sinc(x / np.pi)

    # Define the integrand function for I_n.
    def integrand(x, n):
        if n < 1:
            return 0
        
        product = 1.0
        for k in range(1, int(n) + 1):
            product *= sinc(x / k)
        return product

    pi_half = np.pi / 2
    print(f"Analyzing the Borwein Integral I_n for n = 1 to 9")
    print(f"Reference value π/2 ≈ {pi_half:.15f}\n")
    print("-" * 50)

    for n in range(1, 10):
        # The quad function from scipy is used for numerical integration.
        # It can handle integrals to infinity.
        # The 'args' parameter passes n to the integrand function.
        result, error = quad(integrand, 0, np.inf, args=(n,))
        
        print(f"For n = {n}:")
        print(f"I_{n} = ∫[0 to ∞] Π[k=1 to {n}] sinc(x/k) dx")
        print(f"Computed Value: I_{n} ≈ {result:.15f}")
        print(f"Difference |I_{n} - π/2| ≈ {abs(result - pi_half):.2e}")
        print("-" * 50)

if __name__ == '__main__':
    borwein_integral_evaluation()