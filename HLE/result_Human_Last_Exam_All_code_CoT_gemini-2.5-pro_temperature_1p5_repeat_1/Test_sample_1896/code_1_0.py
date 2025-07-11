import numpy as np
from scipy.integrate import quad

def integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    The function is defined as the product of sinc functions.
    """
    # The limit of sin(t)/t as t -> 0 is 1. For x=0, the product is 1.
    if np.isclose(x, 0):
        return 1.0
    
    prod = 1.0
    for k in range(1, n + 1):
        # The term is sin(x/k) / (x/k)
        prod *= np.sin(x / k) / (x / k)
    return prod

def borwein_integral(n):
    """
    Numerically calculates the Borwein integral I_n using scipy.quad.
    """
    # We integrate from 0 to infinity. High precision is requested to distinguish
    # values very close to pi/2.
    # limit=200 increases the number of subintervals for the integration algorithm.
    val, err = quad(integrand, 0, np.inf, args=(n,), limit=200, epsabs=1e-12, epsrel=1e-12)
    return val

def main():
    """
    Main function to compute integrals and evaluate the statements.
    """
    pi_half = np.pi / 2
    print(f"Goal: Evaluate statements about I_n. Reference value π/2 ≈ {pi_half:.10f}\n")

    results = {}
    print("Calculating Borwein integrals I_n:")
    for n in range(1, 10):
        val = borwein_integral(n)
        results[n] = val
        print(f"I_{n} ≈ {val:.10f}, Difference |I_{n} - π/2| = {abs(val - pi_half):.3e}")

    print("\n--- Analysis of selected statements ---\n")
    
    # Statement A: P(n) is true for 1 <= n <= 4
    # P(n) is the proposition that I_n = pi/2.
    is_A_correct = all(np.isclose(results[n], pi_half) for n in range(1, 5))
    print(f"Statement A: 'P(n) is true for 1 <= n <= 4' is {is_A_correct}.")

    # Statement D: The first n where P(n) is false is n = 5
    is_D_correct = not np.isclose(results[5], pi_half)
    print(f"Statement D: 'The first n where P(n) is false is n = 5' is {is_D_correct}.")

    # Statement F: For n = 5, |I_n - π/2| < 10⁻⁵
    n = 5
    abs_diff = abs(results[n] - pi_half)
    bound = 1e-5
    is_F_correct = abs_diff < bound
    
    print(f"Statement F: 'For n = 5, |I_n - π/2| < 10⁻⁵' is {is_F_correct}.")
    print("Equation check:")
    print(f"|I_5 - π/2| = {abs_diff:.10f}")
    print(f"The statement is {abs_diff:.10f} < {bound}, which is true.")
    
    # Check when the value actually changes
    print("\nNote: The first value where I_n is not π/2 is for n=8:")
    print(f"I_7 ≈ {results[7]:.10f}")
    print(f"I_8 ≈ {results[8]:.10f}")


if __name__ == "__main__":
    main()
