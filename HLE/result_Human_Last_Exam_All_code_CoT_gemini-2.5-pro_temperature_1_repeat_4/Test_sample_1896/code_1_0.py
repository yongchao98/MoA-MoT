import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    ∏_{k=1}^n sin(x/k)/(x/k)
    """
    if x == 0.0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        # The term is sin(t)/t where t = x/k
        t = x / k
        if t == 0.0:
            term = 1.0
        else:
            term = np.sin(t) / t
        product *= term
    return product

def main():
    """
    Calculates the first few Borwein integrals and evaluates statement F.
    """
    pi_half = np.pi / 2
    print(f"The value of π/2 is approximately: {pi_half:.10f}\n")

    # Calculate and print I_n for n=1 to 5 for context
    for n in range(1, 6):
        # We pass n to the integrand function using the 'args' parameter of quad
        result, error = quad(borwein_integrand, 0, np.inf, args=(n,))
        print(f"I_{n} = {result:.10f}")
        # P(n) is the proposition that I_n = π/2
        is_p_n_true = np.isclose(result, pi_half)
        print(f"Is P({n}) true? {is_p_n_true}")

        if n == 5:
            # Specifically check statement F for n=5
            # F) For n = 5, |I_n - π/2| < 10⁻⁵
            abs_diff = abs(result - pi_half)
            is_f_true = abs_diff < 1e-5
            print(f"\n--- Checking Statement F ---")
            print(f"For n = 5, the absolute difference |I_5 - π/2| is: {abs_diff:.2e}")
            print(f"Is this difference less than 10⁻⁵? {is_f_true}")
            print(f"Therefore, statement (F) is {is_f_true}.")

if __name__ == "__main__":
    main()
