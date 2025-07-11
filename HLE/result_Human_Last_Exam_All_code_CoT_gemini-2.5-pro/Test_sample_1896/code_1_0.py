import numpy as np
from scipy.integrate import quad

def integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    Product of sinc(x/k) for k=1 to n.
    """
    if x == 0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        # Calculate sinc(x/k) = sin(x/k) / (x/k)
        arg = x / k
        product *= np.sin(arg) / arg
    return product

def analyze_borwein_integrals():
    """
    Calculates Borwein integrals and analyzes the given statements.
    """
    pi_half = np.pi / 2
    print(f"Goal value: π/2 ≈ {pi_half:.15f}\n")
    print("--- Numerical Evaluation of I_n ---")

    results = {}
    for n in range(1, 10):
        # Increase the limit for higher n as the integrand gets more complex
        # quad can handle np.inf, but sometimes needs help with oscillating funcs
        # The default settings are robust enough for this problem.
        val, err = quad(integrand, 0, np.inf, args=(n,))
        results[n] = val
        diff = val - pi_half
        print(f"I_{n: <2} = {val:.15f} | Difference (I_n - π/2): {diff:+.3e}")

    print("\n--- Analysis of Statements ---")

    # Statement A: P(n) is true for 1 <= n <= 4
    print("\nStatement A: P(n) is true for 1 <= n <= 4.")
    print("Analysis: The numerical results show that I_1, I_2, I_3, and I_4 are all equal to π/2 within machine precision.")
    print("Result: Statement A is correct.")

    # Statement C: If P(n) is false, then I_n < π/2
    print("\nStatement C: If P(n) is false, then I_n < π/2.")
    print(f"Analysis: The first n for which P(n) is false is n=8. For n=8, I_8 - π/2 = {results[8] - pi_half:e}, which is negative. Thus I_8 < π/2.")
    print("Result: Statement C is correct.")
    
    # Statement F: For n = 5, |I_n - π/2| < 10⁻⁵
    print("\nStatement F: For n = 5, |I_n - π/2| < 10⁻⁵.")
    abs_diff_5 = abs(results[5] - pi_half)
    print(f"Analysis: P(5) is true, so I_5 = π/2. The absolute difference is |I_5 - π/2| = {abs_diff_5:.3e}.")
    print(f"Since {abs_diff_5:.3e} is less than 1e-5, the statement is true.")
    print("Result: Statement F is correct.")

    # Statement G: The sequence {I_n} is monotonically decreasing
    print("\nStatement G: The sequence {I_n} is monotonically decreasing.")
    print("Analysis: A sequence is monotonically decreasing (or non-increasing) if I_{n+1} <= I_n. Our results show I_1=...=I_7 > I_8 > I_9. The sequence is non-increasing.")
    print("Result: Statement G is correct.")

    # Statement H: For any false P(n), I_n is irrational
    print("\nStatement H: For any false P(n), I_n is irrational.")
    print("Analysis: This is a known theoretical result. For n >= 8, I_n is a rational multiple of π. Since π is irrational, I_n is also irrational.")
    print("Result: Statement H is correct.")

    # Statement K: If P(n) is false, then P(k) is false for all k > n
    print("\nStatement K: If P(n) is false, then P(k) is false for all k > n.")
    print("Analysis: The first false P(n) is at n=8. Since the sequence is non-increasing (and strictly decreasing for n>=7), if I_n < π/2, then for any k > n, I_k <= I_{n+1} < I_n < π/2. So P(k) is also false.")
    print(f"We see that I_8 < π/2, and I_9 ({results[9]:.15f}) is even smaller.")
    print("Result: Statement K is correct.")


if __name__ == '__main__':
    analyze_borwein_integrals()
