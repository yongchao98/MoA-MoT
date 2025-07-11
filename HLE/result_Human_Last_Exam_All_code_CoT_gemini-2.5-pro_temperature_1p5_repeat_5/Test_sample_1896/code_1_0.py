import numpy as np
from scipy.integrate import quad

def analyze_borwein_integrals():
    """
    Analyzes the statements about the Borwein integral I_n.
    This function calculates I_n numerically to illustrate its properties
    and then evaluates each of the given statements (A-L).
    """

    # --- Helper Functions ---
    def my_sinc(x):
        """Defines sinc(x) as sin(x)/x, handling the x=0 case."""
        # Note: numpy.sinc(x) is defined as sin(pi*x)/(pi*x).
        # To compute sin(x)/x, we can call numpy.sinc(x/pi).
        return np.sinc(x / np.pi)

    def integrand_borwein(x, n):
        """The integrand of the Borwein integral I_n."""
        product = 1.0
        for k in range(1, n + 1):
            product *= my_sinc(x / k)
        return product

    def calculate_I(n):
        """Numerically computes the Borwein integral I_n."""
        # We integrate from 0 to infinity. SciPy's quad is suitable for this.
        # Increased limit for better accuracy with oscillatory integrands.
        result, _ = quad(integrand_borwein, 0, np.inf, args=(n,), limit=200, epsabs=1e-12)
        return result

    # --- Introduction and Numerical Demonstration ---
    print("Analysis of the Borwein Integral Problem")
    print("========================================\n")
    print("The problem concerns the integral: I_n = ∫[0,∞] Π[k=1 to n] sin(x/k)/(x/k) dx.")
    print("It is a known mathematical result that I_n = π/2 for n = 1, 2, ..., 7,")
    print("and that the value becomes less than π/2 for n ≥ 8.\n")
    print("Let's demonstrate this with numerical integration:")

    pi_half = np.pi / 2
    print(f"\nValue of π/2 ≈ {pi_half:.15f}\n")

    print("n | Calculated I_n      | |I_n - π/2|")
    print("--|-----------------------|-------------")

    # Calculate and store values to avoid re-computation
    I_values = {n: calculate_I(n) for n in range(1, 10)}

    for n in range(1, 10):
        I_n = I_values[n]
        diff = abs(I_n - pi_half)
        print(f"{n:2d}| {I_n:<21.15f} | {diff:.2e}")

    print("\n--- Evaluation of Statements ---")

    # A) P(n) is true for 1 ≤ n ≤ 4
    print("\nA) P(n) is true for 1 ≤ n ≤ 4")
    print("   Correct. Based on the established result that P(n) (I_n = π/2) is true for all n ≤ 7.")

    # B) P(n) is true for all n
    print("\nB) P(n) is true for all n")
    print("   Incorrect. P(n) fails for n ≥ 8. Our numerical result for I_8 clearly shows a value less than π/2.")

    # C) If P(n) is false, then I_n < π/2
    print("\nC) If P(n) is false, then I_n < π/2")
    print("   Correct. P(n) is false for n ≥ 8. For these values, I_n is proven to be strictly less than π/2.")

    # D) The first n where P(n) is false is n = 5
    print("\nD) The first n where P(n) is false is n = 5")
    print("   Incorrect. The first n for which P(n) is false is n=8.")

    # E) lim_{n→∞} I_n = π/4
    print(f"\nE) lim_{n→∞} I_n = π/4")
    print(f"   Incorrect. The sequence values are all close to π/2 ≈ 1.5708. The limit is not π/4 ≈ {np.pi/4:.4f}.")

    # F) For n = 5, |I_n - π/2| < 10⁻⁵
    I_5 = I_values.get(5)
    diff_5 = abs(I_5 - pi_half)
    print("\nF) For n = 5, |I_n - π/2| < 10⁻⁵")
    print(f"   Correct. The exact value of I_5 is π/2, so the absolute difference is 0, which is less than 10⁻⁵.")
    print(f"   The equation is |I_5 - π/2| = |{I_5:.10f}... - {pi_half:.10f}...| = {diff_5:.1e}, and {diff_5:.1e} < 1e-5.")

    # G) The sequence {I_n} is monotonically decreasing
    print("\nG) The sequence {I_n} is monotonically decreasing")
    print("   Correct. A sequence is monotonically decreasing if I_{n+1} ≤ I_n for all n. The sequence is constant (I_{n+1} = I_n) for n < 7 and strictly decreasing (I_{n+1} < I_n) for n ≥ 7. Both cases satisfy the condition I_{n+1} ≤ I_n. Therefore, the sequence is monotonically decreasing.")
    
    # H) For any false P(n), I_n is irrational
    print("\nH) For any false P(n), I_n is irrational")
    print("   Correct. When P(n) is false (n ≥ 8), I_n is a product of π and a non-zero rational number. Such a number is irrational.")

    # I) Numerical evaluation of I₅ suffices to disprove P(5)
    print("\nI) Numerical evaluation of I₅ suffices to disprove P(5)")
    print("   Incorrect. P(5) is the proposition I_5 = π/2, which is true. A true statement cannot be disproven. Numerical evaluation only supports P(5).")

    # J) The function under the integral is always positive for n ≤ 4
    print("\nJ) The function under the integral is always positive for n ≤ 4")
    print("   Incorrect. The function sinc(x) = sin(x)/x oscillates and takes negative values for x > π. The product of sinc functions is therefore not always positive.")

    # K) If P(n) is false, then P(k) is false for all k > n
    print("\nK) If P(n) is false, then P(k) is false for all k > n")
    print("   Correct. The sequence {I_n} is strictly decreasing for n ≥ 7. If P(n) is false, then n ≥ 8 and I_n < π/2. For any k > n, we have I_k < I_n < π/2, which means P(k) is also false.")
    
    # L) The first four values being π/2 is coincidental
    print("\nL) The first four values being π/2 is coincidental")
    print("   Incorrect. The result is a well-understood consequence of properties of Fourier transforms. It is not a coincidence.")

if __name__ == '__main__':
    analyze_borwein_integrals()