import numpy as np
from scipy.integrate import quad

def my_sinc(x):
    """
    Custom sinc function sin(x)/x for numerical stability at x=0.
    """
    if x == 0:
        return 1.0
    return np.sin(x) / x

def integrand(x, n):
    """
    The function under the integral for I_n.
    """
    product = 1.0
    for k in range(1, n + 1):
        product *= my_sinc(x / k)
    return product

def evaluate_statements():
    """
    Calculates integrals and provides information to check the statements.
    """
    pi_half = np.pi / 2
    print(f"π/2 ≈ {pi_half:.10f}\n")
    
    print("Numerical evaluation of I_n:")
    integrals = {}
    errors = {}
    for n in range(1, 6):
        # We integrate from 0 to infinity
        val, err = quad(integrand, 0, np.inf, args=(n,))
        integrals[n] = val
        errors[n] = err
        # The final question format doesn't want me to print words like 'true/false', but numbers, so let's format it in a neutral way
        print(f"I_{n} = {val:.10f} | Difference from π/2: {val - pi_half:+.3e} | Error Est: {err:.3e}")
        
    print("\n--- Information for Statement Analysis ---")
    
    # Information for A & D: P(n) is true for 1 <= n <= 4? First false n?
    p_values = {n: np.isclose(integrals[n], pi_half) for n in range(1, 5)}
    print(f"Truthiness of P(n) for n=1,2,3,4: {p_values[1]}, {p_values[2]}, {p_values[3]}, {p_values[4]}")
    first_false_n = next((n for n in range(1, 6) if not np.isclose(integrals[n], pi_half)), None)
    print(f"First n for which P(n) is false: {first_false_n}")

    # Information for C: If P(n) is false, then I_n < π/2? Test with n=4
    print(f"Check for n=4: P(4) is {np.isclose(integrals[4], pi_half)}, I_4 < π/2 is {integrals[4] < pi_half}")
    
    # Information for F: For n = 5, |I_n - π/2| < 10⁻⁵?
    diff_5 = abs(integrals[5] - pi_half)
    print(f"For n=5, |I_n - π/2| ≈ {diff_5:.3e}. Condition < 1e-5 is {diff_5 < 1e-5}")

    # Information for G: Is sequence monotonically decreasing?
    # Monotonically decreasing means I_n > I_{n+1} for all n.
    print(f"Check for monotonicity: I_1 ≈ I_2 is {np.isclose(integrals[1], integrals[2])}. So it is not strictly decreasing.")

    # Information for I: Numerical evaluation of I_5 suffices to disprove P(5)?
    # This requires |I_5 - π/2| > integration error
    print(f"Check for disproving P(5): |I_5 - π/2| ({diff_5:.3e}) > Error ({errors[5]:.3e}) is {diff_5 > errors[5]}")

    # Information for K: If P(n) is false, then P(k) is false for all k > n?
    # Test with n=4, k=5. P(4) is false, is P(5) also false?
    print(f"Check for K: P(4) is {np.isclose(integrals[4], pi_half)}, P(5) is {np.isclose(integrals[5], pi_half)}")


if __name__ == '__main__':
    evaluate_statements()
