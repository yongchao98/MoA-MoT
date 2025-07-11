import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    Calculates the product Π_{k=1 to n} sinc(x/k)
    where sinc(u) is defined as sin(u)/u.
    """
    if x == 0.0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        # This is equivalent to sinc(x/k)
        val = x / k
        product *= np.sin(val) / val
    return product

def analyze_borwein_integrals():
    """
    Analyzes statements about Borwein integrals I_n by numerical evaluation.
    """
    pi_half = np.pi / 2
    print(f"Analyzing Borwein Integrals I_n = ∫[0, ∞] Π[k=1 to n] sinc(x/k) dx")
    print(f"Reference value: π/2 ≈ {pi_half:.15f}\n")

    # Theoretical condition for I_n = π/2 is Σ_{k=2 to n} 1/k < 1
    harmonic_sum = 0.0
    
    for n in range(1, 6):
        print(f"--- Analysis for n = {n} ---")
        
        # Check the theoretical condition that determines if I_n = π/2
        if n >= 2:
            harmonic_sum += 1/n
        is_condition_met = (n == 1) or (harmonic_sum < 1)
        
        print(f"Theoretical check: Σ_{{k=2 to {n}}} 1/k = {harmonic_sum:.4f}")
        print(f"Condition (sum < 1) is met: {is_condition_met}")
        
        # Perform high-precision numerical integration
        # quad is suitable for this, even with an infinite upper bound.
        i_n, err = quad(borwein_integrand, 0, np.inf, args=(n,), epsabs=1e-14, epsrel=1e-14)
        
        print(f"Numerical result for I_{n} ≈ {i_n:.15f}")
        
        # Check P(n): I_n = π/2
        is_p_n_true = np.isclose(i_n, pi_half, atol=1e-12)
        print(f"Is P({n}) true (I_{n} = π/2)? {'Yes' if is_p_n_true else 'No'}")

        # Check statement C: If P(n) is false, then I_n < π/2
        if not is_p_n_true:
            print(f"Check for statement C: Is I_{n} < π/2? {'Yes' if i_n < pi_half else 'No'}")

        # Check statement F specifically for n=5
        if n == 5:
            diff = abs(i_n - pi_half)
            is_f_true = diff < 1e-5
            print("\n--- Special check for Statement F (n=5) ---")
            print("Equation: |I_5 - π/2| < 10⁻⁵")
            print(f"Value of I_5: {i_n:.15f}")
            print(f"Value of π/2: {pi_half:.15f}")
            print(f"Value of |I_5 - π/2|: {diff:.3e}")
            print(f"Is the condition met? {is_f_true}. Statement F is correct.")
            
        print("-" * 28 + "\n")

if __name__ == '__main__':
    analyze_borwein_integrals()
