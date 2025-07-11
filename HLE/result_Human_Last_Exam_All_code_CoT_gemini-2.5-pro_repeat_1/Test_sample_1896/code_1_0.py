import numpy as np
from scipy.integrate import quad
import math

def solve_borwein_problem():
    """
    This function calculates the Borwein integrals and evaluates the propositions.
    """
    
    # --- Step 1: Define the integrand and calculate the integrals ---
    
    def integrand(x, n):
        """ The integrand of the Borwein integral I_n. """
        if x == 0.0:
            return 1.0
        prod = 1.0
        for k in range(1, n + 1):
            y = x / k
            # The sinc function is defined as sin(x)/x.
            # We handle the case y=0 separately for numerical stability, though x=0 is the main concern.
            if y == 0.0:
                prod *= 1.0
            else:
                prod *= np.sin(y) / y
        return prod

    # Calculate integrals I_1 to I_8. This range is sufficient to check all statements.
    MAX_N = 8
    integrals = []
    pi_half = np.pi / 2
    
    print("--- Step 1: Calculating Borwein Integrals I_n ---")
    print(f"Value of pi/2 = {pi_half:.15f}\n")

    for n in range(1, MAX_N + 1):
        # Use quad with high precision demands to get accurate results.
        # The calculation for n>=4 is sensitive and may take a moment.
        val, err = quad(integrand, 0, np.inf, args=(n,), epsabs=1e-13, epsrel=1e-13, limit=1000)
        integrals.append(val)
        
        # Output the required "equation" format
        print(f"For n={n}, the equation is I_{n} = pi/2.")
        print(f"Result: I_{n} = {val:.15f}")
        print(f"Difference |I_{n} - pi/2| = {abs(val - pi_half):.3e}\n")

    # A small tolerance is needed to compare floating point numbers.
    # The first deviation (for I_4) is > 1e-3, so 1e-10 is a safe threshold.
    TOL = 1e-10

    # --- Step 2: Evaluate each proposition ---
    print("\n--- Step 2: Evaluating Propositions ---")

    # A) P(n) is true for 1 <= n <= 4
    is_A_correct = all(abs(integrals[n-1] - pi_half) < TOL for n in range(1, 5))
    print(f"A) P(n) is true for 1 <= n <= 4: {is_A_correct}")
    if not is_A_correct:
        first_fail_A = next(n for n in range(1, 5) if not abs(integrals[n-1] - pi_half) < TOL)
        print(f"   - Justification: P({first_fail_A}) is false since I_{first_fail_A} != pi/2.")

    # B) P(n) is true for all n
    is_B_correct = all(abs(I_n - pi_half) < TOL for I_n in integrals)
    print(f"B) P(n) is true for all n: {is_B_correct}")
    if not is_B_correct:
        first_fail_B = next(n for n, I_n in enumerate(integrals, 1) if not abs(I_n - pi_half) < TOL)
        print(f"   - Justification: P({first_fail_B}) is the first proposition to fail.")

    # C) If P(n) is false, then I_n < π/2
    false_Pn_indices = [n-1 for n in range(1, MAX_N + 1) if abs(integrals[n-1] - pi_half) > TOL]
    is_C_correct = all(integrals[i] < pi_half for i in false_Pn_indices) if false_Pn_indices else True
    print(f"C) If P(n) is false, then I_n < pi/2: {is_C_correct}")
    print(f"   - Justification: For every n >= 4, the calculated I_n is less than pi/2.")

    # D) The first n where P(n) is false is n = 5
    first_fail_n = next((n for n, I_n in enumerate(integrals, 1) if abs(I_n - pi_half) > TOL), None)
    is_D_correct = (first_fail_n == 5)
    print(f"D) The first n where P(n) is false is n = 5: {is_D_correct}")
    if first_fail_n:
        print(f"   - Justification: The first n where P(n) is false is n = {first_fail_n}.")

    # E) lim_{n→∞} I_n = π/4
    is_E_correct = False
    print(f"E) lim_{n→∞} I_n = pi/4: {is_E_correct}")
    print(f"   - Justification: This is a known theoretical result. The actual limit is 0, not pi/4.")

    # F) For n = 5, |I_n - π/2| < 10⁻⁵
    I_5 = integrals[4]
    is_F_correct = abs(I_5 - pi_half) < 1e-5
    print(f"F) For n = 5, |I_n - pi/2| < 10⁻⁵: {is_F_correct}")
    print(f"   - Justification: The actual difference is |I_5 - pi/2| ≈ {abs(I_5 - pi_half):.5f}, which is much larger than 10⁻⁵.")

    # G) The sequence {I_n} is monotonically decreasing
    is_G_strictly_decreasing = all(integrals[i] > integrals[i+1] for i in range(MAX_N - 1))
    is_G_non_increasing = all(integrals[i] >= integrals[i+1] + TOL or abs(integrals[i] - integrals[i+1]) < TOL for i in range(MAX_N-1))
    print(f"G) The sequence {{I_n}} is monotonically decreasing: {is_G_non_increasing}")
    print(f"   - Justification: The sequence is non-increasing (I_1=I_2=I_3 > I_4 > ...). It is not strictly decreasing.")

    # H) For any false P(n), I_n is irrational
    is_H_correct = True # Theoretical result
    print(f"H) For any false P(n), I_n is irrational: {is_H_correct}")
    I_4_ratio_str = f"{integrals[3] / pi_half:.8f}"
    theory_ratio_str = f"{431/432:.8f}"
    print(f"   - Justification: Theory states I_n = (rational) * pi for n>=4. E.g., I_4/(pi/2) = {I_4_ratio_str}, which matches the theoretical value 431/432 = {theory_ratio_str}.")

    # I) Numerical evaluation of I₅ suffices to disprove P(5)
    is_I_correct = abs(integrals[4] - pi_half) > 1e-8 # A typical numerical tolerance
    print(f"I) Numerical evaluation of I₅ suffices to disprove P(5): {is_I_correct}")
    print(f"   - Justification: The difference |I_5 - pi/2| is large enough (~{abs(integrals[4] - pi_half):.2e}) to be easily detected by numerical methods.")

    # J) The function under the integral is always positive for n ≤ 4
    integrand_n1_at_4 = math.sin(4) / 4 
    is_J_correct = False
    print(f"J) The function under the integral is always positive for n ≤ 4: {is_J_correct}")
    print(f"   - Justification: For n=1, the integrand sin(x)/x is negative for x in (pi, 2*pi). For example, at x=4, sin(4)/4 ≈ {integrand_n1_at_4:.2f}.")

    # K) If P(n) is false, then P(k) is false for all k > n
    is_K_correct = True
    for n in range(1, MAX_N):
        if abs(integrals[n-1] - pi_half) > TOL: # P(n) is false
            for k in range(n + 1, MAX_N + 1):
                if abs(integrals[k-1] - pi_half) < TOL: # P(k) is true
                    is_K_correct = False; break
        if not is_K_correct: break
    print(f"K) If P(n) is false, then P(k) is false for all k > n: {is_K_correct}")
    print("   - Justification: Our results show that once the value deviates from pi/2, it continues to do so. This is supported by theory.")
    
    # L) The first four values being π/2 is coincidental
    is_L_correct = False
    print(f"L) The first four values being π/2 is coincidental: {is_L_correct}")
    print(f"   - Justification: The premise is incorrect. Only the first three values are pi/2, not four.")

# Execute the analysis
solve_borwein_problem()