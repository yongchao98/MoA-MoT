import numpy as np
from scipy.integrate import quad

def borwein_problem_solver():
    """
    Analyzes the Borwein Integral problem, finds the correct statements,
    and prints a detailed evaluation.
    """

    print("--- Borwein Integral Analysis ---")
    
    # 1. State the condition for I_n = pi/2
    print("The value of I_n is pi/2 if and only if Sum_{k=2 to n} (1/k) <= 1.")
    print("Let's find the first 'n' for which this condition fails.\n")

    # 2. Verify the condition for n = 2, 3, 4
    for n_check in range(2, 6):
        terms = [f"1/{k}" for k in range(2, n_check + 1)]
        the_sum = sum(1.0/k for k in range(2, n_check + 1))
        holds = the_sum <= 1
        
        print(f"Checking for n = {n_check}:")
        print(f"Sum = {' + '.join(terms)} = {the_sum:.4f}")
        if holds:
            print(f"Result: {the_sum:.4f} <= 1. The condition holds. So, I_{n_check} should be pi/2.\n")
        else:
            print(f"Result: {the_sum:.4f} > 1. The condition fails. So, I_{n_check} should be less than pi/2.\n")
            first_failure_n = n_check
            break
            
    print(f"The first value where P(n) is false is n = {first_failure_n}.\n")

    # 3. Perform numerical integration to verify
    print("--- Numerical Integration Verification ---")

    # The problem's sinc is sin(x)/x. np.sinc(x) is sin(pi*x)/(pi*x).
    # To compute sin(y)/y, use np.sinc(y/np.pi).
    def my_sinc(x):
        return np.sinc(x / np.pi)

    def integrand(x, n):
        if x == 0.0:
            return 1.0
        prod = 1.0
        for k in range(1, n + 1):
            prod *= my_sinc(x / k)
        return prod

    pi_half = np.pi / 2.0
    results = {}
    print(f"pi/2 = {pi_half:.10f}\n")
    for n in range(1, 7):
        # Integrate from 0 to infinity
        val, err = quad(integrand, 0, np.inf, args=(n,))
        results[n] = val
        print(f"I_{n} = {val:.10f}, Difference from pi/2: {val - pi_half:e}")

    # 4. Evaluate all statements
    print("\n--- Statement Evaluation ---")
    
    # A) P(n) is true for 1 <= n <= 4
    is_A = abs(results[4] - pi_half) < 1e-9
    print(f"A) P(n) is true for 1 <= n <= 4?  {'Correct' if is_A else 'Incorrect'}. (Fails for n=4)")

    # B) P(n) is true for all n
    print(f"B) P(n) is true for all n?  Incorrect. (Fails for n>=4)")

    # C) If P(n) is false, then I_n < pi/2
    # The theory and numerical results (for n=4, 5, 6) confirm this.
    is_C = results[4] < pi_half and results[5] < pi_half
    print(f"C) If P(n) is false, then I_n < pi/2?  {'Correct' if is_C else 'Incorrect'}.")
    
    # D) The first n where P(n) is false is n = 5
    is_D = first_failure_n == 5
    print(f"D) The first n where P(n) is false is n = 5?  {'Correct' if is_D else 'Incorrect'}. (It's n=4)")
    
    # E) lim_{n→∞} I_n = π/4
    print(f"E) lim_{'n→∞'} I_n = π/4?  Incorrect. (This is a known result for a different integral.)")
    
    # F) For n = 5, |I_n - π/2| < 10⁻⁵
    diff_5 = abs(results[5] - pi_half)
    is_F = diff_5 < 1e-5
    print(f"F) For n = 5, |I_n - π/2| < 10⁻⁵?  {'Correct' if is_F else 'Incorrect'}. (Difference is ~{diff_5:.2e})")

    # G) The sequence {I_n} is monotonically decreasing
    is_G = results[2] < results[1] # Check for strict decrease
    print(f"G) The sequence {{I_n}} is monotonically decreasing?  {'Correct' if is_G else 'Incorrect'}. (I_1=I_2=I_3, so not strictly decreasing.)")

    # H) For any false P(n), I_n is irrational
    print(f"H) For any false P(n), I_n is irrational?  Likely true, but a deep result beyond the scope of this problem to verify.")

    # I) Numerical evaluation of I₅ suffices to disprove P(5)
    # P(5) is false. Numerical evaluation shows I_5 is not pi/2.
    is_I = abs(results[5] - pi_half) > 1e-9 
    print(f"I) Numerical evaluation of I₅ suffices to disprove P(5)?  {'Correct' if is_I else 'Incorrect'}.")

    # J) The function under the integral is always positive for n ≤ 4
    print(f"J) The function under the integral is always positive for n <= 4?  Incorrect. (sinc(x) is negative for x > pi).")

    # K) If P(n) is false, then P(k) is false for all k > n
    # Based on the sum condition: if Sum > 1 for n, it will be > 1 for k > n.
    is_K = True 
    print(f"K) If P(n) is false, then P(k) is false for all k > n?  {'Correct' if is_K else 'Incorrect'}.")

    # L) The first four values being π/2 is coincidental
    print(f"L) The first four values being π/2 is coincidental?  Incorrect. (The values being pi/2 up to n=3 is a direct consequence of a mathematical theorem, not coincidence.)")
    
if __name__ == '__main__':
    borwein_problem_solver()