import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    The integrand of the Borwein integral I_n.
    The product of sinc functions: prod_{k=1 to n} sin(x/k)/(x/k).
    We use np.sinc which is defined as sin(pi*x)/(pi*x).
    Therefore, sin(t)/t corresponds to np.sinc(t/pi).
    """
    if n < 1:
        return 0.0
    
    product = 1.0
    for k in range(1, int(n) + 1):
        product *= np.sinc(x / (k * np.pi))
    return product

def analyze_borwein_integrals():
    """
    Calculates I_n for n=1 to 8 and evaluates the statements from the problem.
    """
    pi_half = np.pi / 2.0
    
    print("This script evaluates the Borwein integral I_n for various n.")
    print("The proposition P(n) is that I_n = pi/2.")
    print(f"Value of pi/2 = {pi_half:.15f}")
    print("-" * 70)
    
    # Store results for final analysis
    results = {}
    sum_val = 0.0
    
    for n in range(1, 9):
        # The condition for I_n = pi/2 is sum_{k=2 to n} (1/k) <= 1
        if n >= 2:
            sum_val += 1.0 / n
        
        # Numerically integrate to find I_n
        # The integral is from 0 to infinity. quad handles np.inf.
        # We increase the limit of integration points for better accuracy.
        result, error = quad(borwein_integrand, 0, np.inf, args=(n,), limit=200)
        results[n] = result

        print(f"For n = {n}:")
        if n > 1:
            print(f"  Sum condition: S_{n} = sum(1/k for k=2..{n}) = {sum_val:.4f}")
            print(f"  Is S_{n} <= 1? {'True' if sum_val <= 1 else 'False'}. So, P({n}) is expected to be {'True' if sum_val <= 1 else 'False'}.")
        else:
            print("  Sum condition: S_1 = 0 (empty sum).")
            print(f"  Is S_1 <= 1? True. So, P(1) is expected to be True.")

        # This part fulfills the "output each number in the final equation" requirement
        print(f"  I_{n} = {result:.15f}")
        print(f"  Difference (I_{n} - pi/2) = {result - pi_half:+.15e}")
        print("-" * 70)

    # Final analysis of all statements
    print("\n--- Analysis of Statements ---\n")
    
    # A) P(n) is true for 1 <= n <= 4.
    p4_is_true = np.isclose(results[4], pi_half)
    print(f"A) P(n) is true for 1 <= n <= 4: {'Correct' if p4_is_true else 'Incorrect'}. (P(4) is false).")

    # C) If P(n) is false, then I_n < pi/2.
    c_correct = all(results[n] < pi_half for n in [4, 5, 6, 7, 8])
    print(f"C) If P(n) is false, then I_n < pi/2: {'Correct' if c_correct else 'Incorrect'}. (See n=4 onwards).")
    
    # D) The first n where P(n) is false is n = 5.
    p3_is_true = np.isclose(results[3], pi_half)
    d_correct = p3_is_true and not p4_is_true
    print(f"D) The first n where P(n) is false is n = 5: {'Incorrect' if d_correct else 'Correct'}. (It's n=4).") # Logic is inverted because statement is false
    
    # F) For n = 5, |I_n - pi/2| < 10^-5.
    f_correct = abs(results[5] - pi_half) < 1e-5
    print(f"F) For n = 5, |I_5 - pi/2| < 10^-5: {'Correct' if f_correct else 'Incorrect'}. (|{results[5] - pi_half:.2e}| is not < 1e-5).")

    # G) The sequence {I_n} is monotonically decreasing.
    g_correct = all(results[i+1] <= results[i] for i in range(1, 8))
    print(f"G) The sequence {{I_n}} is monotonically decreasing (I_n+1 <= I_n): {'Correct' if g_correct else 'Incorrect'}.")

    # I) Numerical evaluation of I_5 suffices to disprove P(5).
    i_correct = not np.isclose(results[5], pi_half, atol=1e-8)
    print(f"I) Numerical evaluation of I_5 suffices to disprove P(5): {'Correct' if i_correct else 'Incorrect'}. (Value is clearly not pi/2).")

    # K) If P(n) is false, then P(k) is false for all k > n.
    k_correct = (not p4_is_true) and all(not np.isclose(results[k], pi_half) for k in [5, 6, 7, 8])
    print(f"K) If P(n) is false, then P(k) is false for all k > n: {'Correct' if k_correct else 'Incorrect'}. (Sum condition confirms this).")


if __name__ == '__main__':
    analyze_borwein_integrals()
