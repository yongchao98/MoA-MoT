import numpy as np
from scipy.integrate import quad

def borwein_integrand(x, n):
    """
    Calculates the product of sinc functions for the Borwein integral I_n.
    The integrand is prod_{k=1 to n} sin(x/k)/(x/k).
    """
    # The value of the integrand at x=0 is 1, as lim_{t->0} sin(t)/t = 1.
    if x == 0.0:
        return 1.0
    
    product = 1.0
    for k in range(1, n + 1):
        arg = x / k
        # The sinc function sin(t)/t is 1 at t=0.
        if arg == 0.0:
            product *= 1.0
        else:
            product *= np.sin(arg) / arg
    return product

def solve_borwein_problem():
    """
    Calculates I_n for n=1 to 8 and checks the condition for I_n = pi/2.
    """
    pi_half = np.pi / 2
    
    print("Analysis of Borwein Integrals I_n")
    print(f"The reference value is pi/2 ≈ {pi_half:.10f}\n")
    print("The condition for I_n = pi/2 is: Sum_{k=2 to n} (1/k) < 1")
    print("-" * 85)
    print(f"{'n':<3} | {'Sum_{k=2 to n} 1/k':<20} | {'Condition Met?':<16} | {'Value of I_n':<20} | {'I_n - pi/2':<20}")
    print("-" * 85)
    
    # The sum in the condition starts from k=2.
    # We initialize it to 0 for n=1 (empty sum).
    current_sum = 0.0
    
    for n in range(1, 9):
        # Update the sum for the current n (for n>=2)
        if n >= 2:
            current_sum += 1.0 / n
            
        condition_met = current_sum < 1
        
        # Calculate the integral I_n using numerical quadrature.
        # We increase the limit of subdivisions for better accuracy.
        I_n, error = quad(borwein_integrand, 0, np.inf, args=(n,), limit=200)
        
        # Format for printing
        sum_str = f"{current_sum:.4f}" if n >= 2 else "0.0000"
        cond_str = "Yes" if condition_met else "No"
        val_str = f"{I_n:.10f}"
        diff_str = f"{I_n - pi_half:+.3e}"
        
        print(f"{n:<3} | {sum_str:<20} | {cond_str:<16} | {val_str:<20} | {diff_str:<20}")

    print("-" * 85)
    print("\nAnalysis of Statements:")
    print("A) P(n) for 1<=n<=4: False, P(4) is false.")
    print("C) If P(n) false, I_n < pi/2: True, as seen for n>=4.")
    print("D) First false P(n) is n=5: False, it's n=4.")
    print("G) Sequence is monotonically decreasing: True. I_1=I_2=I_3 > I_4 > I_5...")
    print("K) If P(n) false, P(k) for k>n is false: True, the sum is always increasing.")

if __name__ == '__main__':
    solve_borwein_problem()
