import math
from scipy.special import gamma, binom

def solve_task(n_max=10):
    """
    Calculates the sum S_n for n from 0 to n_max.
    S_n = sum_{m=0 to n} (-1)^m * (n choose m) * Gamma(1.5+n) / Gamma(1.5+n-m)
    
    The asymptotic behavior of this sum |S_n| is bounded by C * f(n).
    The function f(n) with the lowest complexity is derived from the identity
    S_n = (-1)^n * n! * L_n^(1/2)(1), where L is the generalized Laguerre polynomial.
    The asymptotic analysis shows |S_n| ~ C * n! * exp(2*sqrt(n)) / sqrt(n).
    So, the function f(n) is n! * exp(2*sqrt(n)) / sqrt(n).
    
    This code computes the sum directly to show its values.
    """
    print("This code calculates the sum S_n for several values of n.")
    print("The function f(n) that bounds this sum has been determined to be: n! * exp(2*sqrt(n)) / n^(1/2)")
    print("-" * 30)
    print("n\tS_n")
    print("-" * 30)

    for n in range(n_max + 1):
        # The term (n+1/2)_m can be calculated more stably than the gamma ratio for large n
        # (n+1/2)_m = (n+0.5)*(n-0.5)*...*(n+0.5-m+1)
        total_sum = 0
        for m in range(n + 1):
            falling_factorial = 1.0
            for i in range(m):
                falling_factorial *= (n + 0.5 - i)
            
            term = ((-1)**m) * binom(n, m) * falling_factorial
            total_sum += term
        
        print(f"{n}\t{total_sum:.4f}")

if __name__ == '__main__':
    # Set the maximum value of n to compute
    n_value = 10
    solve_task(n_value)
