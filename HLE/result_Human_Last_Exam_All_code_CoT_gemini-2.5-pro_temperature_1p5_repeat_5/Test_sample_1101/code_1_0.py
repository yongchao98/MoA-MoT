import math

def get_sum_val(n):
    """
    Calculates the value of the sum for a given integer n.
    S_n = \sum_{m=0}^n (-1)^m {n\choose m} \frac{\Gamma(\frac{3}{2}+n)}{\Gamma(\frac{3}{2}+n-m)}
    """
    if n < 0:
        return 0

    total = 0.0
    # The gamma ratio can be computed as a product of m terms: (n+1/2)(n-1/2)...
    for m in range(n + 1):
        # Calculate binomial coefficient
        binom = math.comb(n, m)
        
        # Calculate the ratio of Gamma functions
        gamma_ratio = 1.0
        for i in range(m):
            gamma_ratio *= (n + 0.5 - i)
            
        term = binom * gamma_ratio
        
        if m % 2 == 1:
            total -= term
        else:
            total += term
            
    return total

def solve():
    """
    Calculates S_n for n=0 to 9 and prints the results.
    The analysis of the growth of these values suggests f(n) = n!
    is the bounding function with the lowest complexity.
    """
    print("This script calculates the sum for n=0 to 9.")
    print("Based on the growth of the sum, f(n)=n! is the function with the lowest complexity.")
    print("Let's look at the values of S_n and |S_n/n!|:")
    print("n\tS_n\t\t|S_n/n!|")
    print("-" * 40)
    for n in range(10):
        sn = get_sum_val(n)
        n_factorial = math.factorial(n)
        ratio = abs(sn / n_factorial)
        print(f"{n}\t{sn:<15.4f}\t{ratio:.4f}")
        
solve()
