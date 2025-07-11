import math

def falling_factorial(x, m):
    """Computes the falling factorial P(x, m) = x * (x-1) * ... * (x-m+1)."""
    if m < 0:
        raise ValueError("m must be a non-negative integer.")
    if m == 0:
        return 1.0
    res = 1.0
    for i in range(m):
        res *= (x - i)
    return res

def binomial_coefficient(n, k):
    """Computes the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    # Using math.comb for precision and efficiency (available in Python 3.8+)
    if hasattr(math, 'comb'):
        return math.comb(n, k)
    # Manual implementation for older Python versions
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_sum(n):
    """Calculates the sum S_n for a given n."""
    x = n + 1.5
    total = 0.0
    for m in range(n + 1):
        term = ((-1)**m) * binomial_coefficient(n, m) * falling_factorial(x - 1, m)
        total += term
    return total

def main():
    """
    Main function to compute the sum for n from 0 to 10 and analyze its growth.
    """
    print("This script analyzes the sum S_n = \\sum_{m=0}^n (-1)^m {n\\choose m} \\frac{\\Gamma(\\frac{3}{2}+n)}{\\Gamma(\\frac{3}{2}+n-m)}")
    print("It calculates S_n for several values of n and compares its growth to f(n)=(1.5)^n * n!")
    print("-" * 80)
    print(f"{'n':>3} {'S_n':>20} {'f(n) = (1.5)^n * n!':>25} {'|S_n| / f(n)':>15}")
    print("-" * 80)
    
    for n in range(11):
        sn = calculate_sum(n)
        fn = (1.5)**n * math.factorial(n) if n > 0 else 1.0
        ratio = abs(sn) / fn if fn != 0 else 0
        
        print(f"{n:>3} {sn:>20.6f} {fn:>25.6f} {ratio:>15.6f}")

if __name__ == "__main__":
    main()
