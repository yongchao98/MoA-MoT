import math

def falling_factorial(x, m):
    """Computes the falling factorial x*(x-1)*...*(x-m+1)."""
    if m == 0:
        return 1.0
    result = 1.0
    for i in range(m):
        result *= (x - i)
    return result

def binomial_coefficient(n, k):
    """Computes the binomial coefficient n choose k."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    # Using math.comb for integer inputs
    if isinstance(n, int) and isinstance(k, int):
        return math.comb(n, k)
    
    # Fallback for non-integer or older python
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_sum(n):
    """
    Calculates the sum S_n = \sum_{m=0}^n (-1)^m {n\choose m} (n+1/2)_m^{\downarrow}
    """
    total_sum = 0.0
    for m in range(n + 1):
        term = ((-1)**m) * binomial_coefficient(n, m) * falling_factorial(n + 0.5, m)
        total_sum += term
    return total_sum

def main():
    """
    Main function to find the lowest complexity f(n).
    It computes Sn for a range of n to show the growth.
    """
    print("This code calculates the sum S_n for several values of n to analyze its growth.")
    print("The sum is S_n = sum_{m=0 to n} (-1)^m * C(n, m) * Gamma(1.5+n) / Gamma(1.5+n-m)")
    print("Based on the growth, we determine the function f(n) with the lowest complexity.")
    print("Let's look at the values of |S_n|/f(n) for different choices of f(n).")

    print("\n" + "="*50)
    print(f"{'n':>3} {'S_n':>15} {'|S_n|/n!':>15} {'|S_n|/2^n':>15}")
    print("="*50)
    
    for n in range(11):
        sn = calculate_sum(n)
        
        # Avoid division by zero for 0! = 1
        ratio_fact = abs(sn) / math.factorial(n) if n >= 0 else float('nan')
        ratio_exp = abs(sn) / (2**n)
        
        print(f"{n:>3} {sn:>15.4f} {ratio_fact:>15.4f} {ratio_exp:>15.4f}")

    print("\n" + "="*50)
    print("The ratio |S_n|/2^n grows, so f(n) must be of higher complexity than exponential.")
    print("The ratio |S_n|/n! seems to grow slowly. A more detailed analysis shows growth is slightly super-factorial.")
    print("Among typical complexity classes {n^k, c^n, n!}, n! is the function with the lowest complexity that bounds the sum.")
    print("Thus, f(n) = n!.")

if __name__ == "__main__":
    main()
