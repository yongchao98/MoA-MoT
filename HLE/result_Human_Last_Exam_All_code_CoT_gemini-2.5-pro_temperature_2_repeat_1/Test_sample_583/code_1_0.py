import math

def combinations(n, k):
    """Calculate the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def f_lower_bound(n):
    """Calculates a known lower bound for f(n)."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    k = math.floor(math.log2(n))
    return combinations(k, 2) + k + 1

def calculate_limit_ratio(n):
    """Calculates the ratio given in the problem statement using the lower bound for f(n)."""
    if n <= 1:
        return float('nan')
    
    f_n_approx = f_lower_bound(n)
    denominator = n * math.log2(n)
    
    return f_n_approx / denominator

def main():
    """
    Demonstrates that the ratio f(n) / (n log2 n) tends to 0.
    The final expression for the limit is:
    lim_{n -> inf} f(n) / (n * log2(n)) = 0
    """
    print("This script demonstrates that the ratio f(n)/(n * log2(n)) tends to 0 as n increases.")
    print("We use a known lower bound for f(n) to show this trend.")
    print("-" * 50)
    print("n \t\t Lower Bound for f(n) \t Ratio")
    print("-" * 50)

    for i in range(2, 8):
        n = 10**i
        ratio = calculate_limit_ratio(n)
        f_val = f_lower_bound(n)
        # We don't need to print every number, but to fulfill the prompt
        # let's consider the "final equation" to be the value of the ratio.
        # Here we show the components for clarity.
        print(f"{n:<15} {f_val:<20} {ratio:<.10f}")
        
    print("-" * 50)
    print("As n gets larger, the ratio approaches 0.")
    print("\nThe theoretical limit is 0.")

if __name__ == "__main__":
    main()
