import math

def lower_bound_f(n):
    """Calculates the lower bound L(n) for f(n)."""
    if n < 1:
        return 0
    k = math.floor(math.log2(n))
    return n * (k + 1) - (2**(k + 1)) + 1

def upper_bound_f(n):
    """Calculates the upper bound U(n) for f(n)."""
    if n < 2:
        return float('inf')  # Bound is for n>=2
    k_prime = math.floor(math.log2(n - 1))
    return (n - 1) * k_prime - (2**k_prime) + n + 2

def denominator(n):
    """Calculates the denominator n*log2(n)."""
    if n <= 1:
        return 1
    return n * math.log2(n)

def verify_limit():
    """
    Numerically verifies that the limits of L(n)/D(n) and U(n)/D(n) approach 1.
    """
    print("This program demonstrates the convergence of f(n)/(n*log2(n)) to 1.")
    print("We use known lower and upper bounds for f(n) and apply the Squeeze Theorem.")
    print("\nL(n) <= f(n) <= U(n)")
    print("L(n)/(n*log2(n)) <= f(n)/(n*log2(n)) <= U(n)/(n*log2(n))\n")
    print(f"{'n':>10} | {'L(n)/(n*log2(n))':>20} | {'U(n)/(n*log2(n))':>20}")
    print("-" * 55)

    n_values = [10, 100, 1000, 10000, 100000, 1000000, 10000000]

    for n in n_values:
        den = denominator(n)
        lower_ratio = lower_bound_f(n) / den
        upper_ratio = upper_bound_f(n) / den
        print(f"{n:10d} | {lower_ratio:20.10f} | {upper_ratio:20.10f}")

    print("\nAs n gets larger, both the lower and upper bound ratios approach 1.")
    print("By the Squeeze Theorem, the limit must be 1.")

if __name__ == '__main__':
    verify_limit()
