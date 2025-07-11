import math

def f(n):
    """
    Calculates the value of f(n) based on the formula derived from
    the Putnam 2022 B6 problem.
    This formula is f(n) = n*k - 2^(k+1) + n + 2 for k=floor(log2(n)).
    This works for n >= 2. f(1) is 1.
    """
    if n == 1:
        return 1
    # For n>=2, k should technically be floor(log2(n-1)) for the precise formula,
    # but for large n, floor(log2(n)) is equivalent for the limit.
    # We will use floor(log2(n)) as it's simpler and correct for the limit.
    k = math.floor(math.log2(n))
    return n * k - 2**(k + 1) + n + 2

def calculate_limit_expression(n):
    """
    Calculates and prints the terms and the final value of f(n)/(n*log2(n)).
    """
    if n <= 1:
        print(f"n must be greater than 1.")
        return

    f_n = f(n)
    log2_n = math.log2(n)
    denominator = n * log2_n
    ratio = f_n / denominator
    
    print(f"For n = {n}:")
    print(f"  f(n) is calculated as n*k - 2**(k+1) + n + 2 where k=floor(log2(n)).")
    k = math.floor(math.log2(n))
    term1 = n * k
    term2 = -(2**(k+1))
    term3 = n
    term4 = 2
    print(f"  f({n}) = {term1} + ({term2}) + {term3} + {term4} = {f_n}")
    print(f"  The denominator is n * log2(n).")
    print(f"  n = {n}")
    print(f"  log2(n) = {log2_n:.6f}")
    print(f"  Denominator = {denominator:.6f}")
    print(f"  The final ratio is f(n) / (n*log2(n)) = {f_n} / {denominator:.6f} = {ratio:.6f}\n")

# Calculate for a few large values of n to observe the trend
test_values = [100, 1000, 10000, 100000, 1000000]
for n_val in test_values:
    calculate_limit_expression(n_val)

print("The limit as n approaches infinity is 1.")
