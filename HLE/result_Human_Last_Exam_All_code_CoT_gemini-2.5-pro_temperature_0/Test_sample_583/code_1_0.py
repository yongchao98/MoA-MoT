import math

def f(n):
    """
    Calculates f(n) based on the formula f(n) = 1 + sum_{k=1 to n} floor(log2(k)).
    """
    if n <= 0:
        return 0
    # The sum can be calculated more efficiently than a direct loop for very large n,
    # but a loop is clear and sufficient for demonstration purposes.
    log_sum = 0
    for k in range(1, n + 1):
        # log2(k) is undefined for k=0, but our range starts from 1.
        # For k=1, log2(1)=0.
        log_sum += math.floor(math.log2(k))
    return 1 + log_sum

def calculate_limit_approximation(n):
    """
    Calculates the ratio f(n) / (n * log2(n)) for a given n.
    """
    if n <= 1:
        return float('nan') # Denominator is 0 or negative

    fn_val = f(n)
    denominator = n * math.log2(n)
    ratio = fn_val / denominator
    
    print(f"For n = {n}:")
    print(f"f(n) = {fn_val}")
    print(f"n * log2(n) = {denominator}")
    print(f"The final equation is: {fn_val} / {denominator} = {ratio}")
    print(f"This ratio numerically approaches the limit.")

# We choose a large value for n to see the trend.
n_large = 1000000
calculate_limit_approximation(n_large)
