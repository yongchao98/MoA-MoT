import math

def f_fast(n):
    """
    Calculates f(n) using the derived closed-form formula.
    f(n) = n(floor(log2(n-1))+1) - 2^(floor(log2(n-1))+1) + 2
    """
    if n <= 0:
        raise ValueError("n must be a positive integer")
    if n == 1:
        return 1
    
    m = n - 1
    L = math.floor(math.log2(m))
    
    val = n * (L + 1) - 2**(L + 1) + 2
    return val

def calculate_limit_ratio(n):
    """
    Calculates the ratio f(n) / (n * log2(n)).
    """
    if n <= 1:
        return float('nan')
        
    fn_val = f_fast(n)
    denominator = n * math.log2(n)
    
    ratio = fn_val / denominator
    return ratio

if __name__ == "__main__":
    print("This script calculates the ratio f(n)/(n*log2(n)) for increasing values of n to observe the limit.")
    
    ns_to_test = [10**2, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
    
    for n_val in ns_to_test:
        ratio_val = calculate_limit_ratio(n_val)
        fn_val = f_fast(n_val)
        n_log2_n = n_val * math.log2(n_val)

        # Per the instruction: "output each number in the final equation"
        # We interpret this as showing the components of the calculation.
        print(f"\nFor n = {n_val}:")
        print(f"f(n) = {fn_val}")
        print(f"n*log2(n) = {n_log2_n:.2f}")
        print(f"The ratio f(n) / (n * log2(n)) is: {ratio_val:.8f}")
