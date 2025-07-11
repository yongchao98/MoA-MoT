import math

def evaluate_limit_expression(n):
    """
    Calculates the value of the expression f(n) / (n * log2(n)) for a given n.
    
    Args:
        n: A positive integer.
    """
    # Step 2: f(n) = C(n, 2) + 1
    # C(n, 2) = n * (n - 1) / 2
    f_n = n * (n - 1) // 2 + 1
    
    # Denominator
    n_log2_n = n * math.log2(n)
    
    # The expression value
    value = f_n / n_log2_n
    
    print(f"For n = {n}:")
    # The instruction "output each number in the final equation" is interpreted
    # as showing the components of the calculation.
    print(f"f({n}) = C({n}, 2) + 1 = {n}({n}-1)/2 + 1 = {f_n}")
    print(f"n * log_2(n) = {n} * log_2({n}) = {n_log2_n:.4f}")
    print(f"The value of the expression is {f_n} / {n_log2_n:.4f} = {value:.4f}")
    print("-" * 30)

# Demonstrate the growth for increasing values of n
evaluate_limit_expression(100)
evaluate_limit_expression(1000)
evaluate_limit_expression(10000)
evaluate_limit_expression(100000)
