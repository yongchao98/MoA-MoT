import math

def p_adic_valuation(n, p):
    """Computes the p-adic valuation of n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count

def check_conditions(n, p, k):
    """
    Checks if n satisfies the conditions for the K-group to be non-zero.
    """
    # Condition 1: n must be a multiple of p for the group to be non-zero.
    if n % p != 0:
        return False
        
    # Condition 2: The formula validity condition.
    val_p_n = p_adic_valuation(n, p)
    log_val = math.floor(math.log(2 * n) / math.log(p))
    
    if k >= val_p_n + log_val:
        return True
    else:
        return False

def solve_problem():
    """
    Finds the largest natural number n such that the (2n)th K-group
    of Z/27 is nonzero.
    """
    p = 3
    k = 3
    
    # We will search for the largest n satisfying the conditions.
    # The function f(n) = nu_p(n) + floor(log_p(2n)) is increasing with n,
    # so we can find the point where it exceeds k and the largest n
    # will be just below that. Let's test values of n.
    
    largest_n = 0
    
    # Let's test a reasonable range of n. Since log(2n) grows slowly,
    # we don't need to check very far.
    for n in range(1, 100):
        if check_conditions(n, p, k):
            if n > largest_n:
                largest_n = n

    print(f"Let's analyze the conditions for n:")
    print(f"1. n must be a multiple of p = {p}.")
    print(f"2. k >= nu_p(n) + floor(log_p(2*n)), with k = {k}.")
    print("\nLet's test some values of n that are multiples of 3:")
    
    test_cases = [3, 6, 9, 12, 15, 18]
    for n in test_cases:
        nu_3_n = p_adic_valuation(n, p)
        log_term = math.floor(math.log(2*n)/math.log(p))
        sum_val = nu_3_n + log_term
        is_valid = k >= sum_val
        print(f"n = {n}:")
        print(f"  nu_3({n}) = {nu_3_n}")
        print(f"  floor(log_3(2*{n})) = {log_term}")
        print(f"  Sum = {nu_3_n} + {log_term} = {sum_val}")
        print(f"  Condition {k} >= {sum_val} is {is_valid}.")
        if is_valid:
            order_log = p_adic_valuation(n,p)
            print(f"  --> K_{2*n}(Z/27) is non-zero (order is 3^{order_log}).")
        else:
            print(f"  --> Formula for K_{2*n}(Z/27) is not valid, group is likely zero.")
        print("-" * 20)
        
    print(f"\nBased on the analysis, the largest value of n satisfying the conditions is {largest_n}.")
    print(f"For n={largest_n}, the condition is: {k} >= {p_adic_valuation(largest_n, p)} + {math.floor(math.log(2*largest_n)/math.log(p))}, which is {k} >= {p_adic_valuation(largest_n, p) + math.floor(math.log(2*largest_n)/math.log(p))}")
    print(f"For n={largest_n+3}, the condition is: {k} >= {p_adic_valuation(largest_n+3, p)} + {math.floor(math.log(2*(largest_n+3))/math.log(p))}, which is {k} >= {p_adic_valuation(largest_n+3, p) + math.floor(math.log(2*(largest_n+3))/math.log(p))}")

solve_problem()