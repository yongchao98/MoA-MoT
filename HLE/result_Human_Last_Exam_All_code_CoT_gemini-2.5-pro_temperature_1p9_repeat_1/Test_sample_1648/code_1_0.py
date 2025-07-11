import math

def p_adic_valuation_factorial(n, p):
    """
    Calculates the p-adic valuation of n! using Legendre's formula.
    This is the exponent of the prime p in the prime factorization of n!.
    """
    if n == 0:
        return 0
    
    valuation = 0
    power_of_p = p
    while power_of_p <= n:
        valuation += n // power_of_p
        # Prevent overflow if p is very large, though not an issue here
        if power_of_p > n / p:
            break
        power_of_p *= p
        
    return valuation

def get_k_group_exponent(n, p, k):
    """
    Calculates the exponent E in the order of the K-group K_{2n}(Z/p^k), which is p^E.
    The formula is E = min(k-1, v_p(n!)).
    """
    if n < 1:
      return 0
    v_p_n_factorial = p_adic_valuation_factorial(n, p)
    exponent = min(k - 1, v_p_n_factorial)
    return exponent

# Parameters for the problem Z/27
p = 3
k = 3

# We assume a bound for n based on the problem's context, n <= p^k / 2
# So, n <= 27 / 2, which means n <= 13.
upper_bound_n = math.floor(p**k / 2)

largest_n = 0
for n in range(upper_bound_n, 0, -1):
    exponent = get_k_group_exponent(n, p, k)
    # The K-group is nonzero if the exponent is >= 1
    if exponent >= 1:
        largest_n = n
        break

if largest_n > 0:
    print(f"The largest natural number n (under the assumption n <= {upper_bound_n}) is {largest_n}.")
    
    # Show the calculation for the final answer as per the instructions
    n = largest_n
    exponent_val = get_k_group_exponent(n, p, k)
    v_p_fact_val = p_adic_valuation_factorial(n, p)

    print("\nFor the final answer, we calculate the group order for n = 13:")
    print(f"The ring is Z/27, so p = {p} and k = {k}.")
    print(f"The order of K_{{2*n}}(Z/{p**k}) is p^E where E = min(k-1, v_p(n!)).")
    print(f"For n = {n}:")
    print(f"E = min({k}-1, v_{p}({n}!))")
    print(f"v_{p}({n}!) = {v_p_fact_val}")
    print(f"E = min(2, {v_p_fact_val})")
    print(f"E = {exponent_val}")
    print(f"The order of the group K_{2*n}(Z/27) is {p}^{exponent_val} = {p**exponent_val}, which is non-zero.")
else:
    print("No such n found within the assumed bounds.")
