import math

def get_radon_hurwitz_params(n):
    """
    Decomposes n to find the parameters for the Radon-Hurwitz number calculation.
    n = (2k+1) * 2^m
    m = 4a + b
    Returns m, a, b, and the odd part of n.
    """
    if n == 0:
        return 0, 0, 0, 0
    m = 0
    while n % 2 == 0:
        n //= 2
        m += 1
    odd_part = n
    
    a = m // 4
    b = m % 4
    
    return m, a, b, odd_part

def radon_hurwitz(n, m, a, b):
    """Calculates the Radon-Hurwitz number rho(n)."""
    return 8 * a + 2**b

def solve():
    """
    Finds the number of natural numbers n for which there exist n real n-by-n 
    matrices A_1,...,A_n such that for all nonzero x in R^n, A_1x,...,A_nx are 
    linearly independent.
    This is equivalent to finding the number of natural numbers n such that rho(n) >= n,
    where rho(n) is the Radon-Hurwitz number.
    """
    solutions = []
    # We only need to check up to a certain limit, as it can be proven that
    # for n > 8, rho(n) < n. We check up to n=20 for demonstration.
    limit = 20
    
    print(f"Checking for solutions for n from 1 to {limit}:\n")

    for n in range(1, limit + 1):
        m, a, b, odd_part = get_radon_hurwitz_params(n)
        rho_n = radon_hurwitz(n, m, a, b)
        
        print(f"n = {n}:")
        print(f"  Decomposition: n = {odd_part} * 2^{m}")
        print(f"  Exponent decomposition: m = {m} = 4*{a} + {b}")
        equation_str = f"rho({n}) = 8*{a} + 2^{b} = {rho_n}"
        print(f"  Radon-Hurwitz number: {equation_str}")
        
        condition_met = rho_n >= n
        comparison_symbol = '>=' if condition_met else '<'
        print(f"  Checking condition: rho({n}) >= {n}  ->  {rho_n} {comparison_symbol} {n}")
        
        if condition_met:
            solutions.append(n)
            print("  Result: n is a solution.\n")
        else:
            print("  Result: n is not a solution.\n")
            
    print("-" * 30)
    print(f"The set of natural numbers n that satisfy the condition is: {solutions}")
    print(f"The number of such natural numbers is: {len(solutions)}")

solve()