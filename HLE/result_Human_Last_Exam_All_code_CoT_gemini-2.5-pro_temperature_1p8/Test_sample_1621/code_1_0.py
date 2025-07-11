import math

def get_n_decomposition(n):
    """
    Decomposes n into m * 2^b where m is odd.
    Returns (m, b).
    """
    if n == 0:
        return 0, 0
    b = 0
    while n % 2 == 0 and n != 0:
        n //= 2
        b += 1
    m = n
    return m, b

def radon_hurwitz(n):
    """
    Calculates the Radon-Hurwitz number rho(n).
    """
    if n <= 0:
        return 0
    
    # Step 1: Decompose n = m * 2^b
    m, b = get_n_decomposition(n)
    
    # Step 2: Decompose b = c + 4d
    c = b % 4
    d = b // 4
    
    # Step 3: Calculate rho(n) = 2^c + 8d
    rho = (2**c) + (8 * d)
    return rho, m, b, c, d

def solve():
    """
    Finds all natural numbers n for which n <= rho(n).
    """
    solutions = []
    print("The problem requires finding all natural numbers n such that n <= rho(n).")
    print("Searching for solutions...\n")
    # We only need to check a limited range, as n grows much faster than rho(n).
    for n in range(1, 30):
        rho, m, b, c, d = radon_hurwitz(n)
        
        is_solution = (n <= rho)
        
        if is_solution:
            solutions.append(n)
            print(f"Found a solution for n = {n}:")
            print(f"  Decomposition: {n} = {m} * 2^{b}")
            print(f"  Exponent b={b} gives: c = {c}, d = {d}")
            print(f"  rho({n}) = 2^{c} + 8*{d} = {2**c} + {8*d} = {rho}")
            print(f"  Checking condition: {n} <= {rho} is True.\n")
        
    print("="*40)
    print(f"The possible values for n are: {solutions}")
    print(f"In total, there are {len(solutions)} such natural numbers.")

solve()