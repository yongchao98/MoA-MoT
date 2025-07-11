import math

# Memoization cache for Fibonacci numbers
fib_cache = {}

def fib(k):
    """Computes the k-th Fibonacci number with memoization."""
    if k in fib_cache:
        return fib_cache[k]
    if k < 0:
        raise ValueError("Index cannot be negative")
    if k == 0:
        return 0
    if k == 1:
        return 1
    result = fib(k - 1) + fib(k - 2)
    fib_cache[k] = result
    return result

def solve():
    """
    Calculates the number of solutions for a fixed large m, assuming
    higher-order coefficients (a,b,c,d,e) are zero.
    """
    m = 30  # Choose a sufficiently large m
    f_m = fib(m)
    
    # We need to find the number of integer tuples (n, f, g) such that:
    # F_n + f * F_m + g = 0
    # with f, g in [-25, 25] and n >= 0.

    solution_count = 0
    
    # The range of n to check needs to be determined.
    # If F_n + f*F_m is outside [-25, 25], there's no solution for g.
    # |f*F_m| <= F_n + 25. For f!=0, this means F_n is large.
    # F_n is approx phi^n/sqrt(5). F_n/F_m is approx phi^(n-m).
    # phi^(n-m) must be around |f| <= 25.
    # n-m < log(25)/log(phi) approx 6.7. So we check n up to m+10.
    n_limit = m + 10
    
    print(f"Checking for m = {m} (F_{m} = {f_m})")
    print("Solutions are tuples (n, f, g) satisfying F_n + f*F_m + g = 0")
    print("-" * 50)

    # Iterate through possible n and f values
    for n in range(n_limit):
        f_n = fib(n)
        for f in range(-25, 26):
            # From F_n + f*F_m + g = 0, we get g = -F_n - f*F_m
            g = -f_n - f * f_m
            
            # Check if g is within the allowed range
            if -25 <= g <= 25:
                solution_count += 1
                # Print the details of the found solution
                print(f"Solution {solution_count}: (n={n}, f={f}, g={g})")
                print(f"  Equation: F_{n} + {f}*F_{m} + {g} = 0")
                print(f"  Values:   {f_n} + {f}*{f_m} + {g} = {f_n + f*f_m + g}")
    
    print("-" * 50)
    print(f"Total number of solutions found for m={m}: {solution_count}")

solve()
