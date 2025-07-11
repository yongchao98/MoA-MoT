import math

def solve_problem():
    """
    Solves for the largest p such that the function I is not in L^p(R^9).
    """

    # The dimension of the parameter space for the function I(a_1, ..., a_9)
    n = 9
    print(f"The dimension of the parameter space is n = {n}.")

    # The monomials in the phase of the integral correspond to the exponents (k, l)
    # P(x,y) = a1*x + a2*y + a3*x^2 + a4*xy + ... + a9*y^3
    monomial_exponents = {
        'a1 (x)': (1, 0),
        'a2 (y)': (0, 1),
        'a3 (x^2)': (2, 0),
        'a4 (xy)': (1, 1),
        'a5 (y^2)': (0, 2),
        'a6 (x^3)': (3, 0),
        'a7 (x^2*y)': (2, 1),
        'a8 (x*y^2)': (1, 2),
        'a9 (y^3)': (0, 3),
    }

    decay_exponents = {}
    print("\nCalculating decay exponents sigma for each monomial term x^k * y^l:")

    for name, (k, l) in monomial_exponents.items():
        if k == 0 and l == 0:
            # Should not happen for this problem
            sigma = float('inf')
        elif k == 0:
            sigma = 1.0 / l
        elif l == 0:
            sigma = 1.0 / k
        else:
            # For a monomial x^k*y^l, the decay exponent is 1/max(k,l)
            sigma = 1.0 / max(k, l)
        
        decay_exponents[name] = sigma
        print(f"For {name}: sigma = 1/max({k},{l}) = {sigma:.4f}")

    # The slowest decay rate is the minimum of these exponents
    sigma_min = min(decay_exponents.values())
    
    # Find which monomial gives this slowest decay
    slowest_terms = [name for name, sigma in decay_exponents.items() if sigma == sigma_min]
    print(f"\nThe minimum decay exponent (slowest decay) is sigma_min = {sigma_min:.4f}, coming from term(s): {', '.join(slowest_terms)}.")
    
    # The condition for the integral of |I|^p over R^n to diverge is p * sigma <= n.
    # The critical case, defining the boundary of integrability, is p * sigma_min = n.
    print(f"\nThe largest p for which I is not in L^p is found by solving the equation:")
    print(f"p * sigma_min = n")
    
    # Solving for p
    p_critical = n / sigma_min
    
    # The final equation with numbers plugged in
    sigma_min_k = 3 # The denominator for the minimal sigma
    print(f"p * (1/{sigma_min_k}) = {n}")
    print(f"p / {sigma_min_k} = {n}")
    print(f"p = {n} * {sigma_min_k}")
    
    print(f"\nThe solution is p = {p_critical}.")

solve_problem()