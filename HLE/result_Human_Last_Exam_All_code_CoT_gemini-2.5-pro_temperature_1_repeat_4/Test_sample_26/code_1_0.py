import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve_pi3_rank():
    """
    Calculates the rank of the third homotopy group of a smooth quintic hypersurface in CP^3.
    """
    # n: dimension of the complex projective space
    n = 3
    # d: degree of the hypersurface
    d = 5

    print(f"The rank of pi_3(X) for a simply connected complex surface X is given by the Hodge number h^(2,0)(X).")
    print(f"For a hypersurface X of degree d in CP^n, h^(2,0)(X) = dim H^0(X, K_X), where K_X is the canonical bundle.")
    print(f"By the adjunction formula, K_X = O_X(d - n - 1).")
    
    # k: the degree of the canonical bundle
    k = d - n - 1
    
    print(f"\nGiven values:")
    print(f"Dimension of projective space, n = {n}")
    print(f"Degree of hypersurface, d = {d}")
    
    print(f"\nStep 1: Calculate the degree of the canonical bundle O_X(k).")
    print(f"k = d - n - 1 = {d} - {n} - 1 = {k}")
    
    print(f"\nStep 2: Calculate h^(2,0)(X) = dim H^0(CP^n, O(k)).")
    print(f"The dimension is given by the binomial coefficient C(n + k, k).")
    
    # Calculate the rank
    rank = combinations(n + k, k)
    
    print(f"\nFinal calculation:")
    print(f"rank = C({n} + {k}, {k}) = C({n+k}, {k}) = {rank}")

solve_pi3_rank()