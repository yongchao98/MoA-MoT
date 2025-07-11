import math

def solve_continuum_decomposition():
    """
    Calculates the largest number n for the given topology problem.
    
    The problem describes a space X which is a "k-junction" for k=5 points.
    The largest number n of subcontinua in an irreducible decomposition of such a space
    is given by the combinatorial formula: n = C(k-1, floor((k-1)/2)).
    """
    
    # The number of special points in the continuum
    k = 5
    
    # According to the theorem, the value n is C(k-1, floor((k-1)/2))
    n_val = k - 1
    k_val = math.floor((k - 1) / 2)
    
    # Calculate the binomial coefficient C(n_val, k_val)
    result = math.comb(n_val, k_val)
    
    # Output the explanation and the final equation step by step
    print(f"The number of special points is k = {k}.")
    print("The largest number n is given by the formula C(k-1, floor((k-1)/2)).")
    print(f"For k = {k}, this evaluates to C({n_val}, floor(({k}-1)/2)) = C({n_val}, {k_val}).")
    print(f"The final equation is C({n_val}, {k_val}) = {result}.")

solve_continuum_decomposition()
