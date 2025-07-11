import math

def solve():
    """
    Calculates the largest possible dimension for R/I.
    """
    n = 10
    d = 10000
    
    # The dimension is given by the binomial coefficient C(n + d - 1, n)
    # which is C(10 + 10000 - 1, 10) = C(10009, 10)
    comb_n = n + d - 1
    comb_k = n

    numerator_str = " * ".join(map(str, range(comb_n, comb_n - comb_k, -1)))
    denominator_str = " * ".join(map(str, range(comb_k, 0, -1)))

    print(f"The dimension is C({comb_n}, {comb_k}), which is calculated as:")
    print(f"({numerator_str}) / ({denominator_str})")

    result = math.comb(comb_n, comb_k)
    
    # As requested, printing the numbers that form the final equation
    print("\nFinal equation:")
    print(f"Dimension = C({comb_n}, {comb_k})")
    print(f"C({comb_n}, {comb_k}) = {result}")

solve()
