import math

def solve():
    """
    Calculates the value of av_{333}^3(1324) based on the combinatorial pattern
    that for large n, av_n^k(1324) = C_{k+1}, where C is the Catalan number.
    """
    n = 333
    k = 3
    
    # For large n and small k, av_n^k(1324) is conjectured to be C_{k+1}.
    # We need to calculate C_{k+1} = C_4.
    catalan_index = k + 1
    
    # Formula for Catalan numbers: C_m = (1/(m+1)) * binomial(2m, m)
    m = catalan_index
    
    # Calculate binomial coefficient (2m choose m)
    # Here, 2m = 2 * 4 = 8, and m = 4
    term1 = 2 * m
    term2 = m
    binomial_coefficient = math.comb(term1, term2)
    
    # Calculate C_m
    result = binomial_coefficient // (m + 1)
    
    print(f"The number of 1324-avoiding permutations of length n={n} with k={k} inversions, av_{n}^{k}(1324),")
    print(f"for large n is given by the Catalan number C_{k+1} = C_{catalan_index}.")
    print("\nCalculation:")
    print(f"C_{m} = (1 / ({m} + 1)) * binomial(2 * {m}, {m})")
    print(f"C_{m} = (1 / {m + 1}) * binomial({term1}, {term2})")
    print(f"C_{m} = (1 / {m + 1}) * {binomial_coefficient}")
    print(f"C_{m} = {result}")
    
    print(f"\nThus, av_{n}^{k}(1324) = {result}")

solve()