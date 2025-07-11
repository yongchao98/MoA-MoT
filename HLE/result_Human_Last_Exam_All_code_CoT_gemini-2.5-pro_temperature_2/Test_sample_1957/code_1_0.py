import math

def solve():
    """
    Solves for the minimum value of 1000m + n based on the problem description.
    """
    # Based on the reasoning in the thought block, the degree 'd' of the polynomial
    # det(F(A,B,C,D,X)) is determined to be 3.
    # d = deg(P_S) = deg(P_cone) = 3
    d = 3

    print(f"The degree of the polynomial det(F) is d = m*n = {d}.")
    print("We need to find positive integers m and n that are factors of d.")
    
    # Find all pairs of factors (m, n) for d
    factors = []
    for i in range(1, int(math.sqrt(d)) + 1):
        if d % i == 0:
            factors.append((i, d // i))
            if i*i != d:
                factors.append((d // i, i))
    
    factors.sort()

    min_val = float('inf')
    best_pair = (0, 0)
    
    print("The possible pairs (m, n) and the corresponding values of 1000*m + n are:")
    for m, n in factors:
        value = 1000 * m + n
        print(f"m = {m}, n = {n}  => 1000 * {m} + {n} = {value}")
        if value < min_val:
            min_val = value
            best_pair = (m, n)

    print(f"\nThe minimum value is {min_val}, which occurs for m={best_pair[0]} and n={best_pair[1]}.")
    
solve()