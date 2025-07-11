def solve():
    """
    This function finds the minimum value of 1000m + n based on the problem description.
    """
    # The geometric condition corresponds to the vanishing of a single irreducible polynomial P.
    # The degree of this polynomial, d, is 6.
    # We need to find positive integers m and n such that det(F) is a polynomial of degree d,
    # where F is an n x n matrix of polynomials of degree m.
    # This gives the constraint n * m = d.
    d = 6

    # We need to minimize 1000m + n subject to m*n = d, where m and n are positive integers.
    # We can list all possible pairs (m, n) since d is small.
    possible_pairs = []
    for m in range(1, d + 1):
        if d % m == 0:
            n = d // m
            possible_pairs.append((m, n))

    min_val = float('inf')
    best_m = 0
    best_n = 0

    # Calculate 1000m + n for each pair and find the minimum.
    for m, n in possible_pairs:
        val = 1000 * m + n
        if val < min_val:
            min_val = val
            best_m = m
            best_n = n
    
    # The user wants the final equation printed.
    print(f"The degree of the polynomial defining the variety is d = {d}.")
    print("We need to minimize 1000m + n, subject to m * n = d.")
    print("The possible pairs (m, n) and the corresponding values are:")
    for m, n in possible_pairs:
        print(f"  (m={m}, n={n}) -> 1000 * {m} + {n} = {1000*m + n}")
    
    print("\nThe minimum value is obtained with the pair (m, n) that gives the smallest result.")
    print(f"The optimal values are m = {best_m} and n = {best_n}.")
    print(f"The final calculation is: {1000} * {best_m} + {best_n} = {min_val}")

solve()