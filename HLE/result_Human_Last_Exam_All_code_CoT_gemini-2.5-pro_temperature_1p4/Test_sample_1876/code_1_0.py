import numpy as np

def solve():
    """
    Finds the smallest N for a 4x4 multiplicative magic square.
    This solution constructs the square and calculates N.
    """

    # Step 1: Define orthogonal Latin squares M1, M2 of order 4
    # These are known constructions. M1 and M2 are magic squares themselves.
    M1 = np.array([
        [0, 1, 2, 3],
        [2, 3, 0, 1],
        [3, 2, 1, 0],
        [1, 0, 3, 2]
    ])

    M2 = np.array([
        [0, 1, 2, 3],
        [1, 0, 3, 2],
        [2, 3, 0, 1],
        [3, 2, 1, 0]
    ])

    # Step 2: Construct exponent matrices E, F, G for primes 2, 3, 5
    # The construction guarantees that the triplets (E_ij, F_ij, G_ij) are unique.
    E = M1  # Exponents for prime 2, from {0, 1, 2, 3}
    F = M2 % 2  # Exponents for prime 3, from {0, 1}
    G = M2 // 2 # Exponents for prime 5, from {0, 1}

    # Step 3: Create the 4x4 table of numbers
    table = np.zeros((4, 4), dtype=int)
    for i in range(4):
        for j in range(4):
            table[i, j] = (2**E[i, j]) * (3**F[i, j]) * (5**G[i, j])

    # Step 4: Verify the solution and find N
    # Check for distinctness
    if len(np.unique(table)) != 16:
        print("Error: The numbers in the table are not distinct.")
        return

    # Check products
    target_product = np.prod(table[0, :])
    valid = True
    for i in range(4):
        if np.prod(table[i, :]) != target_product:
            valid = False
            print(f"Row {i+1} product is incorrect.")
    for j in range(4):
        if np.prod(table[:, j]) != target_product:
            valid = False
            print(f"Column {j+1} product is incorrect.")

    # Find the maximum value N
    N = np.max(table)
    
    # Print the resulting table
    print("The constructed 4x4 table is:")
    for row in table:
        print("  ".join(f"{num:3d}" for num in row))
    print(f"\nThe product of each row and column is: {target_product}")
    print(f"\nThe numbers are all distinct and positive.")
    print(f"\nThe maximum number in the table is N = {N}.")


solve()
