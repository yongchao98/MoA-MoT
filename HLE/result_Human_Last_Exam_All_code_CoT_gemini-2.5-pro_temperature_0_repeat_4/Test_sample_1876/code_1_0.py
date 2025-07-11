import numpy as np

def solve_magic_product_square():
    """
    Constructs a 4x4 multiplicative semi-magic square to find the smallest N.
    The construction is based on creating two orthogonal additive semi-magic squares
    for the exponents of primes 2 and 3.
    """
    # Exponent matrix for the prime 2.
    # This is a semi-magic square with sum 6.
    X = np.array([
        [0, 1, 2, 3],
        [2, 3, 0, 1],
        [3, 2, 1, 0],
        [1, 0, 3, 2]
    ])

    # Exponent matrix for the prime 3.
    # This is also a semi-magic square with sum 6.
    # The pair of matrices (X, Y) is "orthogonal" meaning they produce 16 unique pairs (x,y).
    Y = np.array([
        [0, 1, 2, 3],
        [1, 0, 3, 2],
        [2, 3, 0, 1],
        [3, 2, 1, 0]
    ])

    # Construct the 4x4 table with distinct integers
    # A[i,j] = 2^X[i,j] * 3^Y[i,j]
    A = np.zeros((4, 4), dtype=int)
    for i in range(4):
        for j in range(4):
            A[i, j] = (2**X[i, j]) * (3**Y[i, j])

    print("Constructed 4x4 Table:")
    for row in A:
        print(" ".join(f"{num:4}" for num in row))
    print("-" * 20)

    # Verify the products of rows and columns
    print("Verifying products:")
    # The theoretical product K = 2^6 * 3^6 = 46656
    is_correct = True
    target_product = np.prod(A[0])
    print(f"Target Product (from first row): {target_product}")

    for i in range(4):
        row_product = np.prod(A[i, :])
        print(f"Row {i+1} product: {row_product}")
        if row_product != target_product:
            is_correct = False

    for j in range(4):
        col_product = np.prod(A[:, j])
        print(f"Column {j+1} product: {col_product}")
        if col_product != target_product:
            is_correct = False
    
    print("-" * 20)
    if is_correct:
        print("All row and column products are the same.")
    else:
        print("Error: Products are not consistent.")

    # Find the largest integer N in the table
    N = np.max(A)
    
    # The set of numbers is {2^x * 3^y | x,y in {0,1,2,3}}
    # The largest number N corresponds to the largest exponents, x=3 and y=3.
    # So, N = 2^3 * 3^3
    print("\nThe smallest N is the maximum value in this table.")
    print(f"N = 2^{np.max(X)} * 3^{np.max(Y)} = 2^3 * 3^3 = {2**3} * {3**3} = {N}")

solve_magic_product_square()