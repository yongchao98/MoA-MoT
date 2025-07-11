import numpy as np

def construct_and_verify():
    """
    This function constructs a 4x4 table with distinct positive integers
    such that the product of numbers in each row and column is the same.
    It finds the smallest N, which is the maximum number in this table.
    """

    # Exponent Matrix E1. Entries are {0..7} each twice. Magic Sum = 14.
    E1 = np.array([
        [0, 7, 1, 6],
        [2, 5, 3, 4],
        [7, 0, 6, 1],
        [5, 2, 4, 3]
    ])

    # Exponent Matrix E2. Entries are {0,1} eight times each. Magic Sum = 2.
    E2 = np.array([
        [1, 0, 1, 0],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [0, 1, 0, 1]
    ])

    # Generate the table of numbers A_ij = 2^E1_ij * 3^E2_ij
    A = np.zeros((4, 4), dtype=int)
    for i in range(4):
        for j in range(4):
            A[i, j] = (2**E1[i, j]) * (3**E2[i, j])

    # Find the maximum value N in the table
    N = np.max(A)

    print(f"The constructed 4x4 table is:")
    for row in A:
        print(" ".join(f"{num:4d}" for num in row))
    print(f"\nThe largest number in the table is N = {N}.")

    # Verify the properties
    print("\nVerifying properties:")
    
    # Check for distinct integers
    is_distinct = len(np.unique(A)) == 16
    print(f"1. All integers are distinct: {is_distinct}")

    # Check row and column products
    row_products = [np.prod(A[i, :]) for i in range(4)]
    col_products = [np.prod(A[:, j]) for j in range(4)]
    
    P = row_products[0]
    all_products_equal = all(p == P for p in row_products) and all(p == P for p in col_products)
    
    print(f"2. The common product P is: {P}")
    print(f"   All row and column products are equal to P: {all_products_equal}")
    
    print("\nRow products:")
    for i, p in enumerate(row_products):
        print(f"   Row {i+1}: {' * '.join(map(str, A[i, :]))} = {p}")
        
    print("\nColumn products:")
    for j, p in enumerate(col_products):
        col_str = ' * '.join(map(str, A[:, j]))
        print(f"   Col {j+1}: {col_str} = {p}")

    print(f"\nThe smallest N is {N}.")


construct_and_verify()