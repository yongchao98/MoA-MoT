import math

def solve():
    """
    This function constructs a 4x4 multiplicative magic square with the smallest
    possible maximum entry N, and verifies its properties.
    """
    # Step 1: Define the two sets of numbers for the construction.
    # These sets are chosen to minimize the maximum element in the resulting square.
    Sa = [1, 2, 3, 4]
    Sb = [1, 5, 6, 7]

    # Step 2: Define the Latin squares for the indices of Sa and Sb.
    # These ensure that each row and column product is the same.
    a_indices = [
        [1, 2, 3, 4],
        [2, 1, 4, 3],
        [3, 4, 1, 2],
        [4, 3, 2, 1]
    ]

    b_indices = [
        [1, 2, 3, 4],
        [3, 4, 1, 2],
        [4, 3, 2, 1],
        [2, 1, 4, 3]
    ]

    # Step 3: Construct the 4x4 magic square.
    matrix = [[0 for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            # The elements are products of numbers from Sa and Sb based on the indices
            matrix[i][j] = Sa[a_indices[i][j] - 1] * Sb[b_indices[i][j] - 1]

    print("Constructed 4x4 Multiplicative Magic Square:")
    for row in matrix:
        print("  ".join(map(str, row)))

    print("\nVerifying the magic property:")
    
    # Step 4: Calculate the magic product P.
    magic_product = 1
    for x in Sa:
        magic_product *= x
    for y in Sb:
        magic_product *= y

    # Step 5: Verify row products and print the equations.
    print("\nRow products:")
    for i in range(4):
        product = 1
        for j in range(4):
            product *= matrix[i][j]
        
        equation = " * ".join(map(str, matrix[i]))
        print(f"Row {i+1}: {equation} = {product}")
        if product != magic_product:
            print(f"Verification failed for row {i+1}")

    # Step 6: Verify column products and print the equations.
    print("\nColumn products:")
    for j in range(4):
        product = 1
        col_elements = []
        for i in range(4):
            product *= matrix[i][j]
            col_elements.append(matrix[i][j])
        
        equation = " * ".join(map(str, col_elements))
        print(f"Column {j+1}: {equation} = {product}")
        if product != magic_product:
            print(f"Verification failed for column {j+1}")
            
    # Step 7: Find the maximum value N in the matrix.
    max_n = 0
    all_elements = []
    for i in range(4):
        for j in range(4):
            all_elements.append(matrix[i][j])
            if matrix[i][j] > max_n:
                max_n = matrix[i][j]

    print(f"\nThe 16 distinct integers are: {sorted(all_elements)}")
    print(f"\nThe smallest N is the maximum element in this square, which is {max_n}.")

solve()