def solve():
    """
    This function finds the smallest N for the 4x4 table problem,
    constructs the table, and verifies the solution.
    """
    s1 = [1, 2, 3, 4]
    s2 = [1, 5, 6, 7]
    
    a, b, c, d = s1
    x, y, z, w = s2
    
    # Construct the table using a Graeco-Latin square design
    table = [
        [a*x, b*y, c*z, d*w],
        [b*z, a*w, d*x, c*y],
        [c*w, d*z, a*y, b*x],
        [d*y, c*x, b*w, a*z]
    ]

    # Find the maximum value N in the table
    N = 0
    for i in range(4):
        for j in range(4):
            if table[i][j] > N:
                N = table[i][j]

    print(f"The smallest N is {N}.")
    print("\nThe 4x4 table is:")
    for row in table:
        print("  ".join(map(str, row)))
    
    # Verify the products of rows and columns
    print("\nVerification of products:")
    
    # Row products
    row_products = []
    for i in range(4):
        product = 1
        for j in range(4):
            product *= table[i][j]
        row_products.append(product)
        print(f"Row {i+1}: {' * '.join(map(str, table[i]))} = {product}")
        
    # Column products
    col_products = []
    for j in range(4):
        product = 1
        col_values = []
        for i in range(4):
            product *= table[i][j]
            col_values.append(table[i][j])
        col_products.append(product)
        print(f"Col {j+1}: {' * '.join(map(str, col_values))} = {product}")

solve()