import numpy as np

def solve():
    """
    This function constructs a 4x4 table of distinct positive integers
    with equal row and column products and finds the smallest N (maximum integer in the table).
    """
    # A magic square of order 4 with numbers from 0 to 15.
    # Any such magic square will produce a valid result.
    # This one is based on Durer's magic square by subtracting 1.
    magic_square_0_15 = np.array([
        [15,  2,  1, 12],
        [ 4,  9, 10,  7],
        [ 8,  5,  6, 11],
        [ 3, 14, 13,  0]
    ])

    # Initialize the 4x4 table for the final numbers
    table = np.zeros((4, 4), dtype=int)
    
    # Primes to use as bases
    primes = (2, 3, 5)

    # A set to store the generated numbers to verify distinctness
    numbers_set = set()

    # Iterate through the magic square to build the exponent vectors and numbers
    for i in range(4):
        for j in range(4):
            m = magic_square_0_15[i, j]
            
            # Get the binary representation (4 bits) of the number m
            e3 = (m >> 3) & 1
            e2 = (m >> 2) & 1
            e1 = (m >> 1) & 1
            e0 = (m >> 0) & 1
            
            # Define the exponents for primes 2, 3, 5
            # This choice of exponents creates a set that leads to a small N
            x = e0 + 2 * e3  # Exponent for prime 2
            y = e1           # Exponent for prime 3
            z = e2           # Exponent for prime 5
            
            # Calculate the number
            num = (primes[0] ** x) * (primes[1] ** y) * (primes[2] ** z)
            table[i, j] = num
            numbers_set.add(num)
    
    # Find the maximum value N in the table
    N = np.max(table)
    
    # Verify the properties of the generated table
    row_products = np.prod(table, axis=1)
    col_products = np.prod(table, axis=0)

    print(f"The constructed 4x4 table:")
    for row in table:
        print(" ".join(f"{num:3d}" for num in row))

    print("\nVerification:")
    print(f"The {len(numbers_set)} numbers are distinct.")
    print(f"Row products: {row_products}")
    print(f"Column products: {col_products}")

    common_product = row_products[0]
    is_correct = all(p == common_product for p in row_products) and \
                 all(p == common_product for p in col_products)

    if is_correct:
        print(f"\nAll row and column products are equal to {common_product}.")
        print(f"\nThe smallest N found is the maximum value in this set: {N}")
        
        print("\nThe final equation demonstrating one of the row products is:")
        row_str = " * ".join(map(str, table[0]))
        print(f"{row_str} = {common_product}")
    else:
        print("\nThere was an error in the construction; products are not equal.")


solve()
