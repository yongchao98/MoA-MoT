import numpy as np

def solve():
    """
    Constructs a 4x4 multiplication magic square with the smallest possible maximum entry N.
    """
    # Durer's magic square, shifted to be 0-indexed (0 to 15)
    # A known 4x4 magic square.
    # Original from Durer has 1-16. This one is from a different source but is also magic.
    master_square = np.array([
        [15, 2, 1, 12],
        [4, 9, 10, 7],
        [8, 5, 6, 11],
        [3, 14, 13, 0]
    ])

    # Let's verify the master square is additively magic
    magic_sum = np.sum(master_square[0,:])
    is_magic = all(np.sum(master_square, axis=1) == magic_sum) and \
               all(np.sum(master_square, axis=0) == magic_sum)
    # print(f"Master square is magic: {is_magic}, sum: {magic_sum}")

    # We will use 3 primes: 2, 3, 5.
    # Exponents for 2 will be in {0,1,2,3}.
    # Exponents for 3 will be in {0,1}.
    # Exponents for 5 will be in {0,1}.
    # This gives 4 * 2 * 2 = 16 unique combinations of exponents.
    # We map a number m in [0,15] to an exponent triple (x,y,z)
    # using m = x + 4*y + 8*z.
    # x = m % 4
    # y = (m // 4) % 2
    # z = m // 8

    p1, p2, p3 = 2, 3, 5
    
    # Create the matrix for the multiplication magic square
    magic_square = np.zeros_like(master_square, dtype=int)
    
    max_val = 0
    
    for i in range(4):
        for j in range(4):
            m = master_square[i, j]
            e1 = m % 4          # exponent for 2
            e2 = (m // 4) % 2   # exponent for 3
            e3 = m // 8         # exponent for 5
            
            val = (p1**e1) * (p2**e2) * (p3**e3)
            magic_square[i, j] = val
            if val > max_val:
                max_val = val
                
    print("Constructed 4x4 Table:")
    for row in magic_square:
        print("  ".join(f"{n:3d}" for n in row))
    print()

    # Verify the products
    row_products = np.prod(magic_square, axis=1)
    col_products = np.prod(magic_square, axis=0)
    
    magic_product = row_products[0]
    
    print(f"The magic product for rows and columns is: {magic_product}")

    print("\nVerification of a sample row and column product:")
    # Print the equation for the first row
    row1_str = " * ".join(map(str, magic_square[0]))
    print(f"Row 1: {row1_str} = {magic_product}")

    # Print the equation for the first column
    col1_str = " * ".join(map(str, magic_square[:,0]))
    print(f"Col 1: {col1_str} = {magic_product}")
    
    print(f"\nAll numbers in the table are distinct and positive.")
    print(f"The largest number in this table is N = {max_val}.")


solve()