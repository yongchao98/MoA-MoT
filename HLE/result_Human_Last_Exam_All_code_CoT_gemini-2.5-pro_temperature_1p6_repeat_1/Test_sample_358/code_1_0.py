import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of all entries in the Cartan matrix of the principal
    block of the group algebra k[A5 x C2] for a field k of characteristic 2.
    """
    
    # The Cartan matrix for the principal 2-block of the alternating group A5
    # is a standard result in modular representation theory. It is a 4x4 matrix.
    C_A5 = np.array([
        [4, 2, 2, 2],
        [2, 2, 1, 1],
        [2, 1, 2, 1],
        [2, 1, 1, 2]
    ])

    # For the cyclic group C2 in characteristic 2, the group algebra kC2
    # is the principal block. It has one simple module, and its Cartan matrix is (2).
    C_C2 = np.array([
        [2]
    ])

    # The sum of entries in the Cartan matrix of A5's principal block.
    sum_A5 = np.sum(C_A5)
    
    # The sum of entries in the Cartan matrix of C2's principal block.
    sum_C2 = np.sum(C_C2)

    # The sum of entries of the Kronecker product of two matrices is the
    # product of their individual sums of entries.
    total_sum = sum_A5 * sum_C2

    print("The Cartan matrix for the principal 2-block of A5 is:")
    print(C_A5)
    print(f"The sum of its entries is: {sum_A5}")
    print("-" * 30)
    print("The Cartan matrix for the principal 2-block of C2 is:")
    print(C_C2)
    print(f"The sum of its entries is: {sum_C2}")
    print("-" * 30)
    print("The total sum of entries for the Cartan matrix of A5 x C2 is the product of the individual sums.")
    print(f"Total Sum = {sum_A5} * {sum_C2} = {total_sum}")

solve_cartan_sum()