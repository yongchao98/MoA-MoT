import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of all entries in the Cartan matrix of the principal
    block of the group algebra kG, where G = A5 x C2 and char(k)=2.
    """

    # Part 1: Analyze the principal block of A_5 at p=2.

    # The decomposition matrix for the principal 2-block of A5 is known.
    # Rows correspond to ordinary characters (dims 1, 3, 3, 5).
    # Columns correspond to simple modules (dims 1, 2, 2).
    D_A5 = np.array([
        [1, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [1, 1, 1]
    ])

    # The Cartan matrix C is calculated as D_transpose * D.
    C_A5 = D_A5.T @ D_A5

    # The sum of entries of the Cartan matrix for A5's principal block.
    sum_A5 = np.sum(C_A5)

    # Part 2: Analyze the principal block of C_2 at p=2.

    # The group algebra kC_2 is a 2-dimensional local algebra.
    # Its Cartan matrix is simply [[2]].
    sum_C2 = 2

    # Part 3: Combine the results for G = A5 x C2.

    # The sum of entries of the Cartan matrix for the tensor product of blocks
    # is the product of the individual sums.
    total_sum = sum_A5 * sum_C2

    print("The sum of entries for the Cartan matrix of the principal block of A5 is calculated as follows:")
    print(f"Decomposition Matrix D_A5:\n{D_A5}")
    print(f"Cartan Matrix C_A5 = D_A5^T * D_A5:\n{C_A5}")
    print(f"Sum of entries for A5 block = {int(sum_A5)}")
    print("-" * 20)
    print("The sum of entries for the Cartan matrix of the principal block of C2 is 2.")
    print("-" * 20)
    print("The total sum for the principal block of A5 x C2 is the product of these sums.")
    print(f"Final Equation: {int(sum_A5)} * {sum_C2} = {int(total_sum)}")

solve_cartan_sum()