import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of entries in the Cartan matrix for the principal block
    of the group algebra k(A_5 x C_2) where char(k) = 2.
    """

    # Step 1: Define the Cartan matrix for the principal block of A_5 at p=2.
    # The ordinary characters in the principal block of A_5 are of degrees 1, 3, 3, 5.
    # The simple modules in this block have dimensions 1, 2, 2.
    # The decomposition matrix D, which expresses ordinary characters as sums of Brauer
    # characters (simple modules), is known to be:
    #       L1 L2 L3
    # chi1  [1, 0, 0]
    # chi2  [1, 1, 0]
    # chi3  [1, 0, 1]
    # chi5  [1, 1, 1]
    # The Cartan matrix C_A5 is the product of the transpose of D with D.
    
    D_A5 = np.array([
        [1, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [1, 1, 1]
    ])
    C_A5 = D_A5.T @ D_A5

    # Step 2: Define the Cartan matrix for the principal block of C_2 at p=2.
    # kC_2 has one simple module (trivial). The projective cover of the trivial
    # module is kC_2 itself, which has composition length 2.
    # Thus, the Cartan matrix is a 1x1 matrix with entry 2.
    C_C2 = np.array([[2]])

    # Step 3: The Cartan matrix for the direct product G = A_5 x C_2 is the
    # Kronecker product of the individual Cartan matrices.
    C_G = np.kron(C_A5, C_C2)

    # Step 4: Calculate the sum of all entries in the resulting matrix.
    total_sum = int(np.sum(C_G))

    # Step 5: Print the calculation and the result as requested.
    print("The Cartan matrix for the principal block of k(A_5 x C_2) is:")
    print(C_G)
    
    elements = C_G.flatten()
    equation_str = " + ".join(map(str, elements))
    print(f"\nThe sum of all the entries in this matrix is:\n{equation_str} = {total_sum}")

solve_cartan_sum()