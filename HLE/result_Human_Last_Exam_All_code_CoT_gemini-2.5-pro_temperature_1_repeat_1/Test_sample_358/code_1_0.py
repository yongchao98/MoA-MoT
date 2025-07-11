import numpy as np

def solve_cartan_sum():
    """
    Calculates the sum of entries in the Cartan matrix for the principal block of k[A5 x C2] at char 2.
    """
    # Step 3: Define the decomposition matrices for the principal 2-blocks of A5 and C2.
    # D_A5_0 corresponds to the principal block of A5.
    D_A5_0 = np.array([
        [1, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [1, 1, 1]
    ])

    # D_C2_0 corresponds to the principal block of C2.
    D_C2_0 = np.array([
        [1],
        [1]
    ])

    # Step 4.1: Compute the decomposition matrix for the principal block of G = A5 x C2
    # using the Kronecker product.
    D_G = np.kron(D_A5_0, D_C2_0)

    # Step 4.2: Compute the Cartan matrix C_G = D_G^T * D_G.
    C_G = D_G.T @ D_G

    # Step 5: Sum all the entries of the Cartan matrix.
    total_sum = np.sum(C_G)

    # Format the output equation string.
    equation_parts = [str(entry) for entry in C_G.flatten()]
    equation_str = " + ".join(equation_parts)
    
    print("The Cartan matrix for the principal block is:")
    print(C_G)
    print("\nThe sum of its entries is calculated as follows:")
    print(f"{equation_str} = {total_sum}")

solve_cartan_sum()