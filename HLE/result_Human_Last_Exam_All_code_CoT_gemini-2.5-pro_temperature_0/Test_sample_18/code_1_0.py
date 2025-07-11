import numpy as np

def solve_k_matrix():
    """
    Calculates the K-matrix of a fractional state based on a series of transformations.
    """
    # Step 1: Define the K-matrix for the Bosonic Integer Quantum Hall state at nu=2.
    K_B = np.array([[0, 1],
                    [1, 0]])

    # Step 2: Define the transformation matrix 'A' for inter-species Cooper pairing.
    # This matrix transforms from the composite fermion basis to the (pair, relative) basis.
    A = np.array([[1, 1],
                  [1, -1]])

    # Step 3: Calculate the K-matrix for the composite fermions (K_cf).
    # The transformation from K_cf to K_B is K_B = (1/4) * A @ K_cf @ A.T.
    # The inverse transformation is K_cf = A @ K_B @ A.T.
    K_cf = A @ K_B @ A.T

    # Step 4: Calculate the K-matrix for the original fermions (K_f).
    # This reverses the attachment of two fluxes (K_cf = K_f + 2*I).
    I = np.identity(K_B.shape[0])
    flux_matrix = 2 * I
    K_f = K_cf - flux_matrix

    # Step 5: Print the results step-by-step as requested.
    print("The initial K-matrix for the bosonic state is K_B:")
    print(K_B)
    print("\nTo find the composite fermion K-matrix (K_cf), we invert the pairing transformation:")
    print("K_cf = A @ K_B @ A^T")
    print("K_cf =")
    print(K_cf.astype(int))
    print("\nTo find the K-matrix for the original fermions (K_f), we invert the flux attachment:")
    print("K_f = K_cf - 2*I")
    
    print("\nThe final equation is:")
    
    # Format and print the final matrix equation with each number.
    k_f_00, k_f_01 = K_f[0, 0], K_f[0, 1]
    k_f_10, k_f_11 = K_f[1, 0], K_f[1, 1]
    
    k_cf_00, k_cf_01 = K_cf[0, 0], K_cf[0, 1]
    k_cf_10, k_cf_11 = K_cf[1, 0], K_cf[1, 1]

    flux_00, flux_01 = flux_matrix[0, 0], flux_matrix[0, 1]
    flux_10, flux_11 = flux_matrix[1, 0], flux_matrix[1, 1]

    print(f"[[{k_f_00:.0f}  {k_f_01:.0f}]    [[{k_cf_00:.0f}  {k_cf_01:.0f}]    [[{flux_00:.0f}  {flux_01:.0f}]")
    print(f" [{k_f_10:.0f} {k_f_11:.0f}]] =  [{k_cf_10:.0f} {k_cf_11:.0f}]] -  [{flux_10:.0f} {flux_11:.0f}]]")

    print("\nThus, the K-matrix for the resulting fractional state of fermions is:")
    print(K_f.astype(int))

solve_k_matrix()