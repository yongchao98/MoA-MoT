import numpy as np

def solve_k_matrix():
    """
    This function calculates the K-matrix for a fractional quantum Hall state
    derived from a bosonic integer quantum Hall state.
    """
    # Step 1: Define the K-matrix for the initial bosonic state (nu=2 BIQH).
    # K_boson = sigma_x
    K_boson = np.array([[0, 1],
                        [1, 0]], dtype=float)

    # Step 2: The bosons are Cooper pairs (p=2) of composite fermions.
    # To find the K-matrix of the composite fermions (K_cf), we divide by p^2.
    p = 2
    K_cf = K_boson / (p**2)

    # Step 3: The composite fermions are formed by attaching m=2 fluxes to original fermions.
    # To find the K-matrix of the original fermions (K_final), we subtract m*I.
    m = 2
    # The identity matrix must have the same dimension as K_cf.
    identity_matrix = np.identity(K_cf.shape[0])
    K_final = K_cf - m * identity_matrix

    # Step 4: Print the final result in the requested format.
    k11 = K_final[0, 0]
    k12 = K_final[0, 1]
    k21 = K_final[1, 0]
    k22 = K_final[1, 1]

    print("The K-matrix of the resulting fractional state is:")
    print(f"K = [[{k11}, {k12}],")
    print(f"     [{k21}, {k22}]]")

solve_k_matrix()