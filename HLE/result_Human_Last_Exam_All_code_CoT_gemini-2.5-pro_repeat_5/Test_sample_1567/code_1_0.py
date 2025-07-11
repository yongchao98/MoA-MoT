import numpy as np

def solve():
    """
    This function demonstrates that for k=2, it's possible to construct a
    controlled random walk that can be made recurrent in d=3.
    This shows that the answer to the problem must be k < 2.
    """
    d = 3
    print(f"Let d = {d}. We construct two steering matrices M1 and M2.")
    print("A controlled random walk can be made recurrent if the identity matrix I")
    print("is in the convex hull of the steering matrices {M1, M2}.\n")

    # Define a traceless symmetric matrix Delta
    # We choose Delta such that its eigenvalues are between -1 and 1
    # to ensure I-Delta and I+Delta are positive definite.
    delta_val = 0.5
    Delta = np.diag([delta_val, -delta_val, 0])

    # Identity matrix
    I = np.identity(d)

    # Construct two steering matrices
    M1 = I - Delta
    M2 = I + Delta

    print("Matrix M1:")
    print(M1)
    print(f"\nTrace of M1: {np.trace(M1):.1f}")
    print(f"Eigenvalues of M1: {np.linalg.eigvals(M1)}")
    print("All eigenvalues are positive, so M1 is positive definite.\n")

    print("Matrix M2:")
    print(M2)
    print(f"\nTrace of M2: {np.trace(M2):.1f}")
    print(f"Eigenvalues of M2: {np.linalg.eigvals(M2)}")
    print("All eigenvalues are positive, so M2 is positive definite.\n")

    # Show that I is in the convex hull of {M1, M2}
    alpha = 0.5
    convex_combination = alpha * M1 + (1 - alpha) * M2

    print("We can write I as a convex combination of M1 and M2:")
    print(f"{alpha} * M1 + {1 - alpha} * M2 = I")
    
    # Printing the equation with numbers
    print("\n--- The Final Equation ---")
    print(f"{alpha} * \n{M1}\n+ {1 - alpha} * \n{M2}\n= \n{convex_combination}")
    print("\nThis resulting matrix is the identity matrix, fulfilling the condition for recurrence.")
    print("Therefore, for k=2, transience is NOT guaranteed for any choice of measures.")

solve()