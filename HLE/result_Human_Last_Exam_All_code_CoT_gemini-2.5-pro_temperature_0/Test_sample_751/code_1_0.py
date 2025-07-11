import numpy as np

def solve_observer_canonical_form():
    """
    Reduces a discrete-time system to observer canonical form using duality
    and calculates the new C matrix.
    """
    # Given system matrices
    A = np.array([[1, 1, 0],
                  [2, 1, 1],
                  [0, 2, 0]])
    C = np.array([[0, 1, 0],
                  [1, 1, 0]])

    # Step 1: Define the dual system (A_d = A^T, B_d = C^T)
    A_d = A.T
    B_d = C.T

    # Step 2: Use the first input of the dual system, b1, which corresponds to the first row of C.
    # The system (A_d, b1) is controllable, which simplifies the problem.
    b1 = B_d[:, 0].reshape(-1, 1)

    # Construct the controllability matrix for the single-input dual system (A_d, b1)
    # Wc = [b1, A_d*b1, A_d^2*b1]
    Wc = np.hstack([b1, A_d @ b1, A_d @ A_d @ b1])

    # Step 3: Find the transformation matrix T for the dual system's controller canonical form.
    # This requires the coefficients of the characteristic polynomial of A (or A_d).
    # The polynomial is p(z) = z^3 + a1*z^2 + a2*z + a3
    # np.poly returns coefficients for p[0]*z^n + p[1]*z^(n-1) + ...
    coeffs = np.poly(A)
    a1 = coeffs[1]
    a2 = coeffs[2]
    a3 = coeffs[3]

    # Construct the matrix M used in the transformation T = Wc * M
    M = np.array([[a2, a1, 1],
                  [a1, 1,  0],
                  [1,  0,  0]])

    # Calculate the transformation matrix T
    T = Wc @ M

    # Step 4: The transformation for the original system is P = T^-1
    P = np.linalg.inv(T)

    # Step 5: Calculate the new C matrix in the observer canonical form
    C_new = C @ P

    # Print the final result
    print("The original matrix A is:")
    print(A)
    print("\nThe original matrix C is:")
    print(C)
    print("\nThe new matrix C in observer canonical form is:")
    # We need to print each number in the final equation.
    # Printing the numpy array achieves this.
    print(C_new)

solve_observer_canonical_form()