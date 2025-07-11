import numpy as np

def solve_and_print_pole_placement():
    """
    Calculates the state feedback gain F for a given system so that
    the closed-loop system A + BF has specified eigenvalues.
    """
    # Step 0: Define the system matrices and desired eigenvalues (poles)
    A = np.array([[-1, 1], [1, 0]])
    B = np.array([[1, 2], [1, 0]])
    desired_poles = np.array([-1 + 1j, -1 - 1j])

    # Step 1: Determine the desired characteristic polynomial
    # For poles p1, p2, the polynomial is (s-p1)(s-p2) = s^2 - (p1+p2)s + p1*p2
    # Sum of poles: (-1 + 1j) + (-1 - 1j) = -2
    # Product of poles: (-1 + 1j) * (-1 - 1j) = (-1)^2 - (1j)^2 = 1 - (-1) = 2
    # Desired characteristic polynomial: s^2 + 2s + 2 = 0
    poly_coeffs = [1, 2, 2] # Coefficients for s^2, s^1, s^0

    # Step 2: Construct a target closed-loop matrix A_cl (in companion form)
    # The controller companion form for s^2 + a1*s + a0 = 0 is [[0, 1], [-a0, -a1]]
    a1 = poly_coeffs[1]
    a0 = poly_coeffs[2]
    A_cl = np.array([[0, 1], [-a0, -a1]])

    # Step 3: Calculate the feedback gain F using F = B_inv * (A_cl - A)
    try:
        B_inv = np.linalg.inv(B)
    except np.linalg.LinAlgError:
        print("Matrix B is singular, cannot compute a unique F with this method.")
        return

    F = B_inv @ (A_cl - A)

    # Step 4: Print the results and verification
    print("Given system matrices:")
    print(f"A =\n{A}")
    print(f"B =\n{B}\n")

    print("Desired eigenvalues: [-1 + j, -1 - j]\n")

    print("Calculated state feedback gain matrix F:")
    print(f"F =\n{np.round(F, 4)}\n")
    
    # Verify the result by calculating A + B*F
    result_A_cl = A + B @ F

    print("The final equation A + B * F is:")
    print(f"{A} +\n\n{B} @\n\n{np.round(F, 4)} =\n\n{np.round(result_A_cl, 4)}\n")

    # Final check: Calculate eigenvalues of the resulting closed-loop system
    eigenvalues = np.linalg.eigvals(result_A_cl)
    print("Verification: Eigenvalues of the resulting matrix A + BF are:")
    print(np.round(eigenvalues, 4))

solve_and_print_pole_placement()
<<<[[-3.0, -2.0], [2.0, 1.0]]>>>