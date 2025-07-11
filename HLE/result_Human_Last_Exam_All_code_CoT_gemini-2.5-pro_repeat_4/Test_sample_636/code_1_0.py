import numpy as np

def solve_pole_placement():
    """
    Calculates the state feedback gain F for a given system (A, B)
    to place the closed-loop poles at desired locations.
    """
    # Given system matrices
    A = np.array([[-1, 1],
                  [1, 0]])
    B = np.array([[1, 2],
                  [1, 0]])

    # Desired eigenvalues (poles)
    p = np.array([-1 + 1j, -1 - 1j])

    # The desired characteristic polynomial is (s-p1)(s-p2) = s^2 - (p1+p2)s + p1*p2
    # For p = -1 +- j, the polynomial is s^2 + 2s + 2 = 0
    # The coefficients are [1, 2, 2]
    desired_coeffs = np.poly(p)
    # alpha1 is the coefficient of s, alpha0 is the constant term
    alpha1 = desired_coeffs[1]
    alpha0 = desired_coeffs[2]

    # The characteristic polynomial of A_cl = A + B*F is s^2 - trace(A_cl)s + det(A_cl) = 0
    # This gives: s^2 + (1 - f11 - 2*f21 - f12)s + (-1 - f11 - 2*f12 + 2*f21*f12 - 2*f22 - 2*f11*f22) = 0
    #
    # Equating coefficients with s^2 + alpha1*s + alpha0 = 0 gives two equations for the four unknowns in F.
    # The system is underdetermined. We can find one solution by setting f21=0 and f22=0.
    # The equations simplify to:
    # 1 - f11 - f12 = alpha1
    # -1 - f11 - 2*f12 = alpha0
    #
    # Rewriting as a linear system M * f = v:
    # f11 + f12 = 1 - alpha1
    # f11 + 2*f12 = -1 - alpha0
    #
    # where f = [f11, f12]^T

    M = np.array([[1, 1],
                  [1, 2]])

    v = np.array([1 - alpha1,
                  -1 - alpha0])

    # Solve for f11 and f12
    f_solved = np.linalg.solve(M, v)
    f11 = f_solved[0]
    f12 = f_solved[1]

    # Construct the full feedback gain matrix F
    F = np.array([[f11, f12],
                  [0.0, 0.0]])

    print("Calculated State Feedback Gain Matrix F:")
    print(F)

    # Verification
    A_cl = A + B @ F
    eigenvalues = np.linalg.eigvals(A_cl)

    print("\nVerification of the solution:")
    print("The final equation is A + B @ F = A_cl")
    print(f"{A} + \n{B} @ \n{F} = \n{A_cl}")
    
    print("\nEigenvalues of the closed-loop system A + BF:")
    # Using np.round to clean up floating point inaccuracies
    print(np.round(eigenvalues, 5))

    print("\nDesired eigenvalues:")
    print(p)

solve_pole_placement()
<<<[[1.0, -2.0], [0.0, 0.0]]>>>