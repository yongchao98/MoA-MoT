import numpy as np

def solve_steady_state_probabilities():
    """
    Solves for the steady-state probabilities of the given Kolmogorov-Chepmen system
    and calculates the sum of P0 and P1 at infinity.
    """
    # Step 1: Define the given lambda parameters
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # Step 2: Set up the system of linear equations for the steady state (P_i'(t) = 0).
    # Let pi_i = P_i(+inf).
    # Equations used:
    # 1) -lambda_01*pi_0 + lambda_10*pi_1 = 0
    # 2)  lambda_12*pi_1 - (lambda_21 + lambda_23)*pi_2 = 0
    # 3)  lambda_23*pi_2 - lambda_31*pi_3 = 0
    # 4)  pi_0 + pi_1 + pi_2 + pi_3 = 1
    # We formulate this as a matrix equation Ax = b.

    # Step 3: Construct the coefficient matrix A and the constant vector b.
    A = np.array([
        [-lambda_01, lambda_10, 0, 0],
        [0, lambda_12, -(lambda_21 + lambda_23), 0],
        [0, 0, lambda_23, -lambda_31],
        [1, 1, 1, 1]
    ])

    b = np.array([0, 0, 0, 1])

    # Step 4: Solve the system Ax = b for x = [pi_0, pi_1, pi_2, pi_3].
    try:
        pi = np.linalg.solve(A, b)
        p0_inf = pi[0]
        p1_inf = pi[1]
    except np.linalg.LinAlgError:
        print("The system of equations is singular and cannot be solved.")
        return

    # Step 5: Calculate the required sum P0(+inf) + P1(+inf).
    result = p0_inf + p1_inf

    # Step 6: Print the individual probabilities and the final equation as requested.
    print("The steady-state probabilities are:")
    print(f"P0(+inf) = {p0_inf:.8f}")
    print(f"P1(+inf) = {p1_inf:.8f}")
    print(f"P2(+inf) = {pi[2]:.8f}")
    print(f"P3(+inf) = {pi[3]:.8f}")
    
    print("\nThe final equation is:")
    print(f"P0(+inf) + P1(+inf) = {p0_inf:.8f} + {p1_inf:.8f} = {result:.8f}")

if __name__ == '__main__':
    solve_steady_state_probabilities()
