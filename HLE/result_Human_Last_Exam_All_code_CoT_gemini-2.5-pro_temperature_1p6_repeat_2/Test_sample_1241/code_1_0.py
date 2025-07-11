import numpy as np

def solve_steady_state_probabilities():
    """
    Solves the steady-state probabilities for the given Kolmogorov-Chepmen system.
    """
    # Given lambda parameters
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # At steady state (t -> +infinity), the derivatives P_i'(t) are zero.
    # This leads to a system of linear equations. We use three independent
    # equations from the rate equations plus the normalization condition.
    # The system is Ax = b, where x = [P0, P1, P2, P3]^T.
    #
    # Equations used:
    # 1) -lambda_01*P0 + lambda_10*P1 = 0
    # 2) lambda_12*P1 - (lambda_21 + lambda_23)*P2 = 0
    # 3) lambda_23*P2 - lambda_31*P3 = 0
    # 4) P0 + P1 + P2 + P3 = 1

    # Define the coefficient matrix A
    A = np.array([
        [-lambda_01, lambda_10, 0, 0],
        [0, lambda_12, -(lambda_21 + lambda_23), 0],
        [0, 0, lambda_23, -lambda_31],
        [1, 1, 1, 1]
    ])

    # Define the constant vector b
    b = np.array([0, 0, 0, 1])

    try:
        # Solve the system for P = [P0, P1, P2, P3]
        P = np.linalg.solve(A, b)
        P0, P1 = P[0], P[1]

        # Calculate the required sum P0 + P1
        result = P0 + P1

        # Print the final result in the requested format
        print("The question asks for the value of P0(+infinity) + P1(+infinity).")
        print("This is calculated by solving the steady-state equations.")
        print("The final equation is:")
        print(f"{P0} + {P1} = {result}")

    except np.linalg.LinAlgError:
        print("Error: The system of equations is singular and cannot be solved.")
        print("This might indicate that the selected equations are not linearly independent.")

solve_steady_state_probabilities()