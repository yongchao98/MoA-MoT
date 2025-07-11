import numpy as np

def solve_steady_state():
    """
    Solves the boundary-value problem for the Kolmogorov-Chapman system
    to find the steady-state probabilities and their sum.
    """
    # The problem describes a continuous-time Markov chain with 4 states (0, 1, 2, 3).
    # We need to find the sum of the steady-state probabilities for states 0 and 1,
    # denoted as P0(+infinity) and P1(+infinity).

    # At steady state (as t -> +infinity), the probabilities are constant,
    # so their time derivatives are zero: P_i'(t) = 0.
    # Let pi_i = P_i(+infinity). The system of differential equations becomes a
    # system of linear algebraic equations. We also know that the sum of
    # probabilities must be 1.
    #
    # The system of equations is:
    # 1) -lambda_01*pi_0 + lambda_10*pi_1 = 0
    # 2) lambda_12*pi_1 - (lambda_21 + lambda_23)*pi_2 = 0
    # 3) lambda_23*pi_2 - lambda_31*pi_3 = 0
    # 4) pi_0 + pi_1 + pi_2 + pi_3 = 1
    #
    # We can solve this system using linear algebra (A * pi = b).

    # Define the given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # Construct the coefficient matrix A for the system of linear equations.
    A = np.array([
        [-lambda_01, lambda_10, 0, 0],
        [0, lambda_12, -(lambda_21 + lambda_23), 0],
        [0, 0, lambda_23, -lambda_31],
        [1, 1, 1, 1]
    ])

    # Construct the constant vector b.
    b = np.array([0, 0, 0, 1])

    # Solve the system A * pi = b for the steady-state probabilities pi.
    try:
        pi = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        print("The system of equations is singular and cannot be solved.")
        return

    # Extract the probabilities for states 0 and 1.
    pi_0 = pi[0]
    pi_1 = pi[1]

    # Calculate the required sum.
    result = pi_0 + pi_1

    # Print the final equation with the calculated values, showing each number.
    print(f"P0(+infinity) + P1(+infinity) = {pi_0} + {pi_1} = {result}")

if __name__ == "__main__":
    solve_steady_state()