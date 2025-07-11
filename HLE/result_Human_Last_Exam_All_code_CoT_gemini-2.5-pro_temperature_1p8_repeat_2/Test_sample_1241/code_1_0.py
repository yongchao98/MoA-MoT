import numpy as np

def solve_steady_state():
    """
    Solves the steady-state problem for the given Kolmogorov-Chepmen system.
    """
    # Given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # At steady state (t -> +inf), all P_i'(t) = 0.
    # Let pi_i = P_i(+inf). The system becomes a set of linear equations:
    # 1) 0 = -lambda_01 * pi_0 + lambda_10 * pi_1
    # 2) 0 = lambda_12 * pi_1 - (lambda_21 + lambda_23) * pi_2
    # 3) 0 = lambda_23 * pi_2 - lambda_31 * pi_3
    # 4) 1 = pi_0 + pi_1 + pi_2 + pi_3
    # Note: The second original equation is redundant if the others hold true,
    # as the sum of derivatives must be zero.

    # From (1), we can express pi_0 in terms of pi_1:
    # pi_0 = (lambda_10 / lambda_01) * pi_1
    c0 = lambda_10 / lambda_01

    # From (2), we can express pi_2 in terms of pi_1:
    # pi_2 = (lambda_12 / (lambda_21 + lambda_23)) * pi_1
    c2 = lambda_12 / (lambda_21 + lambda_23)

    # From (3), we can express pi_3 in terms of pi_2:
    # pi_3 = (lambda_23 / lambda_31) * pi_2
    # Substituting the expression for pi_2:
    # pi_3 = (lambda_23 / lambda_31) * c2 * pi_1
    c3 = (lambda_23 / lambda_31) * c2

    # Now substitute these into the normalization equation (4):
    # c0*pi_1 + pi_1 + c2*pi_1 + c3*pi_1 = 1
    # pi_1 * (c0 + 1 + c2 + c3) = 1
    # So, pi_1 = 1 / (c0 + 1 + c2 + c3)
    denominator = c0 + 1 + c2 + c3
    pi_1 = 1 / denominator

    # Now we can calculate pi_0:
    pi_0 = c0 * pi_1

    # The problem asks for P0(+inf) + P1(+inf) = pi_0 + pi_1
    result = pi_0 + pi_1
    
    # Print the final result as a numerical equation
    print("The goal is to compute P0(+infinity) + P1(+infinity).")
    print(f"The steady-state probability for state 0 is P0(+infinity) = {pi_0}")
    print(f"The steady-state probability for state 1 is P1(+infinity) = {pi_1}")
    print("\nThe final equation is:")
    print(f"{pi_0} + {pi_1} = {result}")

solve_steady_state()