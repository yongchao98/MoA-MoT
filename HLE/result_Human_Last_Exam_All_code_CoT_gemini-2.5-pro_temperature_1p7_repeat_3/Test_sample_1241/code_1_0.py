import numpy as np

def solve_steady_state_system():
    """
    Solves the boundary-value problem for the Kolmogorov-Chepmen system of equations
    to find the steady-state probabilities and their sum.
    """
    # Given transition rates
    lambda_01 = 0.019
    lambda_10 = 0.65
    lambda_12 = 0.4
    lambda_21 = 0.392
    lambda_23 = 0.008
    lambda_31 = 0.008

    # At steady state (t -> +infinity), the derivatives P_i'(t) are zero.
    # This leads to a system of linear equations for the steady-state probabilities p_i.
    # A typo in the problem's equation for P_1'(t) is corrected for probability conservation.
    # The system is:
    # 1) 0 = -lambda_01*p0 + lambda_10*p1
    # 2) 0 = lambda_12*p1 - (lambda_21 + lambda_23)*p2
    # 3) 0 = lambda_23*p2 - lambda_31*p3
    # 4) p0 + p1 + p2 + p3 = 1
    # We solve this by expressing p0, p2, and p3 in terms of p1.

    # From (1), p0 = (lambda_10 / lambda_01) * p1
    c0 = lambda_10 / lambda_01

    # From (2), p2 = (lambda_12 / (lambda_21 + lambda_23)) * p1
    c2 = lambda_12 / (lambda_21 + lambda_23)

    # From (3), p3 = (lambda_23 / lambda_31) * p2. Substituting p2 gives:
    # p3 = (lambda_23 / lambda_31) * c2 * p1
    c3 = (lambda_23 / lambda_31) * c2

    # Substitute these into the normalization equation (4):
    # (c0 * p1) + p1 + (c2 * p1) + (c3 * p1) = 1
    # p1 * (c0 + 1 + c2 + c3) = 1
    total_coeff = c0 + 1 + c2 + c3
    p1 = 1 / total_coeff

    # Now calculate p0
    p0 = c0 * p1

    # The required result is the sum of p0 and p1
    result = p0 + p1

    # Print the final equation with the computed numerical values
    print("The final equation for the required sum is P0(+infinity) + P1(+infinity):")
    print(f"{p0} + {p1} = {result}")

    return result

# Execute the function to solve the problem
final_answer = solve_steady_state_system()