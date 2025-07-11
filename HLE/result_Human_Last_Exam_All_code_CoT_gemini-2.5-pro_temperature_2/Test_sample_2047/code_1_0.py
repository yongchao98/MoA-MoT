import math
import numpy as np

def solve_l(d, lam):
    """
    Calculates the value of l(d, lambda) based on the derived formula.

    The formula is derived assuming the probability density p_{d, lambda}(x)
    is that of a Projected Normal distribution on the sphere S^{d-1},
    p(x) proportional to exp(-arccos(mu^T * x)^2 / (2 * lambda)),
    where mu is the mean direction (1_d / sqrt(d)).

    l(d, lambda) = (1 / (2 * lambda)) * [arccos(mu^T * x_2)^2 - arccos(mu^T * x_1)^2]
    
    Args:
        d (int): The dimension of the space, d >= 4.
        lam (float): The lambda parameter, lambda >= 1.

    Returns:
        float: The calculated value of l(d, lambda).
    """
    if not isinstance(d, int) or d < 4:
        raise ValueError("d must be an integer >= 4")
    if not isinstance(lam, (int, float)) or lam < 1:
        raise ValueError("lambda must be a number >= 1")

    # Step 1: Calculate mu^T * x_1 and mu^T * x_2
    # mu is a vector of 1/sqrt(d)
    # x1 is (e1+e2+e3)/sqrt(3)
    # mu.T * x1 = (1/sqrt(d)) * (1/sqrt(3)) * (1+1+1) = 3 / sqrt(3d) = sqrt(3/d)
    c1 = math.sqrt(3 / d)

    # x2 is (e3+e4)/sqrt(2)
    # mu.T * x2 = (1/sqrt(d)) * (1/sqrt(2)) * (1+1) = 2 / sqrt(2d) = sqrt(2/d)
    c2 = math.sqrt(2 / d)

    # Step 2: Calculate arccos of these values
    theta1 = math.acos(c1)
    theta2 = math.acos(c2)

    # Step 3: Square the angles
    theta1_sq = theta1**2
    theta2_sq = theta2**2

    # Step 4: Calculate the coefficient
    coeff = 1 / (2 * lam)

    # Step 5: Compute the final value
    result = coeff * (theta2_sq - theta1_sq)

    # As requested, printing the numbers in the final equation:
    # l(d, lambda) = (1 / (2 * lambda)) * [arccos(sqrt(2/d))^2 - arccos(sqrt(3/d))^2]
    print(f"For d = {d} and lambda = {lam}:")
    print(f"Equation: l(d, lambda) = (1 / (2 * {lam})) * [arccos(sqrt(2/{d}))^2 - arccos(sqrt(3/{d}))^2]")
    print(f"Let's calculate the terms:")
    print(f"  Term for x_1 (arccos(sqrt(3/{d}))^2):")
    print(f"    sqrt(3/{d}) = {c1}")
    print(f"    arccos({c1}) = {theta1}")
    print(f"    (arccos(...))^2 = {theta1_sq}")
    print(f"  Term for x_2 (arccos(sqrt(2/{d}))^2):")
    print(f"    sqrt(2/{d}) = {c2}")
    print(f"    arccos({c2}) = {theta2}")
    print(f"    (arccos(...))^2 = {theta2_sq}")
    print(f"  Coefficient (1 / (2 * {lam})): {coeff}")
    print(f"Final calculation: {coeff} * ({theta2_sq} - {theta1_sq})")
    print(f"Result: {result}")
    return result

# Example usage of the function.
# You can change these values to compute l(d, lambda) for a specific case.
d_example = 10
lambda_example = 2.0
final_answer = solve_l(d_example, lambda_example)
# <<<l(10, 2.0) is 0.060935359409426424>>>