import math

def solve_star_angles():
    """
    Calculates the value of (1 - cos(theta_14)) / (1 - cos(theta_34))
    in the second frame of reference.
    """

    # As derived in the plan, the final ratio Q can be expressed in terms of
    # the known angles in the second frame of reference.
    # Q = (1 - cos(theta'_12)) / (1 - cos(theta'_13))

    # Given angles in the second frame of reference (in radians)
    theta_12_prime = math.pi / 2
    theta_13_prime = 3 * math.pi / 4

    # Calculate the cosine of these angles
    # Note: math.cos(math.pi / 2) might be a very small number close to 0 due to
    # floating point precision, which is physically correct.
    cos_theta_12_prime = math.cos(theta_12_prime)
    cos_theta_13_prime = math.cos(theta_13_prime)

    # The final equation is:
    # Value = (1 - cos(pi/2)) / (1 - cos(3*pi/4))

    # Numerator of the final equation
    numerator = 1 - cos_theta_12_prime

    # Denominator of the final equation
    denominator = 1 - cos_theta_13_prime

    # Calculate the final value
    value = numerator / denominator

    print("This program solves for the ratio (1 - cos(theta'_14)) / (1 - cos(theta'_34)).")
    print("Based on the principles of special relativity, this ratio simplifies to:")
    print("  Ratio = (1 - cos(theta'_12)) / (1 - cos(theta'_13))\n")

    print("Here are the numbers in the final equation:")
    print(f"theta'_12 = pi/2 radians")
    print(f"cos(theta'_12) = {cos_theta_12_prime:.4f}")
    print(f"theta'_13 = 3*pi/4 radians")
    print(f"cos(theta'_13) = {cos_theta_13_prime:.4f}\n")

    print("Calculating the terms:")
    print(f"Numerator = 1 - {cos_theta_12_prime:.4f} = {numerator:.4f}")
    print(f"Denominator = 1 - ({cos_theta_13_prime:.4f}) = {denominator:.4f}\n")

    print("Final result:")
    print(f"Ratio = {numerator:.4f} / {denominator:.4f} = {value:.4f}")
    # The exact value is 2 - sqrt(2)
    exact_value = 2 - math.sqrt(2)
    print(f"The exact value is 2 - sqrt(2), which is approximately {exact_value:.4f}.")

solve_star_angles()