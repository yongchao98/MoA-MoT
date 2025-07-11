import math

def solve_star_angle_problem():
    """
    Calculates the ratio (1 - cos(theta'_14)) / (1 - cos(theta'_34)) based on principles of special relativity.
    """

    # In Frame 1, the stars form a regular tetrahedron. The cosine of the angle
    # between any two stars is -1/3.
    cos_theta_frame1 = -1.0 / 3.0
    # The term (1 - cos(theta)) is constant for any pair in Frame 1.
    term_frame1 = 1.0 - cos_theta_frame1

    # In Frame 2, we are given the angles for pairs (1,2), (1,3), and (2,3).
    # From these, we can find relationships between the Doppler factors D_i.
    # The key relation is (1 - cos(theta')) = (1 - cos(theta)) / (D_i * D_j).

    # For S1, S3: theta'_13 = 3pi/4, so cos(theta'_13) = -1/sqrt(2)
    cos_theta_prime_13 = -1.0 / math.sqrt(2)
    one_minus_cos_theta_prime_13 = 1.0 - cos_theta_prime_13
    
    # The target quantity is (1 - cos(theta'_14)) / (1 - cos(theta'_34)).
    # Using the relativistic formula, this ratio simplifies to D3 / D1.

    # We can establish two equations for the Doppler factors:
    # Eq A: D1 * D2 = term_frame1 (from theta'_12 = pi/2)
    # Eq B: D1 * D3 = term_frame1 / one_minus_cos_theta_prime_13
    # We also know D1 = D2, so Eq A becomes D1^2 = term_frame1.
    
    # To find D3/D1, we divide Eq B by D1^2:
    # (D1 * D3) / D1^2 = (term_frame1 / one_minus_cos_theta_prime_13) / term_frame1
    # D3 / D1 = 1 / one_minus_cos_theta_prime_13
    
    result_ratio = 1.0 / one_minus_cos_theta_prime_13

    # The result is 1 / (1 + 1/sqrt(2)), which simplifies to 2 - sqrt(2).
    # Here are the numbers in the final equation.
    num1 = 2
    num2 = 2 # from sqrt(2)

    print(f"The simplified expression for the ratio is: {num1} - sqrt({num2})")
    print(f"The numbers in this final equation are {num1} and {num2}.")
    
    # Calculate and print the numerical value.
    numerical_value = num1 - math.sqrt(num2)
    print(f"The numerical value of the ratio is: {numerical_value}")

solve_star_angle_problem()