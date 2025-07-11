import math

def solve_limit_inner_product():
    """
    Calculates the limiting inner product based on the geometric relationship
    between three coplanar vectors.
    """
    # Given values from the problem statement
    c1 = 0.9375  # lim <h_p, b_p>
    c2 = 0.9     # lim <h_p, z_p>

    # The problem asks to find theta = lim <b_p, z_p>.
    # Based on the geometric insight that the three vectors are coplanar,
    # the relationship between their inner products (cosines of the angles) is:
    # theta = c1*c2 + sqrt((1-c1^2)*(1-c2^2))
    # We choose the '+' sign based on the physical argument that the leading eigenvector h_p
    # lies geometrically between the signal vector b_p and the dominant noise vector z_p.

    # Calculate the terms
    c1_squared = c1**2
    c2_squared = c2**2
    term1 = c1 * c2
    term2_radicand = (1 - c1_squared) * (1 - c2_squared)
    term2 = math.sqrt(term2_radicand)
    
    # Calculate the final value of theta
    theta = term1 + term2
    
    # Output the formula with the numbers plugged in
    print("The formula derived from the geometric relationship is:")
    print("theta = c1*c2 + sqrt((1 - c1^2) * (1 - c2^2))")
    print("\nPlugging in the given values:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"theta = {c1} * {c2} + sqrt((1 - {c1}^2) * (1 - {c2}^2))")
    print(f"theta = {term1} + sqrt(({1-c1_squared}) * ({1-c2_squared}))")
    print(f"theta = {term1} + sqrt({term2_radicand})")
    print(f"theta = {term1} + {term2}")
    print(f"\nThe final calculated value is:")
    print(f"theta = {theta}")

solve_limit_inner_product()
<<<0.9954345332763321>>>