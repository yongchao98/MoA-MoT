import math

def solve_limit_inner_product():
    """
    This function solves for the limiting inner product based on the given values.
    """
    # Let A and B be the given limits of the inner products.
    # A = lim <h_p, z_p>
    A = 0.9
    # B = lim <h_p, b_p>
    B = 0.9375

    # Based on the coplanarity of the vectors h_p, b_p, and z_p in the limit,
    # the unknown value x = lim <b_p, z_p> satisfies the quadratic equation:
    # x^2 - 2*A*B*x + (A^2 + B^2 - 1) = 0.

    # Define the coefficients of the quadratic equation a*x^2 + b*x + c = 0
    a_coeff = 1
    b_coeff = -2 * A * B
    c_coeff = A**2 + B**2 - 1

    print("The problem is solved by assuming the vectors h_p, b_p, and z_p are coplanar in the limit.")
    print("This leads to a quadratic equation for x = lim <b_p, z_p> of the form: a*x^2 + b*x + c = 0.")
    print("The equation with the given numbers is:")
    print(f"({a_coeff}) * x^2 + ({b_coeff}) * x + ({c_coeff}) = 0")
    print("")

    # Calculate the discriminant of the quadratic equation
    discriminant = b_coeff**2 - 4 * a_coeff * c_coeff

    if discriminant < 0:
        print("The equation has no real solutions, which indicates an issue with the premises.")
        return

    # Calculate the two possible roots for x
    sqrt_discriminant = math.sqrt(discriminant)
    x1 = (-b_coeff + sqrt_discriminant) / (2 * a_coeff)
    x2 = (-b_coeff - sqrt_discriminant) / (2 * a_coeff)

    # The two roots correspond to different geometric arrangements.
    # Let theta_A, theta_B, theta_x be the angles whose cosines are A, B, and x.
    # The larger root corresponds to x = cos(theta_A - theta_B) or cos(theta_B - theta_A).
    # The smaller root corresponds to x = cos(theta_A + theta_B).
    # The physical context suggests that the perturbed eigenvector h_p lies between the
    # original vector b_p and the noise direction z_p. This implies the angles add up,
    # i.e., theta_x = theta_A + theta_B. Therefore, we must choose the smaller root.
    
    final_answer = min(x1, x2)

    print(f"The two solutions for x are {x1} and {x2}.")
    print("Based on the physical reasoning, the smaller solution is the correct one.")
    print("\nThe final answer is:")
    print(final_answer)

solve_limit_inner_product()