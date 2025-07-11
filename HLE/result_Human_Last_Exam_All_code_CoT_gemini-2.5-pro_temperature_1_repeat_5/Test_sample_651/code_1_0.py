import sympy as sp

def solve_limit_billiards():
    """
    This function calculates the limit of M(theta) as theta approaches 0.
    M(theta) is the supremum of the angle alpha between a billiard trajectory
    and the inner normal to a side of a triangle with vertices 0, 5, and 5*e^(i*theta).

    The calculation proceeds with the following steps, using approximations for small theta:
    1.  Determine the angle of the inner normal vector to side A.
    2.  Determine the angles of the extremal trajectory vectors.
    3.  Calculate the maximum angle M(theta) from these values.
    4.  Find the limit of M(theta) as theta approaches 0.
    """
    
    # For the final output, we demonstrate the logic with symbolic representation
    # and print the numerical coefficients as the logic unfolds.
    
    print("Step-by-step derivation of the limit:")
    
    # 1. Angle of the Inner Normal Vector (n)
    # The inner normal vector 'n' to the side A (from 5 to 5*e^(i*theta)) is
    # perpendicular to the side's direction and points inward. Its direction vector
    # is proportional to (sin(theta), 1-cos(theta)).
    # For small theta, this is approx. (theta, theta^2/2).
    # The angle of the inner normal vector is arg(n) = atan2(theta^2/2, theta) which is approx. theta/2.
    print("1. The angle of the inner normal vector, arg(n), is approximately theta/2.")
    
    # 2. Angles of Extremal Trajectory Vectors (v)
    # The extremal trajectories are from an endpoint of the starting arc {x} to an
    # endpoint of the side A {y}.
    #
    # First extremal vector, v1 (from x=1 to y=5*e^(i*theta)):
    # v1 = 5*e^(i*theta) - 1 approx. 5*(1 + i*theta) - 1 = 4 + 5*i*theta
    # arg(v1) approx. atan(5*theta / 4) approx. 5*theta / 4.
    print("2. The first extremal trajectory has an angle arg(v1) ≈ (5/4)*theta.")

    # Second extremal vector, v2 (from x=e^(i*theta) to y=5):
    # v2 = 5 - e^(i*theta) approx. 5 - (1 + i*theta) = 4 - i*theta
    # arg(v2) approx. atan(-theta / 4) approx. -theta / 4.
    print("3. The second extremal trajectory has an angle arg(v2) ≈ (-1/4)*theta.")

    # 3. Calculation of the Supremum Angle M(theta)
    # The angle alpha is the absolute difference between the trajectory angle and the normal angle.
    # M(theta) is the supremum (maximum) of these possible angles.
    print("4. The supremum angle M(theta) is the maximum of the angles from the extremal trajectories.")
    
    # Angle for v1:
    # alpha1 approx. |arg(v1) - arg(n)|
    print("   - For v1, the angle alpha1 ≈ |(5/4)*theta - (1/2)*theta| = |(5/4)*theta - (2/4)*theta| = |(3/4)*theta| = (3/4)*theta.")
    
    # Angle for v2:
    # alpha2 approx. |arg(v2) - arg(n)|
    print("   - For v2, the angle alpha2 ≈ |(-1/4)*theta - (1/2)*theta| = |(-1/4)*theta - (2/4)*theta| = |(-3/4)*theta| = (3/4)*theta.")
    
    print("5. Thus, for small theta, M(theta) is approximately (3/4)*theta.")
    
    # 4. Final Limit Calculation
    # We want to find the limit of M(theta) as theta -> 0.
    theta_sym = sp.Symbol('theta')
    # Using the approximation M(theta) ≈ (3/4)*theta
    M_theta_approx = (sp.S(3)/4) * theta_sym
    limit_M = sp.limit(M_theta_approx, theta_sym, 0)
    
    print(f"6. The limit of M(theta) as theta -> 0 is the limit of (3/4)*theta, which is {limit_M}.")
    
    print("\nFinal Result:")
    print(limit_M)

# Execute the function to print the solution
solve_limit_billiards()