import math

def solve_cone_and_spheres_problem():
    """
    Solves the problem of fitting small spheres in a cone around a larger inscribed sphere.

    The core relationship derived from the geometry of the problem is:
    (sqrt(H**2 + R**2) - R) / H = sin(pi / n)
    
    where:
    H = height of the cone (integer)
    R = base radius of the cone (integer)
    n = number of small spheres (integer, n >= 3)
    
    For H and R to be integers, their ratio R/H must be a rational number.
    Rearranging the equation gives:
    R/H = (1 - sin(pi/n)**2) / (2 * sin(pi/n))
    
    We need to find an integer n for which R/H is rational.
    Analysis shows this only occurs when sin(pi/n) is rational, which by Niven's
    theorem means sin(pi/n) must be 1/2 (for n>=3). This gives n=6.
    """
    
    # Let's verify for n=6
    n = 6
    # sin(pi/6) = 0.5
    c = math.sin(math.pi / n)
    
    # Calculate the required ratio R/H
    # We expect this to be 3/4 = 0.75
    ratio_R_div_H = (1 - c**2) / (2 * c)
    
    # We found the solution n=6, which gives R/H = 3/4.
    # The simplest positive integers for H and R are H=4 and R=3.
    H = 4
    R = 3
    
    print("Yes, it is possible.")
    print(f"The number of smaller spheres that can fit is n = {n}.")
    print("-" * 30)
    print("This is possible for a cone with integer height and radius, for example H = 4 and R = 3.")
    print("Let's verify the solution by plugging these values into the equation:")
    print("(sqrt(H^2 + R^2) - R) / H = sin(pi / n)\n")
    
    # Left Hand Side (LHS) calculation
    lhs_sqrt_val = H**2 + R**2
    lhs_sqrt_res = math.sqrt(lhs_sqrt_val)
    lhs_numerator = lhs_sqrt_res - R
    lhs_final = lhs_numerator / H
    
    # Right Hand Side (RHS) calculation
    rhs_final = math.sin(math.pi / n)

    # Print the verification step-by-step
    print(f"(sqrt({H}^2 + {R}^2) - {R}) / {H} = sin(pi / {n})")
    print(f"(sqrt({H**2} + {R**2}) - {R}) / {H} = {rhs_final}")
    print(f"(sqrt({lhs_sqrt_val}) - {R}) / {H} = {rhs_final}")
    print(f"({lhs_sqrt_res} - {R}) / {H} = {rhs_final}")
    print(f"({lhs_numerator}) / {H} = {rhs_final}")
    print(f"{lhs_final} = {rhs_final}")
    print("\nThe equality holds, so the solution is correct.")

solve_cone_and_spheres_problem()
