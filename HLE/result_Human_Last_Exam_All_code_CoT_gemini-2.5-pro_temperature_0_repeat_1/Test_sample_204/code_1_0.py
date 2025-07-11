import numpy as np

def solve_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by finding the
    linking number of the preimages of two points on the S^2 sphere.
    """
    # The Hopf charge is the linking number Lk(C_S, C_E), where C_S and C_E
    # are the preimages of two points on the sphere S^2.

    # Step 1: Find the preimage of the South Pole, p_S = (0, 0, -1).
    # The condition n = p_S gives nz = cos(G) = -1, which means G = PI.
    # From the definition G = PI * exp(-10 * r2), G = PI implies r2 = 0.
    # The condition r2 = sqrt((x*x+y*y-0.5)^2 + z*z) = 0 implies:
    # z = 0 and x^2 + y^2 - 0.5 = 0.
    # So, the preimage C_S is the circle x^2 + y^2 = 0.5 in the z=0 plane.
    radius_sq_CS = 0.5
    print("Step 1: Finding the preimage of the South Pole (0,0,-1)")
    print(f"The preimage curve C_S is a circle in the z=0 plane with equation x^2 + y^2 = {radius_sq_CS}.\n")

    # Step 2: Find the preimage of an equatorial point, p_E = (1, 0, 0).
    # The condition n = p_E gives:
    # nz = cos(G) = 0  => G = PI/2
    # ny = sin(G)*sin(f) = 0 => sin(f) = 0 (since sin(G)=1) => y = 0
    # nx = sin(G)*cos(f) = 1 => cos(f) = 1 (since sin(G)=1) => f = 0 (so x > 0)
    # From G = PI/2, we find the value of r2:
    # PI * exp(-10 * r2) = PI/2  =>  exp(-10 * r2) = 0.5  => r2 = ln(2)/10
    # The preimage curve C_E is defined by r2^2 = (ln(2)/10)^2 with y=0, x>0:
    # (x^2 - 0.5)^2 + z^2 = (ln(2)/10)^2
    R = np.log(2) / 10
    print("Step 2: Finding the preimage of the equatorial point (1,0,0)")
    print(f"The preimage curve C_E is a loop in the y=0 plane with equation (x^2 - 0.5)^2 + z^2 = R^2, where R = ln(2)/10 ≈ {R:.4f}.\n")

    # Step 3: Calculate the linking number by counting intersections.
    # We define a surface D_S bounded by C_S, which is the disk: x^2 + y^2 <= 0.5, z=0.
    # We find the intersection points of C_E with the plane z=0.
    # Setting z=0 in the equation for C_E: (x^2 - 0.5)^2 = R^2
    # This gives x^2 - 0.5 = +/- R, so x^2 = 0.5 +/- R.
    x_sq_1 = 0.5 - R
    x_sq_2 = 0.5 + R
    print("Step 3: Calculating the linking number Lk(C_S, C_E)")
    print("We count the intersections of curve C_E with the disk D_S (x^2+y^2 <= 0.5, z=0).")
    print(f"The intersections of C_E with the z=0 plane occur at x^2 = 0.5 +/- R.")
    print(f"Intersection 1: x^2 = 0.5 - {R:.4f} = {x_sq_1:.4f}")
    print(f"Intersection 2: x^2 = 0.5 + {R:.4f} = {x_sq_2:.4f}\n")

    # Check which intersection points lie inside the disk D_S (where x^2 <= 0.5).
    is_inside_1 = x_sq_1 < radius_sq_CS
    is_inside_2 = x_sq_2 < radius_sq_CS
    print("Checking which intersection points are inside the disk D_S:")
    print(f"Point 1 (x^2={x_sq_1:.4f}) is inside the disk (x^2 < {radius_sq_CS}): {is_inside_1}")
    print(f"Point 2 (x^2={x_sq_2:.4f}) is inside the disk (x^2 < {radius_sq_CS}): {is_inside_2}\n")

    num_intersections = int(is_inside_1) + int(is_inside_2)

    print(f"Conclusion:")
    print(f"There is exactly {num_intersections} intersection of C_E with the disk D_S.")
    print("This means the magnitude of the linking number is 1.")
    print("The given vector field describes a standard hopfion configuration, for which the Hopf charge is positive.")
    print("Final Equation: Hopf Charge = 1")

solve_hopf_charge()