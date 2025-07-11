import math

def calculate_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by analyzing
    the linking number of the preimages of two points on the sphere S^2.
    """
    print("The Hopf charge H is the linking number Lk(C1, C2) of the preimages")
    print("of two distinct points on the target sphere S^2.\n")

    # We choose two points on S^2:
    # p1 = (0, 0, -1) (South Pole)
    # p2 = (1, 0, 0)  (A point on the equator)

    # --- Step 1: Find the preimage C1 of the South Pole p1 ---
    # The condition n_z = cos(G) = -1 implies G = PI.
    # The field G is defined as G = PI * exp(-10 * r2).
    # Setting G = PI gives exp(-10 * r2) = 1, which means r2 = 0.
    # The term r2 = sqrt((x*x+y*y-0.5)^2 + z*z).
    # For r2 to be 0, both terms under the square root must be zero:
    # z = 0 and x*x+y*y-0.5 = 0.
    # This defines the curve C1.
    c1_radius_sq = 0.5
    c1_z = 0.0
    print("Curve C1 (preimage of the South Pole):")
    print(f"  Equation: x^2 + y^2 = {c1_radius_sq}, z = {c1_z}")
    print(f"  This is a circle in the xy-plane with radius {math.sqrt(c1_radius_sq):.4f}.\n")

    # --- Step 2: Find the preimage C2 of the equatorial point p2 ---
    # The condition for p2 = (1,0,0) is:
    # n_z = cos(G) = 0  => G = PI/2
    # n_x = sin(G)*cos(f) = 1*cos(f) = 1 => f = 0
    # n_y = sin(G)*sin(f) = 1*sin(f) = 0 => f = 0
    # The condition f = atan2(y,x) = 0 implies y=0 and x>0.
    # The condition G = PI/2 implies PI * exp(-10 * r2) = PI/2.
    # This simplifies to exp(-10 * r2) = 0.5, or r2 = ln(2) / 10.
    C = math.log(2) / 10.0
    # The equation for C2 is r2^2 = C^2, with y=0 and x>0:
    # (x^2 - 0.5)^2 + z^2 = C^2
    print("Curve C2 (preimage of the point (1,0,0)):")
    print(f"  Equation: (x^2 - 0.5)^2 + z^2 = {C**2:.6f}, with y=0 and x>0")
    print("  This is a closed loop in the xz-plane.\n")

    # --- Step 3: Determine the linking number of C1 and C2 ---
    # C1 is a circle in the xy-plane. C2 is a loop in the xz-plane.
    # Their linking number is determined by how C2 passes through the disk bounded by C1.
    # Let's find where C2 intersects the xy-plane (where z=0).
    # (x^2 - 0.5)^2 = C^2 => x^2 - 0.5 = +/- C => x^2 = 0.5 +/- C
    x_sq_inside = 0.5 - C
    x_sq_outside = 0.5 + C
    x_inside = math.sqrt(x_sq_inside)
    x_outside = math.sqrt(x_sq_outside)
    r1 = math.sqrt(c1_radius_sq)

    print("Analysis of the linking between C1 and C2:")
    print(f"The loop C2 crosses the xy-plane at two points: x={x_inside:.4f} and x={x_outside:.4f}.")
    print(f"The circle C1 has a radius of {r1:.4f}.")
    print(f"One crossing point ({x_inside:.4f}) is INSIDE C1 (since {x_inside:.4f} < {r1:.4f}).")
    print(f"The other crossing point ({x_outside:.4f}) is OUTSIDE C1 (since {x_outside:.4f} > {r1:.4f}).")
    print("Since the loop C2 passes through the disk bounded by the circle C1 exactly once,")
    print("their linking number is 1.\n")

    # The Hopf charge is also given by the product of two integer winding numbers, k and N.
    # k is the winding in the azimuthal angle f, which is 1.
    # N is the winding of the profile G, which maps the (rho,z) plane to the sphere poles, which is 1.
    k = 1
    N = 1
    H = k * N
    
    print("The Hopf charge H can be understood as the product of two winding numbers.")
    print(f"Azimuthal winding number k = {k}")
    print(f"Polar winding number N = {N}")
    print(f"Final Result: H = {k} * {N} = {H}")

calculate_hopf_charge()
<<<1>>>