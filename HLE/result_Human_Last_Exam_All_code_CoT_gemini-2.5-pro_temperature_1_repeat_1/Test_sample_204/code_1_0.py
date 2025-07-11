import math

def calculate_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by determining the
    linking number of the preimages of two points on the target sphere.
    """
    PI = math.pi

    # Step 1: Define the mapping functions.
    # The vector field n is defined by a map to spherical coordinates (G, f)
    # on the S^2 sphere.
    # f = atan2(y,x)
    # r2 = sqrt((x*x+y*y-0.5)**2 + z*z)
    # G = PI * exp(-10*r2)

    print("Step 1: Analyzing the vector field definition.")
    print("The vector field (nx,ny,nz) is a map from R^3 to S^2.")
    print("The mapping is defined by:")
    print("f = atan2(y,x)")
    print("r2 = sqrt((x*x+y*y-0.5)**2 + z*z)")
    print("G = PI * exp(-10*r2)")
    print("nx = sin(G)*cos(f), ny = sin(G)*sin(f), nz = cos(G)\n")

    # Step 2: Choose two points on the target sphere S^2.
    # We choose the South Pole (v_S) and an equatorial point (v_E).
    v_S = (0, 0, -1)
    v_E = (1, 0, 0)
    print("Step 2: Choose two regular points on the target sphere S^2.")
    print(f"Point 1 (South Pole): v_S = {v_S}")
    print(f"Point 2 (Equatorial): v_E = {v_E}\n")

    # Step 3: Find the preimage of the South Pole (C_S).
    # For n = v_S, we need nz = cos(G) = -1. This implies G = PI.
    # PI = PI * exp(-10*r2)  =>  exp(-10*r2) = 1  =>  -10*r2 = 0  =>  r2 = 0.
    # r2 = sqrt((x*x+y*y-0.5)**2 + z**2) = 0.
    # This holds if and only if z=0 and x*x+y*y-0.5 = 0.
    # So, the preimage C_S is a circle in the xy-plane.
    cs_radius_sq = 0.5
    cs_radius = math.sqrt(cs_radius_sq)
    print("Step 3: Calculate the preimage of the South Pole (C_S).")
    print("n = v_S => cos(G) = -1 => G = PI.")
    print("This implies r2 = 0, which gives z = 0 and x*x + y*y = 0.5.")
    print(f"The preimage C_S is a circle in the z=0 plane with radius sqrt({cs_radius_sq}) = {cs_radius:.4f}.\n")

    # Step 4: Find the preimage of the equatorial point (C_E).
    # For n = v_E, we need nx=1, ny=0, nz=0.
    # nz = cos(G) = 0  =>  G = PI/2.
    # ny = sin(G)*sin(f) = 1 * sin(f) = 0  => f = 0 or PI.
    # nx = sin(G)*cos(f) = 1 * cos(f) = 1  => f = 0.
    # f = atan2(y,x) = 0 implies y=0 and x>0.
    # From G = PI/2, we have PI/2 = PI * exp(-10*r2) => exp(-10*r2) = 0.5.
    # -10*r2 = ln(0.5) = -ln(2) => r2 = ln(2)/10.
    # The preimage C_E is a curve in the y=0, x>0 half-plane (the xz-plane).
    r2_ce = math.log(2) / 10
    r2_ce_sq = r2_ce**2
    print("Step 4: Calculate the preimage of the equatorial point (C_E).")
    print("n = v_E => G = PI/2 and f = 0.")
    print("f=0 implies y=0 and x>0.")
    print(f"G=PI/2 implies r2 = ln(2)/10 = {r2_ce:.4f}.")
    print(f"The preimage C_E is a loop in the xz-plane defined by the equation:")
    print(f"(x*x - {cs_radius_sq})**2 + z*z = {r2_ce_sq:.6f}\n")

    # Step 5: Calculate the linking number Lk(C_S, C_E).
    # We find the intersection points of the loop C_E with the disk D_S
    # bounded by C_S (i.e., x*x+y*y <= 0.5, z=0).
    # C_E intersects the z=0 plane when (x*x - 0.5)**2 = r2_ce_sq.
    # x*x - 0.5 = +/- r2_ce
    # x*x = 0.5 +/- r2_ce
    x_sq_1 = cs_radius_sq - r2_ce
    x_sq_2 = cs_radius_sq + r2_ce
    print("Step 5: Calculate the linking number Lk(C_S, C_E).")
    print("We count how many times the loop C_E pierces the disk bounded by C_S.")
    print("C_E intersects the z=0 plane at two points:")
    
    # Check if intersection points are inside the disk C_S.
    # The disk is x*x+y*y <= 0.5. Since y=0 for C_E, this is x*x <= 0.5.
    
    print(f"Point 1: x^2 = {cs_radius_sq} - {r2_ce:.4f} = {x_sq_1:.4f}")
    if x_sq_1 < cs_radius_sq:
        print(f"  Since {x_sq_1:.4f} < {cs_radius_sq}, this point is INSIDE the disk C_S.")
        num_intersections = 1
    else:
        print(f"  Since {x_sq_1:.4f} > {cs_radius_sq}, this point is OUTSIDE the disk C_S.")
        num_intersections = 0

    print(f"Point 2: x^2 = {cs_radius_sq} + {r2_ce:.4f} = {x_sq_2:.4f}")
    if x_sq_2 < cs_radius_sq:
        print(f"  Since {x_sq_2:.4f} < {cs_radius_sq}, this point is INSIDE the disk C_S.")
        num_intersections += 1
    else:
        print(f"  Since {x_sq_2:.4f} > {cs_radius_sq}, this point is OUTSIDE the disk C_S.")
        
    print(f"\nThe loop C_E pierces the disk D_S exactly {num_intersections} time(s).\n")
    
    # Step 6: Conclude the Hopf charge.
    # The magnitude of the linking number is the number of piercings.
    # By convention for this field configuration, the sign is positive.
    hopf_charge = num_intersections
    print("Step 6: Conclude the Hopf Charge.")
    print("The Hopf charge is an integer equal to the linking number Lk(C_S, C_E).")
    print(f"The magnitude of the linking number is {hopf_charge}.")
    print("By convention, the sign for this standard field configuration is positive.")
    print("\nFinal equation:")
    print(f"Hopf Charge = {hopf_charge}")
    return hopf_charge

# Execute the calculation and print the final answer in the required format.
if __name__ == "__main__":
    final_answer = calculate_hopf_charge()
    print(f"\n<<<{final_answer}>>>")
