import numpy as np

def solve_hopf_charge():
    """
    Calculates the Hopf charge of the given vector field by finding the
    linking number of the preimages of two points on the S^2 sphere.
    """
    # 1. Introduction
    print("Calculating the Hopf charge using the Whitehead formula.")
    print("This is interpreted as the linking number of the preimages of two points on the S^2 sphere.")
    print("The vector field n(x,y,z) is a map from R^3 to S^2.")
    print("Hopf Charge H = Lk(n^-1(p1), n^-1(p2)) for two points p1, p2 on S^2.")
    print("-" * 50)

    # 2. Define the field parameters from the problem description
    # f = atan2(y,x)
    # r2 = sqrt((x*x+y*y-0.5)*(x*x+y*y-0.5)+z*z)
    # G = PI*(exp(-10*r2))
    # n = (sin(G)*cos(f), sin(G)*sin(f), cos(G))
    print("Field Definition:")
    print("n = (sin(G)*cos(f), sin(G)*sin(f), cos(G))")
    print("where G = pi * exp(-10 * sqrt((rho^2 - 0.5)^2 + z^2)) and f = atan2(y,x).")
    print("-" * 50)

    # 3. Choose points p1 and p2
    print("We choose two points on the sphere S^2:")
    p1_name = "South Pole"
    p1_coords = "(0, 0, -1)"
    p2_name = "an Equatorial Point"
    p2_coords = "(1, 0, 0)"
    print(f"p1 = {p1_name} {p1_coords}")
    print(f"p2 = {p2_name} {p2_coords}")
    print("-" * 50)

    # 4. Find the preimage of p1 (South Pole)
    print(f"Step 1: Find the preimage of p1 = {p1_coords}")
    print("n = (0, 0, -1) requires cos(G) = -1.")
    print("This implies G = pi * (2k+1) for some integer k.")
    print("The field is defined with G = pi * exp(-10*r2). The range of G is (0, pi].")
    print("So, we must have G = pi.")
    print("G = pi  =>  pi * exp(-10*r2) = pi  =>  exp(-10*r2) = 1  =>  -10*r2 = 0  =>  r2 = 0.")
    print("r2 = sqrt((x*x+y*y-0.5)^2 + z^2) = 0.")
    print("This holds if and only if z = 0 and x*x+y*y-0.5 = 0.")
    print("So, the preimage of the South Pole is a circle, C1.")
    R_sq = 0.5
    R = np.sqrt(R_sq)
    print(f"Equation of C1: x^2 + y^2 = {R_sq}, z = 0. This is a circle of radius {R:.4f}.")
    print("-" * 50)

    # 5. Find the preimage of p2 (Equatorial Point)
    print(f"Step 2: Find the preimage of p2 = {p2_coords}")
    print("n = (1, 0, 0) requires:")
    print("  nx = sin(G)*cos(f) = 1")
    print("  ny = sin(G)*sin(f) = 0")
    print("  nz = cos(G) = 0")
    print("From nz=0, we get cos(G) = 0, which means G = pi/2 (given the range of G).")
    print("With G=pi/2, sin(G)=1. The equations for nx and ny become:")
    print("  cos(f) = 1  and  sin(f) = 0")
    print("This implies f = 0. From f = atan2(y,x), this means y=0 and x > 0.")
    print("Now we use the condition G = pi/2:")
    print("pi * exp(-10*r2) = pi/2  =>  exp(-10*r2) = 0.5  =>  -10*r2 = ln(0.5) = -ln(2).")
    r2_val = np.log(2) / 10
    r2_val_sq = r2_val**2
    print(f"This gives r2 = ln(2)/10, so r2^2 = (ln(2)/10)^2 = {r2_val_sq:.6f}.")
    print("Substituting the definition of r2, the preimage C2 is a loop defined by:")
    print(f"Equation of C2: (x^2 - {R_sq})^2 + z^2 = {r2_val_sq:.6f}, with y=0 and x>0.")
    print("-" * 50)

    # 6. Calculate the linking number
    print("Step 3: Calculate the linking number Lk(C1, C2).")
    print("We count the signed intersections of loop C2 with the disk D1 spanned by circle C1.")
    print(f"Disk D1 is: x^2 + y^2 <= {R_sq}, z = 0.")
    print("To find where C2 intersects the plane z=0, we set z=0 in the equation for C2:")
    print(f"(x^2 - {R_sq})^2 = {r2_val_sq:.6f}")
    print(f"x^2 - {R_sq} = +/- sqrt({r2_val_sq:.6f}) = +/- {r2_val:.6f}")
    print(f"x^2 = {R_sq} +/- {r2_val:.6f}")
    x_sq_1 = R_sq - r2_val
    x_sq_2 = R_sq + r2_val
    print(f"Two possible values for x^2: {x_sq_1:.6f} and {x_sq_2:.6f}.")
    x1 = np.sqrt(x_sq_1)
    x2 = np.sqrt(x_sq_2)
    print(f"The intersection points with the z=0 plane are at x={x1:.4f} and x={x2:.4f} (since y=0).")
    print("\nNow we check if these points are inside the disk D1 (where x^2 <= 0.5).")
    print(f"Point 1: x^2 = {x_sq_1:.6f}. Since {x_sq_1:.6f} < {R_sq}, this point is INSIDE the disk.")
    print(f"Point 2: x^2 = {x_sq_2:.6f}. Since {x_sq_2:.6f} > {R_sq}, this point is OUTSIDE the disk.")
    print("\nTherefore, the loop C2 pierces the disk D1 exactly once.")
    print("The linking number is the number of net piercings.")
    print("-" * 50)

    # 7. Conclusion
    hopf_charge = 1
    print(f"The Hopf charge is the magnitude of the linking number, which is {hopf_charge}.")
    return hopf_charge

if __name__ == '__main__':
    charge = solve_hopf_charge()
    print(f"\n<<<Hopf Charge = {charge}>>>")
