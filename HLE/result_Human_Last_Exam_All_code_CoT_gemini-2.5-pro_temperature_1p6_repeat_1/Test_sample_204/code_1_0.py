import math

def solve_hopf_charge():
    """
    This script calculates the Hopf charge of the given vector field
    by computing the linking number of the preimages of two points.
    """

    # Step 1: Explain the Hopf charge and the method
    print("The Hopf charge is a topological invariant of a vector field.")
    print("It can be computed as the linking number of the preimages of two distinct points on the target sphere S^2.")
    print("The vector field n = (nx, ny, nz) defines a map from R^3 to S^2, since nx^2 + ny^2 + nz^2 = 1.")
    print("We will calculate the Hopf charge H by finding the linking number Lk(C1, C2), where C1 and C2 are the preimages of two chosen points.")
    print("-" * 30)

    # Step 2: Choose two points on the sphere
    print("Step 1: Choose two points on the target sphere S^2.")
    print("Let p1 = (0, 0, -1) (the South Pole).")
    print("Let p2 = (1, 0, 0) (a point on the equator).")
    print("-" * 30)

    # Step 3: Find the preimage C1 of p1
    print("Step 2: Find the preimage C1 = n^-1(p1).")
    print("The condition n = (0, 0, -1) implies nz = cos(G) = -1.")
    print("This means G must be equal to PI.")
    print("From the definition G = PI * exp(-10 * r2), we get:")
    print("PI = PI * exp(-10 * r2)  =>  exp(-10 * r2) = 1  =>  r2 = 0.")
    print("The definition of r2 is sqrt((x*x + y*y - 0.5)^2 + z*z).")
    print("r2 = 0 implies (x*x + y*y - 0.5)^2 = 0 and z*z = 0.")
    print("This gives the equations: z = 0 and x*x + y*y = 0.5.")
    print("So, the preimage C1 is a circle in the xy-plane with radius sqrt(0.5) centered at the origin.")
    print("-" * 30)

    # Step 4: Find the preimage C2 of p2
    print("Step 3: Find the preimage C2 = n^-1(p2).")
    print("The condition n = (1, 0, 0) implies:")
    print("nx = sin(G)*cos(f) = 1")
    print("ny = sin(G)*sin(f) = 0")
    print("nz = cos(G) = 0")
    print("From nz = cos(G) = 0, we get G = PI / 2.")
    print("From G = PI * exp(-10 * r2), we solve for r2:")
    print("PI/2 = PI * exp(-10 * r2)  =>  exp(-10 * r2) = 0.5  => -10 * r2 = -ln(2)")
    R = math.log(2) / 10
    print(f"This gives r2 = ln(2) / 10, which is approximately {R:.4f}.")
    print("From ny = 0, since sin(G)=sin(PI/2)=1, we need sin(f) = 0.")
    print("From nx = 1, since sin(G)=1, we need cos(f) = 1.")
    print("sin(f)=0 and cos(f)=1 implies f = 0.")
    print("f = atan2(y, x) = 0 implies y = 0 and x > 0.")
    print("Finally, substituting these conditions into the equation for r2^2:")
    print(f"(x*x + y*y - 0.5)^2 + z*z = r2^2 becomes (x^2 - 0.5)^2 + z^2 = (ln(2)/10)^2.")
    print("So, C2 is a closed loop in the half-plane defined by y=0, x>0.")
    print("-" * 30)

    # Step 5: Calculate the linking number
    print("Step 4: Calculate the linking number Lk(C1, C2).")
    print("We can do this by counting the signed intersections of C2 with a disk D1 bounded by C1.")
    print("Let D1 be the disk x^2 + y^2 <= 0.5 in the z=0 plane.")
    print("We orient C1 counter-clockwise (viewed from z>0), so the normal vector to D1 is n_disk = (0, 0, 1).")
    print("To find where C2 intersects the z=0 plane, we set z=0 in its equation:")
    print(f"(x^2 - 0.5)^2 = {R**2:.4f}  =>  x^2 - 0.5 = +/- {R:.4f}.")
    x_sq_1 = 0.5 - R
    print(f"Since x > 0, C2 intersects the x-axis at x=sqrt(0.5 - {R:.4f}) and x=sqrt(0.5 + {R:.4f}).")
    print(f"One point has x^2 = {x_sq_1:.4f} < 0.5, so it is inside the disk D1. This is our single intersection point P.")
    print("To find the sign of the intersection, we need the tangent vector to C2 at P.")
    print("A parameterization of C2 is x(t) = sqrt(0.5 + R*cos(t)), z(t) = R*sin(t) (where R=ln(2)/10).")
    print("The intersection point P corresponds to the parameter t = PI.")
    print("The tangent vector v = (dx/dt, dy/dt, dz/dt) at t=PI is computed as:")
    print(f"v_z = R*cos(PI) = -R = -{R:.4f}")
    print("v_x and v_y are 0.")
    print(f"So the tangent vector v is (0, 0, -{R:.4f}).")
    sign = -1  # from sign(v_z * n_disk_z)
    print(f"The sign of the intersection is sign(v . n_disk) = sign(-R * 1) = {sign}.")
    print("Since there is only one intersection point, the linking number is equal to this sign.")
    print("-" * 30)
    
    # Step 6: Final Answer
    hopf_charge = -1
    print("Final result:")
    print("The Hopf charge is the linking number Lk(C1, C2).")
    print(f"Hopf Charge = {hopf_charge}")
    print(f"\nFinal Equation: Hopf Charge = {hopf_charge}")

if __name__ == '__main__':
    solve_hopf_charge()

<<<-1>>>