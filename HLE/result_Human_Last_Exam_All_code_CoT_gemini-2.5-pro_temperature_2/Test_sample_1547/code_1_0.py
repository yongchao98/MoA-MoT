import math

def illustrate_hitting_probability():
    """
    This function illustrates that for a simple transient set A in Z^3,
    the hitting probability P_x(tau_A < infinity) tends to 0 as |x| increases,
    which contradicts the condition in the problem statement.

    We consider A = {(0,0,0)} and calculate the hitting probability from
    points x = (k,0,0) for increasing k.
    
    The hitting probability is given by P_x(tau_A < inf) = G(x,0)/G(0,0),
    where G is the Green's function. By translation invariance, this is G(0,-x)/G(0,0).
    
    We use the following known values/approximations for Z^3:
    - G(0,0) is a constant, approximately 1.516386.
    - For large |z|, G(0,z) is approximately (3 / (2 * pi * |z|)).
    """

    # Known value for the Green's function at the origin for SRW in Z^3.
    G_0_0 = 1.516386

    # The point set A is just the origin { (0,0,0) }.
    a = (0, 0, 0)
    
    print("Calculating hitting probability P_x(tau_A < infinity) for A = {(0,0,0)} in Z^3.")
    print("The starting points are x = (k, 0, 0) for k = 1, 2, ..., 10.")
    print("-" * 50)
    
    for k in range(1, 11):
        # Starting point x
        x = (k, 0, 0)
        
        # Vector from x to a is z
        z = (a[0] - x[0], a[1] - x[1], a[2] - x[2])
        
        # Euclidean distance of z
        dist_z = math.sqrt(z[0]**2 + z[1]**2 + z[2]**2)

        # Approximate Green's function value G(0, z) for |z| > 0
        if dist_z == 0:
            # This case corresponds to x in A, prob is 1.
            prob = 1.0
            G_0_z_approx = G_0_0
        else:
            G_0_z_approx = (3 / (2 * math.pi * dist_z))
            # The final probability calculation
            prob = G_0_z_approx / G_0_0

        # Output the numbers used in the final equation
        # Final Equation: P_x(hit A) = G(0, a-x) / G(0,0)
        print(f"For x = {x}:")
        print(f"  P(hit A) = G(0, {z}) / G(0, (0,0,0))")
        print(f"           ≈ {G_0_z_approx:.6f} / {G_0_0:.6f}")
        print(f"           = {prob:.6f}\n")
        
    print("As k increases, the distance |x| increases, and the probability clearly tends to 0, not 1.")
    print("This shows that for a finite (and thus transient) set, it's impossible to have")
    print("P_x(tau_A < infinity) = 1 for infinitely many x.")

if __name__ == '__main__':
    illustrate_hitting_probability()
