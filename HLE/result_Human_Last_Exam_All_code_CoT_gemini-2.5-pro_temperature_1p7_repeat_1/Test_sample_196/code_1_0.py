import sympy
from sympy import pi, Symbol, Rational, integrate, sqrt

def calculate_volume():
    """
    Calculates the volume of the space enclosed by a cone S1 and an ellipsoid S2 for y > 0.
    """
    y = Symbol('y')

    # --- Define ellipsoid and cone parameters ---
    a_sq, b_sq = 3, 4
    v_y = 4

    # --- Define cross-sectional radii squared as functions of y ---
    # Ellipsoid: x^2/a^2 + y^2/b^2 + z^2/a^2 = 1 => r^2 = a^2 * (1 - y^2/b^2)
    r_ellipsoid_sq = a_sq * (1 - y**2 / b_sq)

    # Cone: Tangent at y=1, where r_ellipsoid^2 = 3*(1-1/4) = 9/4, so r=3/2.
    # The radius is linear with y, r(y=1)=3/2, r(y=4)=0.
    # This gives r(y) = k*(4-y). 3/2 = k*(4-1) => k=1/2.
    # So r_cone(y) = (4-y)/2.
    r_cone_sq = (v_y - y)**2 / 4

    # --- Set up the integral for the volume ---
    # The volume is the integral of the difference in cross-sectional areas.
    # A(y) = pi * (r_cone^2 - r_ellipsoid^2)
    # The ellipsoid is defined for y in [-2, 2]. With y > 0, the range is [0, 2].
    integrand = pi * (r_cone_sq - r_ellipsoid_sq)
    integrand_simplified = sympy.simplify(integrand)

    y_lower = 0
    y_upper = 2

    # --- Perform the integration ---
    volume = integrate(integrand, (y, y_lower, y_upper))
    
    # --- Print the step-by-step calculation ---
    print("Step 1: Define the volume as an integral of the difference in cross-sectional areas.")
    print(f"Volume V = integral from y={y_lower} to y={y_upper} of [A_cone(y) - A_ellipsoid(y)] dy")
    print(f"A_cone(y) = pi * r_cone^2 = pi * (({v_y}-y)/2)^2")
    print(f"A_ellipsoid(y) = pi * r_ellipsoid^2 = pi * {a_sq}*(1 - y^2/{b_sq})")
    print("\nStep 2: Substitute the areas and simplify the integrand.")
    print(f"V = integral from {y_lower} to {y_upper} of ({integrand}) dy")
    print(f"Simplifying the integrand gives: {integrand_simplified}")
    print("\nStep 3: Evaluate the definite integral.")
    # The integral is pi * Integral[(y-1)^2] dy
    # Antiderivative is pi * (y-1)^3 / 3
    antiderivative_str = f"pi * (y-1)^3 / 3"
    print(f"The integral of {integrand_simplified} dy is {antiderivative_str}")
    val_upper = f"pi * (({y_upper}-1)^3) / 3"
    val_lower = f"pi * (({y_lower}-1)^3) / 3"
    print(f"V = [ {antiderivative_str} ] from {y_lower} to {y_upper}")
    print(f"V = ( {val_upper} ) - ( {val_lower} )")
    val_upper_calc = pi * (y_upper - 1)**3 / 3
    val_lower_calc = pi * (y_lower - 1)**3 / 3
    print(f"V = ( {val_upper_calc} ) - ( {val_lower_calc} )")
    
    print("\nStep 4: Final calculation.")
    final_eq = f"V = {val_upper_calc - val_lower_calc}"
    print(final_eq)

calculate_volume()