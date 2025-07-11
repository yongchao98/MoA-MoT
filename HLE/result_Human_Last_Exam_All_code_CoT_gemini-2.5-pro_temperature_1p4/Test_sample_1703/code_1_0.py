import math

def find_triangle_area_function():
    """
    This function derives and prints the formula for the area of triangle T(t)
    as a function of time t.
    """

    # --- Step 1: Define Given Parameters ---
    # Radius of the circumscribed circle C
    R = 10
    # Speed of the triangle's vertices along the hexagon's sides
    v = 1
    # Note: The angular velocity ω is irrelevant for calculating the area,
    # as area is invariant under rotation.

    # --- Step 2: Determine Geometric Properties ---
    # The side length 's' of a regular hexagon is equal to the radius 'R'
    # of its circumscribed circle.
    s = R

    # --- Step 3: Analyze the Triangle at t=0 ---
    # At t=0, the triangle T(0) connects the midpoints of three alternating
    # sides of the hexagon. This triangle is equilateral due to symmetry.
    # By placing the hexagon's center at the origin, we can calculate the
    # square of the side length of T(0).
    # Let the side length be a_0. It can be shown that a_0 = 1.5 * s.
    # Therefore, the square of the side length at t=0 is:
    a_sq_0 = (1.5 * s)**2

    # --- Step 4: Derive the Side Length Function a_sq(t) ---
    # As the vertices move symmetrically with speed v, the triangle T(t)
    # remains equilateral. The square of its side length, a_sq(t), can be
    # expressed as a function of time t. The derived formula is:
    # a_sq(t) = 3 * (v*t)^2 + a_sq_0
    
    # We define the coefficients for this quadratic equation.
    coeff_t2_side_sq = 3 * v**2
    constant_side_sq = a_sq_0

    # --- Step 5: Derive the Area Function Area(t) ---
    # The area of an equilateral triangle with side 'a' is (sqrt(3)/4) * a^2.
    # Therefore, the area of T(t) is (sqrt(3)/4) * a_sq(t).
    # Area(t) = (sqrt(3)/4) * (3*t^2 + 225)

    # --- Step 6: Print the Final Formula ---
    # The final equation for the area of the triangle T(t) is presented below,
    # showing each numerical component as requested.
    
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print("Area(t) = (√3 / 4) * (a_sq(t))")
    print("\nFirst, the square of the side length of the triangle, a_sq(t), is:")
    print(f"a_sq(t) = {int(coeff_t2_side_sq)} * t^2 + {int(constant_side_sq)}")
    print("\nSubstituting this into the area formula gives the final result:")
    print(f"Area(t) = (√3 / 4) * ({int(coeff_t2_side_sq)} * t^2 + {int(constant_side_sq)})")

# Execute the function to get the solution.
find_triangle_area_function()