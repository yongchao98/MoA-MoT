import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of triangle T(t).
    """

    # --- Problem Parameters ---
    # Radius of the circumscribed circle
    R = 10
    # Speed of the triangle's vertices along the hexagon's sides
    v = 1

    # --- Derivation Steps ---

    # 1. The side length of a regular hexagon is equal to the radius of its circumscribed circle.
    s = R

    # 2. At t=0, the triangle's vertices are the midpoints of three alternating sides.
    # This forms an equilateral triangle with a side length of 1.5 * R.
    L0 = 1.5 * R

    # 3. The square of the initial side length L0.
    L0_squared = L0**2

    # 4. As the vertices move a distance d = v*t = t from the midpoints,
    # the triangle remains equilateral. The square of the new side length, L(t)^2,
    # is given by L(t)^2 = L0^2 + 3*(v*t)^2.
    # The coefficient for the t^2 term is 3 * v^2.
    coeff_t_squared = 3 * v**2

    # 5. The area of an equilateral triangle is (sqrt(3)/4) * L^2.
    # We can now write the final function for the area A(t).

    # --- Output the Final Formula ---
    print("The area of the triangle T(t) as a function of time t is derived below.")
    print("Note: The rotation of the hexagon does not affect the area of the triangle.")
    print("\nThe area A(t) is given by the formula for an equilateral triangle: A(t) = (sqrt(3)/4) * L(t)^2")
    print("where L(t) is the side length of the triangle at time t.")
    print("\nThe square of the side length, L(t)^2, is:")
    print(f"L(t)^2 = {int(L0_squared)} + {int(coeff_t_squared)}*t^2")
    print("\nSubstituting this into the area formula, we get the final answer:")
    print(f"A(t) = (sqrt(3)/4) * ({int(L0_squared)} + {int(coeff_t_squared)}*t^2)")

solve_triangle_area()