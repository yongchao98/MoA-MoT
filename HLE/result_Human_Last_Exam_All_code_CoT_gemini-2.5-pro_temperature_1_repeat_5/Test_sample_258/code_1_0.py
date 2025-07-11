import math

def solve_circle_grid_problem():
    """
    Calculates the minimal and maximal number of grid cells a circle can cross.

    The number of cells a simple closed curve (like a circle) crosses on a grid,
    without passing through any grid vertices, is equal to the number of times
    it intersects the grid lines.

    Let R be the radius of the circle, and (xc, yc) be its center.
    The problem states R = 500.
    The number of unique vertical grid lines the circle crosses is Nv.
    Nv = floor(xc + R) - floor(xc - R).
    The number of unique horizontal grid lines the circle crosses is Nh.
    Nh = floor(yc + R) - floor(yc - R).

    Since the circle is a closed loop, it crosses each of these lines twice.
    Total number of crossings (and thus cells crossed) = 2 * Nv + 2 * Nh.

    Given R=500 (an integer) and the non-tangency condition (which implies
    xc and yc are not integers), Nv and Nh are constant.
    Let xc = I + f, where I is integer and 0 < f < 1.
    Nv = floor(I + f + R) - floor(I + f - R)
       = (I + R + floor(f)) - (I + floor(f - R))
       = (I + R + 0) - (I - R)  (since floor(f-R) = -R for 0<f<1)
       = 2 * R

    So, Nv = 2 * R and Nh = 2 * R.
    Total cells = 2 * (2 * R) + 2 * (2 * R) = 8 * R.
    """
    R = 500

    # Number of unique vertical lines crossed
    # For any non-integer xc, this evaluates to 2*R
    Nv = 2 * R

    # Number of unique horizontal lines crossed
    # For any non-integer yc, this evaluates to 2*R
    Nh = 2 * R

    # Each line is crossed twice by the circle.
    # Total cells crossed = Total number of line crossings
    # For a circle (a simple closed curve not passing through vertices),
    # the number of cells crossed is equal to the number of line crossings.
    # This number is 2*Nv + 2*Nh, which simplifies to 8R.
    # Some interpretations lead to N_cells = Nv + Nh (for monotonic curves)
    # or other variations. The established result for a circle with integer
    # radius R, centered away from vertices and tangents, is 8R.
    num_cells = 8 * R

    # Since this value is constant for any valid circle placement,
    # the minimum and maximum are the same.
    min_cells = num_cells
    max_cells = num_cells

    print(f"The radius of the circle is R = {R}.")
    print(f"The number of unique vertical lines crossed is Nv = 2 * R = {Nv}.")
    print(f"The number of unique horizontal lines crossed is Nh = 2 * R = {Nh}.")
    print("The total number of cells crossed is 2*Nv + 2*Nh = 8*R.")
    print(f"8 * {R} = {num_cells}.")
    print(f"Since this number is constant for any valid circle, the minimal and maximal numbers are the same.")
    print(f"Final Answer: ({min_cells}, {max_cells})")

solve_circle_grid_problem()
<<< (4000, 4000) >>>