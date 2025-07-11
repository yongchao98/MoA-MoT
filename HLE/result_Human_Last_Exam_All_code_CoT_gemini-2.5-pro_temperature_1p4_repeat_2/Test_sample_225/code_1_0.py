def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for a Lissajous curve
    of the form x(t) = cos(nt), y(t) = sin(mt) where n and m are odd and coprime.
    """
    # For the curve (cos(9t), sin(5t)), n=9 and m=5.
    n = 9
    m = 5

    # Check properties (optional for the code, but part of the reasoning)
    # Both n and m are odd, and gcd(9, 5) = 1.

    # Formula for nodes on the y-axis
    on_y_axis = n - 1
    # Formula for nodes on the x-axis
    on_x_axis = m - 1
    # Formula for nodes not on the axes
    off_axis = (n - 1) * (m - 1)

    # Total number of self-intersection points
    total_intersections = on_y_axis + on_x_axis + off_axis

    print(f"The curve is defined by x(t) = cos({n}t), y(t) = sin({m}t).")
    print("Since n and m are odd and coprime, we can calculate the number of self-intersection points.")
    print(f"Number of intersection points on the y-axis: {n} - 1 = {on_y_axis}")
    print(f"Number of intersection points on the x-axis: {m} - 1 = {on_x_axis}")
    print(f"Number of intersection points not on the axes: ({n} - 1) * ({m} - 1) = {off_axis}")
    print("\nThe total number of self-intersection points is the sum of these values.")
    print(f"Total = {on_y_axis} + {on_x_axis} + {off_axis} = {total_intersections}")

solve_lissajous_intersections()