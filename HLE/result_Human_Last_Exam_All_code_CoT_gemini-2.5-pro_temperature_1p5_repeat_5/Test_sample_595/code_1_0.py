import math

def solve_triangle_crossing():
    """
    Calculates the maximum number of grid squares a specific triangle's perimeter can cross.

    The triangle has side lengths 18, 18, and 18*sqrt(2), making it an
    isosceles right-angled triangle.

    The number of squares crossed (k) is the sum of horizontal and vertical grid
    lines crossed. This is maximized when the triangle is oriented such that its
    legs have slopes of 2 and -1/2 with respect to the coordinate axes.

    We place the triangle's vertices (A, B, C) in the coordinate plane based on this
    optimal orientation and calculate the number of crossings. A is the right-angle vertex.
    """
    L = 18

    # The optimal orientation for maximizing grid crossings corresponds to
    # a leg vector v1 making an angle theta with the x-axis where tan(theta) = 2.
    # v1 = (x1, y1) such that x1^2 + y1^2 = L^2 and y1/x1 = 2.
    # This gives x1 = L/sqrt(5) and y1 = 2L/sqrt(5).
    x1 = L / math.sqrt(5)
    y1 = 2 * L / math.sqrt(5)

    # The second leg vector, v2, is perpendicular to v1.
    # v2 = (-y1, x1)
    x2 = -y1
    y2 = x1

    # We place the right-angle vertex A such that all vertex coordinates are positive.
    # A shift by (y1, 0) is sufficient. A tiny offset epsilon is assumed to avoid
    # vertices landing on grid lines, as per the problem.
    # For calculation purposes, we can ignore epsilon as it does not affect the floor value.
    v_a_x = y1
    v_a_y = 0

    # The other two vertices, B and C, are found by adding the leg vectors to A.
    # Perimeter is A -> B -> C -> A
    v_b_x = v_a_x + x1
    v_b_y = v_a_y + y1

    v_c_x = v_a_x + x2
    v_c_y = v_a_y + y2

    # Get the floor of each coordinate, which represents the integer grid cell index.
    f_ax = math.floor(v_a_x)  # floor(36/sqrt(5)) = floor(16.099) = 16
    f_ay = math.floor(v_a_y)  # floor(0) = 0

    f_bx = math.floor(v_b_x)  # floor(54/sqrt(5)) = floor(24.149) = 24
    f_by = math.floor(v_b_y)  # floor(36/sqrt(5)) = floor(16.099) = 16

    f_cx = math.floor(v_c_x)  # floor(0) = 0
    f_cy = math.floor(v_c_y)  # floor(18/sqrt(5)) = floor(8.049) = 8

    # Calculate the number of vertical grid lines crossed for each segment.
    nx_ab = abs(f_bx - f_ax)
    nx_bc = abs(f_cx - f_bx)
    nx_ca = abs(f_ax - f_cx)
    total_nx = nx_ab + nx_bc + nx_ca

    # Calculate the number of horizontal grid lines crossed for each segment.
    ny_ab = abs(f_by - f_ay)
    ny_bc = abs(f_cy - f_by)
    ny_ca = abs(f_ay - f_cy)
    total_ny = ny_ab + ny_bc + ny_ca

    k = total_nx + total_ny

    print(f"The number of vertical lines crossed is {nx_ab} + {nx_bc} + {nx_ca} = {total_nx}.")
    print(f"The number of horizontal lines crossed is {ny_ab} + {ny_bc} + {ny_ca} = {total_ny}.")
    print(f"The largest number k such that the triangle's perimeter can pass through at least k coordinate grid squares is:")
    print(f"k = {total_nx} + {total_ny} = {k}")

solve_triangle_crossing()
<<<80>>>