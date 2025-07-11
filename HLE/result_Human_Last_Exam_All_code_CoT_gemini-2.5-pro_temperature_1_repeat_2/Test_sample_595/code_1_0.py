import math

def calculate_squares_crossed():
    """
    Calculates the maximum number of grid squares the triangle's perimeter can cross.
    This is based on the optimal orientation where the triangle's legs are at 45 degrees
    to the coordinate axes.
    """

    # Side lengths of the right-angled isosceles triangle
    leg_length = 18
    hypotenuse_length = 18 * math.sqrt(2)

    # In the optimal orientation, the leg vectors are rotated by 45 degrees.
    # The change in x and y for the legs is leg_length / sqrt(2).
    delta_leg = leg_length / math.sqrt(2)  # This is 9 * sqrt(2)

    # Define vertex coordinates with a small offset to avoid lattice points.
    # A is the right-angle vertex.
    # We choose offsets delta_x=1/3, delta_y=1/4 to ensure no line passes through lattice points.
    ax, ay = 1/3, 1/4
    # B and C are the other two vertices.
    bx = ax + delta_leg
    by = ay + delta_leg
    cx = ax - delta_leg
    cy = ay + delta_leg

    # --- Calculate squares for side AB ---
    # Number of vertical lines crossed (integers between ax and bx)
    nv_ab = math.floor(bx) - math.ceil(ax) + 1
    # Number of horizontal lines crossed (integers between ay and by)
    nh_ab = math.floor(by) - math.ceil(ay) + 1
    # Number of squares crossed by segment AB
    k_ab = 1 + nv_ab + nh_ab

    # --- Calculate squares for side AC ---
    # Number of vertical lines crossed (integers between cx and ax)
    nv_ac = math.floor(ax) - math.ceil(cx) + 1
    # Number of horizontal lines crossed (integers between ay and cy)
    nh_ac = math.floor(cy) - math.ceil(ay) + 1
    # Number of squares crossed by segment AC
    k_ac = 1 + nv_ac + nh_ac

    # --- Calculate squares for side BC (hypotenuse) ---
    # Number of vertical lines crossed (integers between cx and bx)
    nv_bc = math.floor(bx) - math.ceil(cx) + 1
    # Number of horizontal lines crossed (y is constant, so 0)
    nh_bc = 0
    # Number of squares crossed by segment BC
    k_bc = 1 + nv_bc + nh_bc

    # The squares containing the three vertices are counted twice. We subtract 3 for the overlaps.
    overlaps = 3
    total_squares = k_ab + k_ac + k_bc - overlaps

    print("Calculation for the optimal orientation (legs at 45 degrees):")
    print(f"Side AB (leg): Crosses {nv_ab} vertical and {nh_ab} horizontal lines. Total squares = 1 + {nv_ab} + {nh_ab} = {k_ab}")
    print(f"Side AC (leg): Crosses {nv_ac} vertical and {nh_ac} horizontal lines. Total squares = 1 + {nv_ac} + {nh_ac} = {k_ac}")
    print(f"Side BC (hypotenuse): Crosses {nv_bc} vertical and {nh_bc} horizontal lines. Total squares = 1 + {nv_bc} + {nh_bc} = {k_bc}")
    print(f"\nTotal squares for each side summed = {k_ab} + {k_ac} + {k_bc} = {k_ab + k_ac + k_bc}")
    print(f"Subtracting {overlaps} for the squares at the vertices which are counted twice.")
    print(f"Final number of squares k = {k_ab + k_ac + k_bc} - {overlaps} = {total_squares}")

calculate_squares_crossed()