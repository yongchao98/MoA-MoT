import math

def get_fcc_110_projection_pattern():
    """
    This function calculates the coordinates for a unit cell of an FCC lattice
    projected along the [110] direction, which forms a centered rectangular pattern.

    The projected unit cell is a rectangle with side lengths 'a' and 'a*sqrt(2)'.
    We will assume a lattice constant 'a' = 10 for visualization purposes.
    The pattern consists of atoms at the four corners and one in the center.
    """
    a = 10.0  # Assumed lattice constant for clarity

    # The side lengths of the projected rectangular cell
    side_y = a
    side_x = a * math.sqrt(2)

    # Coordinates of the atoms in the projected 2D unit cell
    corner_bottom_left = (0.0, 0.0)
    corner_bottom_right = (side_x, 0.0)
    corner_top_left = (0.0, side_y)
    corner_top_right = (side_x, side_y)
    center_atom = (side_x / 2, side_y / 2)

    print("The projection of an FCC lattice along the [110] direction results in a centered rectangular pattern.")
    print("This pattern is shown in Image C.")
    print("\nFor a lattice constant 'a', the projected unit cell is a rectangle with side lengths 'a' and 'a*sqrt(2)'.")
    print(f"Assuming a = {a}, the side lengths are {side_y:.2f} and {side_x:.2f}.")
    print("\nThe coordinates of the 5 atoms forming one such centered rectangle are:")
    print(f"Corner 1 (Bottom-Left): ({corner_bottom_left[0]:.2f}, {corner_bottom_left[1]:.2f})")
    print(f"Corner 2 (Bottom-Right): ({corner_bottom_right[0]:.2f}, {corner_bottom_right[1]:.2f})")
    print(f"Corner 3 (Top-Left):    ({corner_top_left[0]:.2f}, {corner_top_left[1]:.2f})")
    print(f"Corner 4 (Top-Right):   ({corner_top_right[0]:.2f}, {corner_top_right[1]:.2f})")
    print(f"Center Atom:            ({center_atom[0]:.2f}, {center_atom[1]:.2f})")
    print("\nThis arrangement matches the topology of the atoms in Image C.")

get_fcc_110_projection_pattern()
